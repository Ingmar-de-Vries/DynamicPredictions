function DynamicPredictions_STATS_ERFdynamicRSA_ROIsource(cfg,~,iroi)

restoredefaultpath
% You can run this function locally (by looping over subjects), or send to a cluster (which I did using the accompanying "cluser_shell.m" function)
% set some directories depending on running this locally or on a cluster:
if isfield(cfg,'cluster')
    addpath('//XXX/ActionPrediction/toolboxes/fieldtrip-20191113');
    addpath(genpath('//XXX/ActionPrediction/code'));
    addpath(genpath('//XXX/matlab_toolboxes/CoSMoMVPA-master'));
    cfg.path = '//XXX/ActionPrediction'; 
else
    addpath('\\XXX\ActionPrediction\toolboxes\fieldtrip-20191113');
    addpath(genpath('\\XXX\ActionPrediction\code'));
    addpath(genpath('\\XXX\matlab_toolboxes\CoSMoMVPA-master'));
    cfg.path = '\\XXX\ActionPrediction';
end
ft_defaults

if cfg.similarity == 0
    corrORglm = 'corr';
elseif cfg.similarity == 1
    corrORglm = ['pcr_' num2str(cfg.nPCAcomps) 'comps'];
end

%% input and output folders
indir = fullfile(cfg.path,'data','MEG',['source_' cfg.atlas],'RSA', [corrORglm '_' num2str(cfg.randshuff(1)) 'iterations_' num2str(ceil(cfg.randshuff(2)*1000)) 'msec']);
outdir = fullfile(cfg.path,'data','MEG',['source_' cfg.atlas],'RSA','statistics', [corrORglm '_' num2str(cfg.randshuff(1)) 'iterations_' num2str(ceil(cfg.randshuff(2)*1000)) 'msec']);
if ~exist(outdir,'dir')
    mkdir(outdir);
end

% if on cluster, only send a single ROI as job to cluster, as indicated by iroi in cluster_shell.m script
if isfield(cfg,'cluster')
    cfg.ROIVec = iroi;
end

% divide ROIs
load(fullfile(cfg.path,'code','neuralDecoding','ROIdefinitions'),'ROIdefinition');

% load dRSA curves resulting from PCR on simulated data to compute the representational spread (RS) index (see methods in article for more information)
load(fullfile(indir,'..','modelautocorr_slopes'),'modelautocorr_slopes');

% remove absolute and relative acceleration, as we won't compute the sustained representation index for those
modelautocorr_slopes(8:9,:) = [];

for iROI = cfg.ROIVec
    
    %% load data, combine and average subjects
    dRSAallperSub = zeros(cfg.subnum,(5-cfg.randshuff(2))*cfg.downsample,(5-cfg.randshuff(2))*cfg.downsample,length(cfg.models2test));
    subcount = 0;
    for isub = cfg.sub4stat
        subcount = subcount+1;

        % load data
        fn = sprintf('%s%cdRSA_SUB%02d_%dHz_%s_smMEG%d_smRDMneu%d_smRDMmod%d',indir, filesep, isub, cfg.downsample, ROIdefinition.names{iROI}{:}, cfg.smoothingMSec, cfg.smoothNeuralRDM, cfg.smoothModelRDM);
        fprintf('loading %s..\n',fn);        
        load(fn,'dRSAall'); 

        ROIname = ROIdefinition.names{iROI}{:};

        % fisher Z transform
        if cfg.fisherz
            dRSAallperSub(subcount,:,:,:) = atanh(dRSAall);
        else
            dRSAallperSub(subcount,:,:,:) = dRSAall;
        end
            
    end% subject loop
    
    % some parameters for plotting:
    maxLag = 1;% maximum lag to plot
    TimeVec = 0:1/cfg.downsample:size(dRSAallperSub,2)/cfg.downsample-1/cfg.downsample;% dRSA time vector
    TimeVec2 = -maxLag:1/cfg.downsample:maxLag;% time vector for dRSA lag plot (i.e., [-maxLag maxLag])
    
    %slice neural vectors per model time point (i.e., locked to on-diagonal / zero lag in 2D dRSA plot), and then stack
    tRange=maxLag*cfg.downsample;% maximum lag in samples
    
    rstack = zeros(cfg.subnum,size(dRSAallperSub,4),length(TimeVec),length(TimeVec2));
    for isub = 1:cfg.subnum
        for iRDM = 1:size(dRSAallperSub,4)
            for iModelTime = 1:length(TimeVec)
                
                timeidx = iModelTime - tRange:iModelTime + tRange;% time indices for current slice
                
                % for early time points, negative lag doesn't fall in video time, while for later time points, positive lag doesn't fall in video time
                % indices that fall before or after video time will be removed, but first give them value of 1 here, and NaN later below
                NotInVid = logical((timeidx < 1)+(timeidx > length(TimeVec)));
                timeidx(NotInVid) = 1;
                
                % extract slice from 2D dRSA matrix
                slice = squeeze(dRSAallperSub(isub,iModelTime,timeidx,iRDM));
                
                % now remove indices that fall before or after video
                slice(NotInVid) = NaN;
                
                % and stack slices, which are now locked to diagonal / zero lag in 2D dRSA matrix
                rstack(isub,iRDM,iModelTime,:) = slice;
                
            end
        end
    end
    
    % Average over all slides (i.e., video time) to get single dRSA curve for lag plot
    rstackLine = squeeze(nanmean(rstack,3));
    
    %% determine peak latency and representational spread (RS) index, using jackknife approach
    % determine index for zero lag
    newzeroID = dsearchn(TimeVec2',0);
    
    % determining the peak latency and RS are posthoc analyses that depend on observing an actual peak in the dRSA curves
    % in the current dataset, for the third model (optical flow vector direction) we observed both a predictive and a lagged peak
    % the following variable sets for which models to compute peak latency and RS separately for the predictive and lagged peaks
    RDMwith2peaks = 3;
    
    % initialize
    peakLatency = zeros(cfg.subnum,size(rstackLine,2),2);% if single peak, its in: peakLatency(:,:,2), if two peaks, predictive is in peakLatency(:,:,1) and lagged is (:,:,2)
    RS = zeros(cfg.subnum,size(rstackLine,2),2,length(TimeVec2));% RS = representational spread index, same structure as peakLatency    
    for iRDM = 1:size(rstackLine,2)
        
        % look for only single peak
        if ~any(RDMwith2peaks == iRDM)
            intervals = [0 0 ; 1 length(TimeVec2)];                       
            dirs = 2;% if dirs = 2, it will be single peak in whole [-1 1] lag plot interval
            
        % look for two peaks: a predictive and a lagged peak
        else
            intervals = [1 newzeroID ; newzeroID-1 length(TimeVec2)];
            dirs = 1:2;% if dirs = 1:2, 1 = predictive peak, 2 = lagged peak            
        end
        
        for isub = 1:cfg.subnum
            
            for idir = dirs
                
                if ~cfg.jackknife% single-participant approach for estimating peak latency and RS
                    idx = isub;
                    curve = squeeze(rstackLine(idx,iRDM,:));
                elseif cfg.jackknife% jackknife approach for estimating peak latency and RS
                    idx = 1:size(rstackLine,1) ~= isub;% leave-one-out
                    curve = squeeze(mean(rstackLine(idx,iRDM,:)));
                end
                
                % find peaks
                [maxVAL,maxID] = findpeaks(curve(intervals(idir,1):intervals(idir,2)));
                
                % if no peak found... select maximum instead
                if isempty(maxVAL)
                    [maxVAL,maxID] = max(curve(intervals(idir,1):intervals(idir,2)));
                end
                
                % if multiple peaks, find highest of those
                [maxVAL,ID] = max(maxVAL);
                maxID = maxID(ID);
                maxID = maxID+intervals(idir,1)-1;% if looking at predictive or lagged peak separately, adjust 
                
                % store jackknifed peak latencies
                peakLatency(isub,iRDM,idir) = TimeVec2(maxID);
                
                %% calculate slope from peak to adjecent negative peaks on either side of peak
                % first scale both data curve and model autocorrelation curve to max = 1
                curve = curve/maxVAL;
                auto = modelautocorr_slopes(iRDM,:);
                auto = auto/max(auto);
                                
                % align curve peak to model autocorrelation peak (i.e., at zero lag), call the rest outside of the interval: NaN
                latency = maxID - newzeroID;
                if latency > 0% if positive (i.e., lagged) latency, shift curve to left
                    curve = [curve(latency+1:end) ; NaN(latency,1)];
                elseif latency < 0% if negative (i.e., predictive) latency, shift curve to right
                    curve = [NaN(-latency,1) ; curve(1:end+latency)]; 
                end% if latency = 0, no realigning needed
                
                % compute representational spread (i.e., relative dRSA value above and beyond what would be expected based on model autocorrelation alone) 
                RS(isub,iRDM,idir,:) = curve - auto';% representational spread index                 
                
            end% idir loop
            
        end% isub loop
        
        % Retrieve individual-subject latencies from subsample jackknife averages (Smulders 2010 Psychophysiology)
        if cfg.jackknife
            n = size(peakLatency,1);
            peakLatency(:,iRDM,:) = repmat(n .* mean(squeeze(peakLatency(:,iRDM,:))),n,1) - (n-1) .* squeeze(peakLatency(:,iRDM,:));
            
            for itime = 1:size(RS,4)
                RS(:,iRDM,:,itime) = repmat(n .* mean(squeeze(RS(:,iRDM,:,itime))),n,1) - (n-1) .* squeeze(RS(:,iRDM,:,itime));
            end
        end
        
    end% iRDM loop  
    
    %% Statistics
    % run cluster-based permutation test using Fieldtrip
    % Here we call a shell script "DynamicPredictions_runFTstats.m" where parameters are set, and where the Fieldtrip functions are called from
    
    % first set some pre-defined parameters and make filename
    pthreshstring = num2str(cfg.pthresh);
    pthreshstring(1:2) = [];
    pthreshcluststring = num2str(cfg.pthreshclust);
    pthreshcluststring(1:2) = [];

    fnSTATS = sprintf('%s%cSTATS_p%s_p%s_fisherz%d_%dHz_%s_smMEG%d_smRDMneu%d_smRDMmod%d',outdir, filesep, pthreshcluststring, pthreshstring, cfg.fisherz, cfg.downsample, ROIname, cfg.smoothingMSec, cfg.smoothNeuralRDM, cfg.smoothModelRDM);        
    
    % compare against a correlation of zero
    chance = 0;
    
    % run stats on 2D dRSA plots (i.e., neural time by model time), or only on dRSA curves below (which is what I did for the article)
    if cfg.timeXtime
        
        %time-frequency (or time-time, as is the case here)
        testtype = 'multivar_timefreq';%univar_timeseries, univar_timefreq, univar_timefreqchan, multivar_timeseries, multivar_timefreq
        signposdRSA = false(size(dRSAallperSub,4),length(TimeVec),length(TimeVec));
        signnegdRSA = false(size(dRSAallperSub,4),length(TimeVec),length(TimeVec));

        for iRDM = 1:size(dRSAallperSub,4)
            addinfo.time = TimeVec;
            addinfo.freq = TimeVec;% for 2D dRSA plot act as if second dimension is also time (instead of freq)
            [~,signposdRSA(iRDM,:,:),signnegdRSA(iRDM,:,:)] = DynamicPredictions_runFTstats(squeeze(dRSAallperSub(:,:,:,iRDM)),addinfo,[],testtype,cfg.pthresh,cfg.pthreshclust,chance);
        end
        signdRSA = logical(signposdRSA+signnegdRSA);
    
    else    
        signdRSA = 'Set cfg.timeXtime = 1 to run stats on 2D dRSA plot';        
    end
    
    % dRSA curve, i.e., after averaging 2D dRSA over stimulus time (after locking to on-diagonal) 
    testtype = 'multivar_timeseries';%univar_timeseries, univar_timefreq, univar_timefreqchan, multivar_timeseries, multivar_timefreq
    signposLine = false(size(rstackLine,2),size(rstackLine,3));
    signnegLine = false(size(rstackLine,2),size(rstackLine,3));
    signposRS = false(size(rstackLine,2),size(rstackLine,3),2);
    signnegRS = false(size(rstackLine,2),size(rstackLine,3),2);
    
    for iRDM = 1:size(dRSAallperSub,4)
        
        addinfo.time = TimeVec2;
        addinfo.tail = 1;% only positive tail
        [~,signposLine(iRDM,:),signnegLine(iRDM,:)] = DynamicPredictions_runFTstats(squeeze(rstackLine(:,iRDM,:)),addinfo,[],testtype,cfg.pthresh,cfg.pthreshclust,chance);
        
        addinfo.tail = 1;% only positive tail
        [~,signposRS(iRDM,:,2),signnegRS(iRDM,:,2)] = DynamicPredictions_runFTstats(squeeze(RS(:,iRDM,2,:)),addinfo,[],testtype,cfg.pthresh,cfg.pthreshclust,chance);
        if iRDM == RDMwith2peaks
            [~,signposRS(iRDM,:,1),signnegRS(iRDM,:,1)] = DynamicPredictions_runFTstats(squeeze(RS(:,iRDM,1,:)),addinfo,[],testtype,cfg.pthresh,cfg.pthreshclust,chance);
        end

    end
    signLine = logical(signposLine+signnegLine);
    signRS = logical(signposRS+signnegRS);
    
    dRSAall = dRSAallperSub;
    save(fnSTATS,'dRSAall','rstack','rstackLine','RS','peakLatency','signLine','signRS','signdRSA','ROIname');
    
end

end