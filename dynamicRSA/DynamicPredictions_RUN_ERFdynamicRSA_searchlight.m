function DynamicPredictions_RUN_ERFdynamicRSA_searchlight(cfg,iSub,~)

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

nstim = cfg.nstim;

%% input and output folders
indir = fullfile(cfg.path,'data','MEG','PreprocessedSensor');
indirModel = fullfile(cfg.path,'data','modelRDMs');
indirEye = fullfile(cfg.path,'data','MEG','ELpreprocessed');
subfilz = dir(fullfile(indir,'preprocessedSensor*'));

if cfg.similarity == 0
    corrORglm = 'corr';
elseif cfg.similarity == 1
    corrORglm = ['pcr_' num2str(cfg.nPCAcomps) 'comps'];
end
outdir = fullfile(cfg.path,'data','MEG',['searchlight_' cfg.atlas],'RSA', [corrORglm '_' num2str(cfg.randshuff(1)) 'iterations_' num2str(ceil(cfg.randshuff(2)*1000)) 'msec']);
if ~exist(outdir,'dir')
    mkdir(outdir);
end

% check if this subject already ran, if so skip
fn2save = sprintf('%s%cdRSA_SUB%02d_%dHz_smMEG%d_smRDMneu%d_smRDMmod%d',outdir, filesep, iSub, cfg.downsample, cfg.smoothingMSec, cfg.smoothNeuralRDM, cfg.smoothModelRDM);

% check if this subject already ran, if so skip
if exist([fn2save '.mat'],'file')
    return
end

% source data
indirSource = fullfile(cfg.path,'data','brainstorm_database','ActionPrediction','data',subfilz(iSub).name(end-7:end-4),[subfilz(iSub).name(end-7:end-4) '_allSensorData']);
indirAtlas = fullfile(cfg.path,'data','brainstorm_database','ActionPrediction','anat',subfilz(iSub).name(end-7:end-4));

%% load and prepare neural data
fn2load = sprintf('%s%c%s',indir, filesep, subfilz(iSub).name);
fprintf('load %s',fn2load);
load(fn2load,'FTdata'); % load in preprocessed data in ft format

% smooth with sliding window
if cfg.smoothingMSec>0
    
    % transform smoothing in msec to smoothing in samples, because that's what ft_preproc_smooth uses
    fsample = round(1/(FTdata.time(2)-FTdata.time(1)));
    smoothingSamples = cfg.smoothingMSec / (1000/fsample);
    
    % now do the smoothing using Fieldtrip function: 
    for i = 1:size(FTdata.trial,1)
        FTdata.trial(i,:,:) = ft_preproc_smooth(squeeze(FTdata.trial(i,:,:)), smoothingSamples);
    end
end

% downsample using Fieldtrip function:
if ~any(cfg.downsample==[0 500])
    ftcfg = [];
    ftcfg.resamplefs = cfg.downsample;
    ftcfg.detrend = 'no';
    FTdata = ft_resampledata(ftcfg, FTdata);
else
    cfg.downsample = round(1/(FTdata.time(2)-FTdata.time(1)));
end

% now convert to cosmoMVPA format:
ds = cosmo_meeg_dataset(FTdata); clear FTdata

% correction trial codes, i.e., for some subjects the usable catch trials still have +100
catchtrials = ds.sa.trialinfo>100;
ds.sa.trialinfo(catchtrials) = ds.sa.trialinfo(catchtrials)-100;

% targets
ds.sa.targets=ds.sa.trialinfo;

% throw away trials without recognizable event / target
cond_idx = find(ismember(ds.sa.targets,1:nstim));
ds=cosmo_slice(ds,cond_idx,1);

% average over all trials beloning to the same unique sequence for correlation-based dissimilarity analysis
ds = cosmo_fx(ds, @(x)mean(x,1),'targets');%average over all trials

% define the measure and arguments.
measure=@cosmo_dissimilarity_matrix_measure;
measure_args.metric='correlation';
measure_args.center_data=true;

% time definitions
fsMod = 100;% sample rate kinematic recordings
tMod = 0:1/fsMod:5-1/fsMod;
tNeural = ds.a.fdim.values{2};
triallengthsec = 5;
newtimelength = (triallengthsec-cfg.randshuff(2))*cfg.downsample;

%% load model RDMs
load([indirModel filesep 'ActionPrediction_dynRDM_graysmooth'],'RDMgraysmooth');
load([indirModel filesep 'ActionPrediction_dynRDM_OFmag'],'RDMoptflow_mag');
load([indirModel filesep 'ActionPrediction_dynRDM_OFdir'],'RDMoptflow_dir');
load([indirModel filesep 'ActionPrediction_dynRDM_posturedep'],'RDMposture_dep');
load([indirModel filesep 'ActionPrediction_dynRDM_postureinvar'],'RDMposture_invar');
load([indirModel filesep 'ActionPrediction_dynRDM_motiondep'],'RDMmotion_dep');
load([indirModel filesep 'ActionPrediction_dynRDM_motioninvar'],'RDMmotion_invar');
load([indirModel filesep 'ActionPrediction_dynRDM_accdep'],'RDMacc_dep');
load([indirModel filesep 'ActionPrediction_dynRDM_accinvar'],'RDMacc_invar');

dynRDM{1} = RDMgraysmooth;
dynRDM{2} = RDMoptflow_mag;
dynRDM{3} = RDMoptflow_dir;
dynRDM{4} = RDMposture_dep;
dynRDM{5} = RDMposture_invar;
dynRDM{6} = RDMmotion_dep;
dynRDM{7} = RDMmotion_invar;
dynRDM{8} = RDMacc_dep;
dynRDM{9} = RDMacc_invar;

% clear up workspace
clear RDMgraysmooth RDMopticalflow_FBmag RDMopticalflow_FB RDMposture_dep RDMposture_invar RDMmotion_dep RDMmotion_invar RDMacc_dep RDMacc_invar

%% load and prepare eyetracker RDMs
fn= sprintf('%s%cSUB%02d_dynRDM_eyeTracker',indirEye, filesep, iSub);
load(fn,'RDMeyePOS');
dynRDM{10} = RDMeyePOS;
clear RDMeyePOS;

%% MNE source reconstruction using sensor data and inversion kernel created in Brainstorm
sourcefile = dir(fullfile(indirSource,'*MEG_GRAD_MEG_MAG_KERNEL_211116*'));
fn2load = sprintf('%s%c%s',indirSource,filesep,sourcefile.name);
kernel = load(fn2load);
chanIdx = kernel.GoodChannel;
kernel = kernel.ImagingKernel;

atlasfile = dir(fullfile(indirAtlas,'*pial_low*'));
fn2load = sprintf('%s%c%s',indirAtlas,filesep,atlasfile(1).name);
atlas = load(fn2load);
if iSub == 12 && contains(cfg.atlas,'400')
    atlasIdx = contains(extractfield(atlas.Atlas,'Name'),'400');
elseif iSub == 12 && contains(cfg.atlas,'200')
    atlasIdx = contains(extractfield(atlas.Atlas,'Name'),'200');
else
    atlasIdx = contains(extractfield(atlas.Atlas,'Name'),cfg.atlas);
end
atlas = atlas.Atlas(atlasIdx).Scouts;

% do some atlas fixing:
% 1) in Schaefer there is also a background, for subject 12 for some reason the atlas has a slightly different name, and so do all the ROIs
% 2) in HCP some of the L and R to indicate hemisphere behind the parcel names are gone
idx2remove = contains(extractfield(atlas,'Label'),'Background');
atlas(idx2remove) = [];

for iparcel = 1:length(atlas)
    
    if ~contains(atlas(iparcel).Label,' L') && ~contains(atlas(iparcel).Label,' R')
        
        if strcmp(atlas(iparcel).Region,'LU')
            atlas(iparcel).Label = [atlas(iparcel).Label ' L'];
        elseif strcmp(atlas(iparcel).Region,'RU')
            atlas(iparcel).Label = [atlas(iparcel).Label ' R'];
        end
        
    end
    
end

targets = ds.sa.targets;
ds=cosmo_map2meeg(ds);% convert data to FT format

% select only gradiometers because MNE is only applied to those
ds.label = ds.label(chanIdx);
ds.trial = ds.trial(:,chanIdx,:);

%% start iteration loop
% create matrix with times for randomly shuffling sequence onsets over X iterations
% shuffleTime = start times size [iterations X stimuli];

% make sure rand numbers are different for each subject, and each time this script is ran
rng((round(sum(clock))+iSub*1000)*100000);

% make sure it's a multiple of 0.01, which is the video sampling rate (50 Hz), such that exact same time point is picked for neural and kinematic data
shuffleTime = round((rand(cfg.randshuff(1),nstim)*(cfg.randshuff(2)-1/cfg.downsample))./0.01).*0.01;

dRSAperIter = zeros(cfg.randshuff(1),length(atlas),size(cfg.peaks2test,1));
for iter = 1:cfg.randshuff(1)
    
    % determine indices for stimulus-specific re-alignment for current iteration
    shuffleIDmodel = dsearchn(tMod',shuffleTime(iter,:)')';
    shuffleIDneural = dsearchn(tNeural',shuffleTime(iter,:)')'; 
    
    %% prepare neural RDM for current iteration
    ds_iter = ds;
    
    temptrial = zeros(nstim,length(ds_iter.label),cfg.downsample*(triallengthsec-cfg.randshuff(2)));
    for istim = 1:nstim
        
        tNewID = shuffleIDneural(istim):shuffleIDneural(istim)+(triallengthsec-cfg.randshuff(2))*cfg.downsample-1;
        
        temptrial(istim,:,:) = squeeze(ds_iter.trial(istim,:,tNewID));
        
    end
    ds_iter.trial = temptrial; clear temptrial
    ds_iter.time = 0:1/cfg.downsample:triallengthsec-cfg.randshuff(2)-1/cfg.downsample;
    ds_iter.trialinfo = 1:nstim;
        
    % now loop over sources for computing neural RDM per source
    neuralRDM = zeros((nstim*nstim-nstim)/2,length(ds_iter.time),length(atlas));
    for isource = 1:length(atlas)
        
        ds_source = ds_iter;
        
        % source reconstruction
        vertices = atlas(isource).Vertices;
        kernelSel = kernel(vertices,:);
        sourcedata = zeros(size(ds_source.trial,1),size(kernelSel,1),size(ds_source.trial,3));
        for itrial = 1:size(ds_source.trial,1)
            temp = squeeze(ds_source.trial(itrial,:,:));
            sourcedata(itrial,:,:) = kernelSel*temp;
        end
        ds_source.trial = sourcedata;
        ds_source.label = vertices';
        
        ds_source = cosmo_meeg_dataset(ds_source);% convert data back to cosmo format
        ds_source.sa.targets = [1:nstim]';% for some reason this field gets lost when coverting back and forth from cosmo to FT format
        
        % define neighborhood over time
        nbrhood=cosmo_interval_neighborhood(ds_source,'time','radius',0);
        
        % compute neural RDM
        ds_res=cosmo_searchlight(ds_source,nbrhood,measure,measure_args);
        
        % unpack the matrices
        MATxTIME = zeros(nstim,nstim,length(ds_res.a.fdim.values{:}));
        for itime = 1:length(ds_res.a.fdim.values{:})
            MATxTIME(:,:,itime) = squareform(ds_res.samples(:,itime));
        end
        
        % smooth neural RDM across time
        if cfg.smoothNeuralRDM > 0
            for i=1:size(MATxTIME,1)
                MATxTIME(i,:,:) = ft_preproc_smooth(squeeze(MATxTIME(i,:,:)),cfg.smoothNeuralRDM);
            end
        end
                
        for iBin = 1:size(MATxTIME,3)
            matrices = squeeze(MATxTIME(:,:,iBin));            
            matrices = matrices / sum(mean(matrices,2)); % normalize, shouldn't do anything to final result
            
            % simplify matrix by making it symmetrical, which it already is anyway...
            DSM = (matrices + matrices')/2;
            
            DSM(logical(eye(size(DSM)))) = 0;
            neuralRDM(:,iBin,isource) = squareform(DSM);
        end
        
    end% source loop
    
    %% prepare model RDM for current iteration
    % re-align model RDMs for current iteration
    modelRDMsquare = zeros(length(dynRDM),nstim,nstim,cfg.downsample*(triallengthsec-cfg.randshuff(2)));
    for istim1 = 1:nstim
        for istim2 = 1:nstim
            
            tNewID1 = shuffleIDmodel(istim1):shuffleIDmodel(istim1)+(triallengthsec-cfg.randshuff(2))*cfg.downsample-1;
            tNewID2 = shuffleIDmodel(istim2):shuffleIDmodel(istim2)+(triallengthsec-cfg.randshuff(2))*cfg.downsample-1;
            
            for iRDM = 1:length(dynRDM)
                
                temptrial = squeeze(dynRDM{iRDM}(istim1,istim2,tNewID1,tNewID2));
                
                temptrial = diag(temptrial);
                
                modelRDMsquare(iRDM,istim1,istim2,:) = temptrial;
                
            end% model loop
        end% second stimulus loop
    end% first stimulus loop
    
    % for relative models take average of (stim1*stim2trans) and (stim1trans*stim2), which is not exactly the same because procrustes is not symmetrical
    for iRDM = [5 7 9]
        for itime = 1:size(modelRDMsquare,4)
            temp = squeeze(modelRDMsquare(iRDM,:,:,itime));
            modelRDMsquare(iRDM,:,:,itime) = (temp + temp')/2;
        end% time loop
    end% model loop
    
    %% some more model RDM preparation
    % extract vector RDM from triangle of square RDM
    modelRDM = zeros(size(modelRDMsquare,1),size(modelRDMsquare,4),(nstim*nstim-nstim)/2);
    for iRDM = 1:size(modelRDMsquare,1)
        for iBin = 1:size(modelRDMsquare,4)
            modelRDM(iRDM,iBin,:) = squareform(tril(squeeze(modelRDMsquare(iRDM,:,:,iBin)),-1));
        end
    end
    
    % smooth model RDM across time
    if cfg.smoothModelRDM
        for i=1:size(modelRDM,3)
            modelRDM(:,:,i) = ft_preproc_smooth(squeeze(modelRDM(:,:,i)),cfg.smoothModelRDM);
        end
    end
      
    %% rescale all RDMs to same [0 2] interval. Unscaled might be problematic for PCA or regression (i.e., larger scale = more variance = higher component)
    % scale RDMs once across all time points, but separately per model and per source to keep internal structure of variance intact
    % we want all models and sources roughly on the same scale for the PCA and least-squares regression below    
    for iRDM = 1:size(modelRDM,1)
        % no need to center individual time points for modelRDM, because matlab's PCA function below does it automatically anyway (except for Xtest,
        % for which we do it below separately before running the regression)
        modelRDM(iRDM,:,:) = reshape(rescale(reshape(squeeze(modelRDM(iRDM,:,:)),size(modelRDM,2)*size(modelRDM,3),1),0,2),size(modelRDM,2),size(modelRDM,3));
    end
    
    for isource = 1:length(atlas)
        neuralRDM(:,:,isource) = reshape(rescale(reshape(squeeze(neuralRDM(:,:,isource)),size(neuralRDM,1)*size(neuralRDM,2),1),0,2),size(neuralRDM,1),size(neuralRDM,2));    
        
        % Additionally, RDMs need to be centered per individual time point for PCA and regression
        % This happens below for the model RDM after the first PCA round, i.e., in Xregressout
        neuralRDM(:,:,isource) = squeeze(neuralRDM(:,:,isource)) - repmat(nanmean(squeeze(neuralRDM(:,:,isource))),size(neuralRDM,1),1);
    end
    
    %% dynamic RSA
    % simple correlation
    if cfg.similarity == 0
        
        % WARNING: SEARCHLIGHT NOT TESTED FOR SIMPLE CORRELATION, NEED TO ADJUST THIS BIT:
%         dRSA = zeros(newtimelength,size(cfg.peaks2test,1),length(atlas));
%         for iRDM = 1:length(cfg.models2test)% only test for our standard 7
%             
%             % neural - model correlation
%             dRSA(:,:,iRDM) = corr(squeeze(modelRDM(cfg.models2test(iRDM),:,:))',squeeze(neuralRDM(:,:,isource)));
%             
%         end
        
    % use regression to regress out principal components of other models
    elseif cfg.similarity == 1
        models2regressout = [2:10; 1 3:10; 1 2 4:10; 1:3 5:10; 1:4 6:10; 1:5 7:10; 1:6 8:10; 1:7 9:10; 1:8 10; 1:9];

        % to attenuate effects of model autocorrelation, we regress out the model itself, a certain distance away from our time point of interest
        % (i.e., ibin1 in the loops below, i.e., the time point of Xtest). This distancse is the regborder (regression border) variable, and it is
        % previously determined based on the simulations (see methods section of article for details), and loaded in here.
        load(fullfile(outdir,'..','regressionBorderPerModel_smRDM30msec'),'regborder');
        regborder.subinvarmods(10) = regborder.subvarmods(iSub);
        regborder = regborder.subinvarmods;
        
        dRSA = zeros(newtimelength,size(cfg.peaks2test,1),length(atlas));
        for ipeak = 1:size(cfg.peaks2test,1)% peaks from dRSA curves to run the searchlight on
            for ibin1 = 1:size(modelRDM,2)% model time
                
                iRDM = cfg.peaks2test(ipeak,1);
                lag2test = ibin1+round(cfg.peaks2test(ipeak,2)*cfg.downsample);% in samples
                                
                % indices of to-be-regressed out model RDMs for regressor selection: regress out in the -1.5 to 1.5 lag interval, and discarding indices outside of video time
                regidx = ibin1-cfg.downsample*1.5:ibin1+cfg.downsample*1.5;
                regidx(logical((regidx<1) + (regidx>size(modelRDM,2)))) = [];
                
                % indices at which the to-be-tested model itself will be regressed out, again discarding indices outside of video time
                test2regidx = [ibin1-cfg.downsample*1.5:ibin1-regborder(iRDM) ibin1+regborder(iRDM):ibin1+cfg.downsample*1.5];
                test2regidx(logical((test2regidx<1) + (test2regidx>size(modelRDM,2)))) = [];
                
                % the to-be-tested model RDM:
                Xtest = squeeze(modelRDM(iRDM,ibin1,:));
                
                % the to-be-tested model RDM at indices where it needs to be regressed out to attenuate model autocorrelation
                % and first PCA for each to-be-regressed out model to reduce dimensionality
                Xtest2regressout = squeeze(modelRDM(iRDM,test2regidx,:))';
                [~, score, ~, ~, exp, ~] = pca(Xtest2regressout);
                imax = sum(exp>.1);% only components with minimum variance of X%
                Xtest2regressout = score(:,1:imax);
                
                % now loop over all other models and also here run the first PCA on each to-be-regressed out model to reduce dimensionality
                Xregressout = zeros(size(modelRDM,3),500);% add some zeros for initialization
                for ireg = 1:nnz(models2regressout(iRDM,:))
                    
                    temp = squeeze(modelRDM(models2regressout(iRDM,ireg),regidx,:))';
                    [~, score, ~, ~, exp, ~] = pca(temp);
                    imax = sum(exp>.1);% only components with minimum variance of X%
                    score = score(:,1:imax);
                    Xregressout(:,nnz(Xregressout(1,:))+1:nnz(Xregressout(1,:))+size(score,2)) = score;
                    
                end
                Xregressout = Xregressout(:,1:nnz(Xregressout(1,:)));% remove additional zeros from the initialization
                    
                % combine to-be-regressed out models with to-be-tested model at indices where it needs to be regressed out
                Xregressout = [Xregressout Xtest2regressout];
                
                % center Xtest (actually not necessary because matlab's PCA implementation below does this automatically)                
                Xtest = Xtest - nanmean(Xtest);
                
                % rescale over all predictors in Xregressout at once to keep their differences in variance while keeping it close to Xtest and neuralRDM
                % additionally, center per individual predictor, again not necessary because PCA will do it below anyway
                Xregressout = reshape(rescale(reshape(Xregressout,size(Xregressout,1)*size(Xregressout,2),1),0,2),size(Xregressout,1),size(Xregressout,2));
                Xregressout = Xregressout-repmat(nanmean(Xregressout),size(Xregressout,1),1);
                                
                % make sure lag2test falls within stimulus window
                if lag2test<1 || lag2test>newtimelength
                    dRSA(ibin1,ipeak,:) = nan(length(atlas),1);% this lag is not defined at the edge of the stimulus window
                else
                    X = [Xtest Xregressout];
                    Y = squeeze(neuralRDM(:,lag2test,:));
                
                    % principal component regression (PCR): i.e., first PCA is run on X, resulting in the PCAscores (components), which are used as
                    % predictor variables in a least-squares regression, with Y as response variable. Last, the principal component regression weights
                    % betaPCR are projected back onto the original variable space using the PCA loadings, to extract a single regression weight
                    % corresponding to the original X
                    [PCALoadings,PCAScores] = pca(X,'NumComponents',75);
                    betaPCR = PCAScores\(Y-mean(Y));
                    temp = PCALoadings*betaPCR;
                
                    % select the first weight, which corresponds to the weight for Xtest. The other weights correspond to Xregressout, which we are
                    % not interested in for this particular analysis. Rather, they are 'regressed out'. 
                    dRSA(ibin1,ipeak,:) = temp(1,:);
                end               
                
            end% ibin1 loop
            
        end% model loop
        
    end
    
    dRSA = squeeze(nanmean(dRSA));% average over stimulus time. On edges there are some NaNs because some lags fall outside of resampled stimulus window. 
    
    % combine data from all iterations
    dRSAperIter(iter,:,:) = dRSA';
    clear dRSA
    
end

% average over batches
dRSAall = squeeze(mean(dRSAperIter));
clear dRSAperIter;

save(fn2save,'dRSAall','atlas');

end