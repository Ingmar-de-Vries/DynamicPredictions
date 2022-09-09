function DynamicPredictions_RUN_ERFdynamicRSA_simulations(cfg,simmod,eyeSub,iterbatch)

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
indirModel = fullfile(cfg.path,'data','modelRDMs');
indirKin = fullfile(cfg.path,'experiment','Stimuli','sequences_final','kinematics100Hz');
indirEye = fullfile(cfg.path,'data','MEG','ELpreprocessed');

if cfg.glmRSA == 0
    corrORglm = 'corr';
elseif cfg.glmRSA == 1
    corrORglm = ['pcr_' num2str(cfg.nPCAcomps) 'comps'];
end
outdir = fullfile(cfg.path,'data','MEG',['source_' cfg.atlas],'RSA', 'simulations', [corrORglm '_' num2str(cfg.lag*1000) 'lag_' num2str(cfg.randshuff(1)) 'iterations_' num2str(ceil(cfg.randshuff(2)*1000)) 'msec']);
if ~exist(outdir,'dir')
    mkdir(outdir);
end

fn2save = sprintf('%s%cdRSA_SUB%02d_implant%s_batch%d_%dHz_smMEG%d_smRDMneu%d_smRDMmod%d',outdir, filesep, eyeSub, num2str(simmod), iterbatch, cfg.downsample, cfg.smoothingMSec, cfg.smoothNeuralRDM, cfg.smoothModelRDM);

% check if this subject already ran, if so skip
if exist([fn2save '.mat'],'file')
    return
end

%% time definitions
TimeVec2 = -1:1/cfg.downsample:1;% time for lag-plot in seconds, used to define lag at which simulated model RDM will be 'implanted'
fsVid = 50;% sample rate videos (for video-based RDMs)
fsKin = 100;% sample rate kinematic recordings
tVid = 0:1/fsVid:5-1/fsVid;
tKin = 0:1/fsKin:5-1/fsKin;
zerolag = ceil(length(TimeVec2)./2);% for the simulations
triallengthsec = 5;
newtimelength = (triallengthsec-cfg.randshuff(2))*cfg.downsample;
% tNew and t4deriv are only used for interpolation of derivatives of kinematic models
tNew = 0:1/cfg.downsample:triallengthsec+cfg.randshuff(2)-1/cfg.downsample;% after resampling
t4deriv = tNew+(1/cfg.downsample)/2;% timesteps in between for any derivative (e.g., position --> motion, or motion --> acceleration)
t4deriv(end) = [];
modelID = 1:triallengthsec*cfg.downsample;

%%  set weights for random neural RDM to which simulated model RDM will be added, if zero there is no random neural RDM
randomNeuralRDMweight = cfg.randweight;
simWeights = zeros(1,length(TimeVec2));
simWeights(dsearchn(TimeVec2',cfg.lag)) = 1;
simWeights = repmat(simWeights,99,1);
if simmod == 99
    simWeights = zeros(size(simWeights));
end

%% load model RDMs
load([indirModel filesep 'DynamicPredictions_dynRDM_graysmooth'],'RDMgraysmooth');
load([indirModel filesep 'DynamicPredictions_dynRDM_opticalFlow'],'RDMopticalflow_dir','RDMopticalflow_mag');

dynRDM{1} = RDMgraysmooth;
dynRDM{2} = RDMopticalflow_mag;
dynRDM{3} = RDMopticalflow_dir;

% clear up workspace
clear RDMgraysmooth RDMopticalflow_mag RDMopticalflow_dir

% Kinematic models need to be computed inside the iteration loop, because
% procrustes alignment needs to be done after randomly re-aligning the data
% Here we first simply load the kinematic data once
kinNames = dir(fullfile(indirKin, '*.mat'));

seqs = zeros(length(kinNames),13,3,500);
for iseq=1:length(kinNames)
    
    load(fullfile(indirKin, kinNames(iseq).name),'Trajs');
    seqs(iseq,:,:,:) = Trajs;
    
end

% linear interpolation of models to meet the neural sample frequency
% kinematic marker data
if fsKin ~= cfg.downsample
    tempall = zeros(size(seqs,1),size(seqs,2),size(seqs,3),length(tKin));
    for istim = 1:size(seqs,1)
        for imarker = 1:size(seqs,2)
            for idim = 1:size(seqs,3)
                
                temp = squeeze(seqs(istim,imarker,idim,:));
                temp = interp1(tKin,temp,tKin,'pchip','extrap');
                tempall(istim,imarker,idim,:) = temp;
                
            end% dimensions loop
        end% marker loop
    end% stimuli loop
    
    seqs = tempall;
end
clear tempall

% video-based data
if fsVid ~= cfg.downsample
    for iRDM = [1 2 3]%10 if semantic is included
        
        tempall = zeros(size(dynRDM{1},1),size(dynRDM{1},1),length(tKin),length(tKin));
        for istim1 = 1:size(dynRDM{1},1)
            for istim2 = 1:size(dynRDM{1},1)
                
                temp = squeeze(dynRDM{iRDM}(istim1,istim2,:,:));
                % bilinear interpolation:
                temp = interp1(tVid,temp,tKin,'pchip','extrap');
                temp = interp1(tVid,temp',tKin,'pchip','extrap')';
                tempall(istim1,istim2,:,:) = temp;
                
            end% second stimuli loop
        end% first stimuli loop
        
        dynRDM{iRDM} = tempall;
    end% model loop
end% if loop
clear tempall

%% load and prepare eyetracker RDMs
fn= sprintf('%s%cSUB%02d_dynRDM_eyeTracker',indirEye, filesep, eyeSub);
load(fn,'RDMeyePOS');
dynRDM{10} = RDMeyePOS;
clear RDMeyePOS;

%% start iteration loop
% create matrix with times for randomly shuffling sequence onsets over X iterations
% shuffleTime = start times size [iterations X stimuli];

% make sure rand numbers are different for each simulated model RDM and each time this script is ran (i.e., for different iteration batches that are
% run in parallel we want to make sure they're not actually exact copies because they use the same random number stream!)
rng((round(sum(clock))+simmod*100+iterbatch)*100000);

% make sure it's a multiple of 0.01, which is the sampling rate, such that exact same time point is picked for neural and kinematic data
shuffleTime = round((rand(cfg.randshuff(1),nstim)*(cfg.randshuff(2)-1/cfg.downsample))./0.01).*0.01;

% loop over iterations in which data are randomly resampled
iterations = 1+(iterbatch-1)*cfg.iterationsPERbatch:iterbatch*cfg.iterationsPERbatch;% iterations in this batch

% initialize
omegaPerIter = zeros(length(iterations),newtimelength,newtimelength,length(cfg.dynRDMnames));
for iter = iterations
    
    % determine indices for stimulus-specific re-alignment for current iteration
    shuffleIDmodel = dsearchn(tKin',shuffleTime(iter,:)')';
    
    % tstart is where the stimulus model should start after shifting it in order to align it with others (relative to t = 0, which is cfg.randshuff(2)*cfg.downsample)
    tstart = cfg.randshuff(2)*cfg.downsample-shuffleIDmodel+1;
    
    %% prepare and add model RDM to random neural RDM to create simulated RDM (or without random neural RDM if cfg.randweight = 0)
    
    % randomize neural RDM as basis for simulated neural RDM, random uniformly distributed numbers between 0 and 2
    neuralRDM = 2*rand((nstim*nstim-nstim)./2,cfg.downsample*(triallengthsec+cfg.randshuff(2)));
    
    % create model RDM that will be implanted into the random neural RDM according to [lag X model] weight matrix in simWeights
    RDM2implant = zeros(nstim,nstim,cfg.downsample*(triallengthsec+cfg.randshuff(2)));
    for istim1 = 1:nstim
        for istim2 = 1:nstim
            
            % realignment happens by creating a new time vector (tKin4extrap) which is used for interpolating data to new time
            % advantage of this instead of simply selecting the 3-sec data, is that you can keep the non-selected data for a little longer as padding,
            % so there's less chance of potential edge artifacts, e.g., when temporally smoothing the RDMs. Padding will be cut later on. This
            % inevitably leads to some time window without data, in which now data will be interpolated. This is meaningless data, but this doesn't
            % matter as it is only part of the padding that will be cut off anyway
            tKin4extrap1 = [-tstart(istim1)/cfg.downsample:1/cfg.downsample:tKin(1) tKin(2:end-1) tKin(end):1/cfg.downsample:tKin(end)+(cfg.randshuff(2)*cfg.downsample-tstart(istim1))/cfg.downsample];
            tKin4extrap2 = [-tstart(istim2)/cfg.downsample:1/cfg.downsample:tKin(1) tKin(2:end-1) tKin(end):1/cfg.downsample:tKin(end)+(cfg.randshuff(2)*cfg.downsample-tstart(istim2))/cfg.downsample];
            
            if simmod < 4 || simmod > 9% dynRDM already computed, just extract correct windows
                
                % loop over different lags at which current model is simulated 
                lags = find(simWeights(simmod,:) ~= 0);
                for ilag = 1:length(lags)
                    
                    currentlag = lags(ilag);
                    weight = simWeights(simmod,currentlag);
                    lag2implant = currentlag - zerolag;
                    lag2implant = lag2implant / cfg.downsample;%change from samples to seconds
                    
                    % select time window for resampling + simulation lag
                    temptrial = interp1(tKin,squeeze(dynRDM{simmod}(istim1,istim2,modelID,modelID)),tKin4extrap1-lag2implant,'pchip','extrap');
                    temptrial = interp1(tKin,temptrial',tKin4extrap2-lag2implant,'pchip','extrap')';
                    
                    % keep extrapolated values within boundaries, because they become extreme sometimes, which could be a problem for rescaling the RDMs below
                    % find max value in the 'real data window', i.e., without the padding which contains extreme extrapolated values
                    maxval = max(max(temptrial(tstart(istim1):tstart(istim1)+length(modelID),tstart(istim2):tstart(istim2)+length(modelID))));
                    % find min value in the 'real data window', i.e., without the padding which contains extreme extrapolated values
                    minval = min(min(temptrial(tstart(istim1):tstart(istim1)+length(modelID),tstart(istim2):tstart(istim2)+length(modelID))));
                    temptrial(temptrial > maxval) = maxval;
                    temptrial(temptrial < minval) = minval;
                    
                    % diagonal is now what we're looking for, i.e., the dissimilarity between stim1 and stim2 at each time point after
                    % re-alignment on this particular iteration
                    temptrial = diag(temptrial);
                    
                    RDM2implant(istim1,istim2,:) = squeeze(RDM2implant(istim1,istim2,:)) + weight*temptrial;
                end% lag loop
                
            elseif simmod > 3 && simmod < 10 && nnz(simWeights(simmod,:)) ~= 0% for if dynRDM still needs to be created, i.e. for kinematic models
                
                lags = find(simWeights(simmod,:) ~= 0);% loop over different lags at which at least one of the kinematic models is implanted
                for ilag = 1:length(lags)
                    
                    currentlag = lags(ilag);
                    weight = simWeights(simmod,currentlag);
                    lag2implant = currentlag - zerolag;
                    lag2implant = lag2implant / cfg.downsample;%change from samples to seconds
                    
                    % kinematic RDMs still need to be computed
                    % again, use interpolation for alignment
                    seq1 = permute(interp1(tKin,permute(squeeze(seqs(istim1,:,:,modelID)),[3 1 2]),tKin4extrap1-lag2implant,'pchip','extrap'),[2 3 1]);
                    seq2 = permute(interp1(tKin,permute(squeeze(seqs(istim2,:,:,modelID)),[3 1 2]),tKin4extrap2-lag2implant,'pchip','extrap'),[2 3 1]);
                    
                    % again, keep extrapolated values within boundaries:
                    for imark = 1:size(seq1,1)
                        for idim = 1:size(seq1,2)
                            maxval = max(seq1(imark,idim,tstart(istim1):tstart(istim1)+length(modelID)));
                            minval = min(seq1(imark,idim,tstart(istim1):tstart(istim1)+length(modelID)));
                            seq1(imark,idim,seq1(imark,idim,:)>maxval) = maxval;
                            seq1(imark,idim,seq1(imark,idim,:)<minval) = minval;
                            maxval = max(seq2(imark,idim,tstart(istim2):tstart(istim2)+length(modelID)));
                            minval = min(seq2(imark,idim,tstart(istim2):tstart(istim2)+length(modelID)));
                            seq2(imark,idim,seq2(imark,idim,:)>maxval) = maxval;
                            seq2(imark,idim,seq2(imark,idim,:)<minval) = minval;
                        end
                    end
                    
                    % Align sequence 2 to sequence 1 using modified version of procrustes
                    % note that my modified version of procrustes is very time inefficient now, as it loops over angles. There should be a much
                    % faster implementation possible, still need to work on this. But currently works.
                    seq2trans = zeros(size(seq2));
                    for iframe = 1:size(seq2,3)
                        [~, seq2trans(:,:,iframe)] = procrustes_constrain_rotationZaxis_IdV(squeeze(seq1(:,:,iframe)),squeeze(seq2(:,:,iframe)),'Reflection',false,'Scaling',false);
                    end
                    
                    % Position vectors
                    pos1 = reshape(seq1,numel(seq1(:,:,1)),size(seq1,3));
                    pos2 = reshape(seq2,numel(seq2(:,:,1)),size(seq2,3));
                    pos2trans = reshape(seq2trans,numel(seq2trans(:,:,1)),size(seq2trans,3));
                    
                    if simmod == 4
                        % absolute position
                        temp = weight*(1 - diag(corr(pos1,pos2)));
                        temp(isnan(temp)) = 0;% if correlation over zero-padding, results in NaN, which will mess everything up. Better zeros here.
                        RDM2implant(istim1,istim2,:) = squeeze(RDM2implant(istim1,istim2,:)) + temp;
                    elseif simmod == 5
                        % relative position
                        temp = weight*(1 - diag(corr(pos1,pos2trans)));
                        temp(isnan(temp)) = 0;% if correlation over zero-padding, results in NaN, which will mess everything up. Better zeros here.
                        RDM2implant(istim1,istim2,:) = squeeze(RDM2implant(istim1,istim2,:)) + temp;
                    end
                    
                    %% Velocity (motion) vectors: original order
                    if simmod > 5
                        vel1 = interp1(t4deriv,diff(pos1,[],2)',tNew,'pchip','extrap')';
                        vel2 = interp1(t4deriv,diff(pos2,[],2)',tNew,'pchip','extrap')';
                        vel2trans = interp1(t4deriv,diff(pos2trans,[],2)',tNew,'pchip','extrap')';
                    end
                    
                    if simmod == 6
                        % absolute motion
                        temp = weight*(1 - diag(corr(vel1,vel2)));
                        temp(isnan(temp)) = 0;% if correlation over zero-padding, results in NaN, which will mess everything up. Better zeros here.
                        RDM2implant(istim1,istim2,:) = squeeze(RDM2implant(istim1,istim2,:)) + temp;
                    elseif simmod == 7
                        % relative motion
                        temp = weight*(1 - diag(corr(vel1,vel2trans)));
                        temp(isnan(temp)) = 0;% if correlation over zero-padding, results in NaN, which will mess everything up. Better zeros here.
                        RDM2implant(istim1,istim2,:) = squeeze(RDM2implant(istim1,istim2,:)) + temp;
                    end
                    
                    %% Acceleration vectors: original order
                    if simmod > 7
                        acc1 = interp1(t4deriv,diff(vel1,[],2)',tNew,'pchip','extrap')';
                        acc2 = interp1(t4deriv,diff(vel2,[],2)',tNew,'pchip','extrap')';
                        acc2trans = interp1(t4deriv,diff(vel2trans,[],2)',tNew,'pchip','extrap')';
                    end
                    
                    if simmod == 8
                        % absolute acceleration
                        temp = weight*(1 - diag(corr(acc1,acc2)));
                        temp(isnan(temp)) = 0;% if correlation over zero-padding, results in NaN, which will mess everything up. Better zeros here.
                        RDM2implant(istim1,istim2,:) = squeeze(RDM2implant(istim1,istim2,:)) + temp;
                    elseif simmod == 9
                        % relative acceleration
                        temp = weight*(1 - diag(corr(acc1,acc2trans)));
                        temp(isnan(temp)) = 0;% if correlation over zero-padding, results in NaN, which will mess everything up. Better zeros here.
                        RDM2implant(istim1,istim2,:) = squeeze(RDM2implant(istim1,istim2,:)) + temp;
                    end
                    
                end% lag loop
            end% if iRDM
            
        end% second stimulus loop
    end% first stimulus loop
    
    % for relative models take average of (stim1*stim2trans) and (stim1trans*stim2), which is not exactly the same because procrustes is not symmetrical
    if simmod == 5 || simmod == 7 || simmod == 9
        for itime = 1:newtimelength
            temp = squeeze(RDM2implant(:,:,itime));
            RDM2implant(:,:,itime) = (temp + temp')/2;
        end% time loop
    end
    
    % extract vector RDM from triangle of square RDM
    temp = zeros((nstim*nstim-nstim)/2,size(RDM2implant,3));
    for iBin = 1:size(RDM2implant,3)
        temp(:,iBin) = squareform(tril(squeeze(RDM2implant(:,:,iBin)),-1));
    end
    RDM2implant = temp;clear temp
    
    % smooth implant RDM across time
    if cfg.smoothModelRDM
        RDM2implant = ft_preproc_smooth(RDM2implant,cfg.smoothModelRDM);
    end
    
    % cut extra padding used for dRSA
    % because outside those bounds some RDM vectors are not defined (i.e., they contain zeros, or interpolated values)
    RDM2implant = RDM2implant(:,cfg.randshuff(2)*cfg.downsample:cfg.randshuff(2)*cfg.downsample+newtimelength-1);
    neuralRDM = neuralRDM(:,cfg.randshuff(2)*cfg.downsample:cfg.randshuff(2)*cfg.downsample+newtimelength-1);
    
    % rescale RDM2implant to same [0 2] interval, so has same scale as random neural data (and as each other, so dRSA values are better comparable)
    for ibin = 1:size(RDM2implant,2)
        RDM2implant(:,ibin) = rescale(squeeze(RDM2implant(:,ibin)),0,2);
    end
    
    % Implant model RDM into randomized neural RDM
    if simmod ~= 99% if simmod == 99, we only simulated random neural data without implanted model RDM, just to confirm we don't get any meaningful patterns in the dRSA results
        neuralRDM = randomNeuralRDMweight*neuralRDM + RDM2implant;
    end
    
    %% prepare and re-align model RDM for current iteration
    modelRDMsquare = zeros(length(cfg.dynRDMnames),nstim,nstim,cfg.downsample*(triallengthsec+cfg.randshuff(2)));
    % 1-3 are based on video data, > 4 are based on kinematic data, 10 is based on eye position data
    for iRDM = [1:4 10]
        for istim1 = 1:nstim
            for istim2 = 1:nstim                
                
                % realignment happens by creating a new time vector (tKin4extrap) which is used for interpolating data to new time
                % advantage of this instead of simply selecting the 3-sec data, is that you can keep the non-selected data for a little longer as padding,
                % so there's less chance of potential edge artifacts, e.g., when temporally smoothing the RDMs. Padding will be cut later on. This
                % inevitably leads to some time window without data, in which now data will be interpolated. This is meaningless data, but this doesn't
                % matter as it is only part of the padding that will be cut off anyway
                tKin4extrap1 = [-tstart(istim1)/cfg.downsample:1/cfg.downsample:tKin(1) tKin(2:end-1) tKin(end):1/cfg.downsample:tKin(end)+(cfg.randshuff(2)*cfg.downsample-tstart(istim1))/cfg.downsample];
                tKin4extrap2 = [-tstart(istim2)/cfg.downsample:1/cfg.downsample:tKin(1) tKin(2:end-1) tKin(end):1/cfg.downsample:tKin(end)+(cfg.randshuff(2)*cfg.downsample-tstart(istim2))/cfg.downsample];
                
                if iRDM < 4 || iRDM > 9% dynRDM already computed, just extract correct windows
                    
                    % select time window
                    temptrial = interp1(tKin,squeeze(dynRDM{iRDM}(istim1,istim2,modelID,modelID)),tKin4extrap1,'pchip','extrap');
                    temptrial = interp1(tKin,temptrial',tKin4extrap2,'pchip','extrap')';
                    
                    % keep extrapolated values within boundaries, because they become extreme sometimes, which could be a problem for rescaling the RDMs below
                    % find max value in the 'real data window', i.e., without the padding which contains extreme extrapolated values
                    maxval = max(max(temptrial(tstart(istim1):tstart(istim1)+length(modelID),tstart(istim2):tstart(istim2)+length(modelID))));
                    % find min value in the 'real data window', i.e., without the padding which contains extreme extrapolated values
                    minval = min(min(temptrial(tstart(istim1):tstart(istim1)+length(modelID),tstart(istim2):tstart(istim2)+length(modelID))));
                    temptrial(temptrial > maxval) = maxval;
                    temptrial(temptrial < minval) = minval;
                    
                    % diagonal is now what we're looking for, i.e., the dissimilarity between stim1 and stim2 at each time point after
                    % re-alignment on this particular iteration
                    temptrial = diag(temptrial);
                    
                    modelRDMsquare(iRDM,istim1,istim2,:) = temptrial;
                    
                else% for kinematic dynRDM is not yet computed (can only be done here due to relative pos/vel/acc/jerk)
                    
                    % kinematic RDMs still need to be computed
                    % again, use interpolation for alignment
                    seq1 = permute(interp1(tKin,permute(squeeze(seqs(istim1,:,:,modelID)),[3 1 2]),tKin4extrap1,'pchip','extrap'),[2 3 1]);
                    seq2 = permute(interp1(tKin,permute(squeeze(seqs(istim2,:,:,modelID)),[3 1 2]),tKin4extrap2,'pchip','extrap'),[2 3 1]);
                    
                    % again, keep extrapolated values within boundaries:
                    for imark = 1:size(seq1,1)
                        for idim = 1:size(seq1,2)
                            maxval = max(seq1(imark,idim,tstart(istim1):tstart(istim1)+length(modelID)));
                            minval = min(seq1(imark,idim,tstart(istim1):tstart(istim1)+length(modelID)));
                            seq1(imark,idim,seq1(imark,idim,:)>maxval) = maxval;
                            seq1(imark,idim,seq1(imark,idim,:)<minval) = minval;
                            maxval = max(seq2(imark,idim,tstart(istim2):tstart(istim2)+length(modelID)));
                            minval = min(seq2(imark,idim,tstart(istim2):tstart(istim2)+length(modelID)));
                            seq2(imark,idim,seq2(imark,idim,:)>maxval) = maxval;
                            seq2(imark,idim,seq2(imark,idim,:)<minval) = minval;
                        end
                    end
                    
                    % Align sequence 2 to sequence 1 using modified version of procrustes
                    % note that my modified version of procrustes is very time inefficient now, as it loops over angles. There should be a much
                    % faster implementation possible, still need to work on this. But currently works.
                    seq2trans = zeros(size(seq2));
                    for iframe = 1:size(seq2,3)
                        [~, seq2trans(:,:,iframe)] = procrustes_constrain_rotationZaxis_IdV(squeeze(seq1(:,:,iframe)),squeeze(seq2(:,:,iframe)),'Reflection',false,'Scaling',false);
                    end
                    
                    % Position vectors
                    pos1 = reshape(seq1,numel(seq1(:,:,1)),size(seq1,3));
                    pos2 = reshape(seq2,numel(seq2(:,:,1)),size(seq2,3));
                    pos2trans = reshape(seq2trans,numel(seq2trans(:,:,1)),size(seq2trans,3));
                    
                    % compute RDMs
                    modelRDMsquare(4,istim1,istim2,:) = 1 - diag(corr(pos1,pos2));% absolute position
                    modelRDMsquare(5,istim1,istim2,:) = 1 - diag(corr(pos1,pos2trans));% relative position
                    
                    %% Velocity (motion) vectors
                    vel1 = interp1(t4deriv,diff(pos1,[],2)',tNew,'pchip','extrap')';
                    vel2 = interp1(t4deriv,diff(pos2,[],2)',tNew,'pchip','extrap')';
                    vel2trans = interp1(t4deriv,diff(pos2trans,[],2)',tNew,'pchip','extrap')';
                    
                    % compute RDMs
                    modelRDMsquare(6,istim1,istim2,:) = 1 - diag(corr(vel1,vel2));% absolute motion
                    modelRDMsquare(7,istim1,istim2,:) = 1 - diag(corr(vel1,vel2trans));% relative motion
                    
                    %% Acceleration vectors
                    acc1 = interp1(t4deriv,diff(vel1,[],2)',tNew,'pchip','extrap')';
                    acc2 = interp1(t4deriv,diff(vel2,[],2)',tNew,'pchip','extrap')';
                    acc2trans = interp1(t4deriv,diff(vel2trans,[],2)',tNew,'pchip','extrap')';
                    
                    % compute RDMs
                    modelRDMsquare(8,istim1,istim2,:) = 1 - diag(corr(acc1,acc2));% absolute
                    modelRDMsquare(9,istim1,istim2,:) = 1 - diag(corr(acc1,acc2trans));% relative
                    
                end
                
            end% second stimulus loop
        end% first stimulus loop
    end% model loop
    
    % for relative models take average of (stim1*stim2trans) and (stim1trans*stim2), which is not exactly the same because procrustes is not symmetrical
    for iRDM = [5 7 9]
        for itime = 1:newtimelength
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
    
    % remove NaNs due to combination of extrapolation with constant value and taking the derivative for motion/acceleration, which results in
    % zeros, which in turn results in NaNs when performing correlation for RDMs (but only at the edges that are not chopped off below anyway)
    modelRDM(isnan(modelRDM)) = 0;
    
    % smooth model RDM across time
    if cfg.smoothModelRDM
        for i=1:size(modelRDM,3)
            modelRDM(:,:,i) = ft_preproc_smooth(squeeze(modelRDM(:,:,i)),cfg.smoothModelRDM);
        end
    end
    
    %% cut extra padding used for dRSA
    % because outside those bounds some RDM vectors are not defined (i.e., they contain zeros, or interpolated values)
    modelRDM = modelRDM(:,cfg.randshuff(2)*cfg.downsample:cfg.randshuff(2)*cfg.downsample+newtimelength-1,:);
    
    %% rescale all RDMs to same [0 2] interval. Unscaled might be problematic for PCA or regression (i.e., larger scale = more variance = higher component)
    % scale once across all time points but per model to keep internal model structure of variance over time intact
    for iRDM = 1:size(modelRDM,1)
        modelRDM(iRDM,:,:) = reshape(rescale(reshape(squeeze(modelRDM(iRDM,:,:)),size(modelRDM,2)*size(modelRDM,3),1),0,2),size(modelRDM,2),size(modelRDM,3));
    end
    neuralRDM = reshape(rescale(reshape(neuralRDM,size(neuralRDM,1)*size(neuralRDM,2),1),0,2),size(neuralRDM,1),size(neuralRDM,2));
    
    % Additionally, RDMs need to be centered per individual time point for PCA and regression
    % This happens below for the model RDM after the first PCA round, i.e., in Xregressout
    neuralRDM = neuralRDM - repmat(mean(neuralRDM),size(neuralRDM,1),1);
    
    %% dynamic RSA
    % simple correlation
    if cfg.glmRSA == 0
        
        omega = zeros(size(neuralRDM,2),size(neuralRDM,2),length(cfg.dynRDMnames));
        for iRDM = 1:length(cfg.dynRDMnames)
            
            % neural - model correlation
            omega(:,:,iRDM) = corr(squeeze(modelRDM(iRDM,:,:))',neuralRDM);
            
        end
        
    elseif cfg.glmRSA == 1% principal component regression (PCR)
            
        % which models to regress out, i.e., all other models:
        models2regressout = [2:10; 1 3:10; 1 2 4:10; 1:3 5:10; 1:4 6:10; 1:5 7:10; 1:6 8:10; 1:7 9:10; 1:8 10; 1:9;];
        
        % to attenuate effects of model autocorrelation, we regress out the model itself, a certain distance away from our time point of interest
        % (i.e., ibin1 in the loops below, i.e., the time point of Xtest). This distancse is the regborder (regression border) variable, and it is
        % previously determined based on the simulations (see methods section of article for details), and loaded in here.
        load(fullfile(outdir,'..','..','regressionBorderPerModel_smRDM20msec'),'regborder');
        
        % The rows correspond to amount of autocorrelation still present at a certain lag, as can be seen in Figure S4a
        % 1 = 0.5 autocorrelation, 2 = 0.25, 3 = 0.32 (which corresponds to 10% shared variance), 4 = 0.71 (which corresponds to 50% shared variance)
        regborder = regborder(:,3);
                
        omega = zeros(size(neuralRDM,2),size(neuralRDM,2),length(cfg.dynRDMnames));
        for iRDM = 1:length(cfg.dynRDMnames)% models to test
            for ibin1 = 1:size(modelRDM,2)% model time
                
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
                
                X = [Xtest Xregressout];
                Y = neuralRDM;
                
                % principal component regression (PCR): i.e., first PCA is run on X, resulting in the PCAscores (components), which are used as
                % predictor variables in a least-squares regression, with Y as response variable. Last, the principal component regression weights
                % betaPCR are projected back onto the original variable space using the PCA loadings, to extract a single regression weight
                % corresponding to the original X
                [PCALoadings,PCAScores] = pca(X,'NumComponents',cfg.nPCAcomps);
                betaPCR = PCAScores\(Y-mean(Y));
                temp = PCALoadings*betaPCR;
                
                % select the first weight, which corresponds to the weight for Xtest. The other weights correspond to Xregressout, which we are
                % not interested in for this particular analysis. Rather, they are 'regressed out'.
                omega(ibin1,:,iRDM) = temp(1,:);
                                
            end% ibin1 loop
            
        end% model loop
        
    end
    
    % combine data from all iterations
    omegaPerIter(iter-(iterbatch-1)*cfg.iterationsPERbatch,:,:,:) = omega;
    clear omega
    
end

% average over iterations within this batch
omegaAll = squeeze(mean(omegaPerIter));
clear omegaPerIter;

save(fn2save,'omegaAll');

end
