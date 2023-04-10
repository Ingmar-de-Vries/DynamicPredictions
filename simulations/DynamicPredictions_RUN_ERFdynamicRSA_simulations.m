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
indirEye = fullfile(cfg.path,'data','MEG','ELpreprocessed');

if cfg.similarity == 0
    corrORglm = 'corr';
elseif cfg.similarity == 1
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
fsKin = 100;% sample rate kinematic recordings
tKin = 0:1/fsKin:5-1/fsKin;
zerolag = ceil(length(TimeVec2)./2);% for the simulations
triallengthsec = 5;
newtimelength = (triallengthsec-cfg.randshuff(2))*cfg.downsample;

%%  set weights for random neural RDM to which simulated model RDM will be added, if zero there is no random neural RDM
randomNeuralRDMweight = cfg.randweight;
simWeights = zeros(1,length(TimeVec2));
simWeights(dsearchn(TimeVec2',cfg.lag)) = 1;
simWeights = repmat(simWeights,99,1);
if simmod == 99
    simWeights = zeros(size(simWeights));
end

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
dRSAperIter = zeros(length(iterations),newtimelength,newtimelength,length(cfg.dynRDMnames));
for iter = iterations
    
    % determine indices for stimulus-specific re-alignment for current iteration
    shuffleIDmodel = dsearchn(tKin',shuffleTime(iter,:)')';
    
    %% prepare and add model RDM to random neural RDM to create simulated RDM (or without random neural RDM if cfg.randweight = 0)
    
    % randomize neural RDM as basis for simulated neural RDM, random uniformly distributed numbers between 0 and 2
    neuralRDM = 2*rand((nstim*nstim-nstim)./2,cfg.downsample*(triallengthsec-cfg.randshuff(2)));
    
    % smooth simulated neural RDM across time
    if cfg.smoothNeuralRDM
        neuralRDM = ft_preproc_smooth(neuralRDM,cfg.smoothNeuralRDM);
    end
    
    % create model RDM that will be implanted into the random neural RDM according to [lag X model] weight matrix in simWeights
    RDM2implant = zeros(nstim,nstim,cfg.downsample*(triallengthsec-cfg.randshuff(2)));
    for istim1 = 1:nstim
        for istim2 = 1:nstim
            
            tNewID1 = shuffleIDmodel(istim1):shuffleIDmodel(istim1)+(triallengthsec-cfg.randshuff(2))*cfg.downsample-1;
            tNewID2 = shuffleIDmodel(istim2):shuffleIDmodel(istim2)+(triallengthsec-cfg.randshuff(2))*cfg.downsample-1;
            
            lags = find(simWeights(implant,:) ~= 0);
            for ilag = 1:length(lags)% loop over different lags at which current model is implanted
                
                currentlag = lags(ilag);
                weight = simWeights(implant,currentlag);
                lag2implant = currentlag - zerolag;
                lag2implant = lag2implant / cfg.downsample;%change from samples to seconds
                
                temptrial = squeeze(dynRDM{implant}(istim1,istim2,tNewID1+lag2implant,tNewID2+lag2implant));
                
                temptrial = diag(temptrial);
                RDM2implant(istim1,istim2,:) = squeeze(RDM2implant(istim1,istim2,:)) + weight*temptrial;
            end% lag loop
            
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
    
    % rescale RDM2implant to same [0 2] interval, so has same scale as random neural data (and as each other, so dRSA values are better comparable)
    RDM2implant = reshape(rescale(reshape(RDM2implant,size(RDM2implant,1)*size(RDM2implant,2),1),0,2),size(RDM2implant,1),size(RDM2implant,2));
    
    % Implant model RDM into randomized neural RDM
    if simmod ~= 99% if simmod == 99, we only simulated random neural data without implanted model RDM, just to confirm we don't get any meaningful patterns in the dRSA results
        neuralRDM = randomNeuralRDMweight*neuralRDM + RDM2implant;
    end
    
    %% prepare model RDM for current iteration
    % re-align model RDMs for current iteration
    modelRDMsquare = zeros(length(dynRDM),nstim,nstim,cfg.downsample*(triallengthsec-cfg.randshuff(2)));%cfg.downsample*(triallengthsec+cfg.randshuff(2)));
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
    
    % smooth model RDM across time
    if cfg.smoothModelRDM
        for i=1:size(modelRDM,3)
            modelRDM(:,:,i) = ft_preproc_smooth(squeeze(modelRDM(:,:,i)),cfg.smoothModelRDM);
        end
    end
    
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
    if cfg.similarity == 0
        
        dRSA = zeros(size(neuralRDM,2),size(neuralRDM,2),length(cfg.dynRDMnames));
        for iRDM = 1:length(cfg.dynRDMnames)
            
            % neural - model correlation
            dRSA(:,:,iRDM) = corr(squeeze(modelRDM(iRDM,:,:))',neuralRDM);
            
        end
        
    elseif cfg.similarity == 1% principal component regression (PCR)
            
        % which models to regress out, i.e., all other models:
        models2regressout = [2:10; 1 3:10; 1 2 4:10; 1:3 5:10; 1:4 6:10; 1:5 7:10; 1:6 8:10; 1:7 9:10; 1:8 10; 1:9;];
        
        % to attenuate effects of model autocorrelation, we regress out the model itself, a certain distance away from our time point of interest
        % (i.e., ibin1 in the loops below, i.e., the time point of Xtest). This distancse is the regborder (regression border) variable, and it is
        % previously determined based on the simulations (see methods section of article for details), and loaded in here.
        load(fullfile(outdir,'..','..','regressionBorderPerModel_smRDM30msec'),'regborder');
        regborder.subinvarmods(10) = regborder.subvarmods(isub);
        regborder = regborder.subinvarmods;
                
        dRSA = zeros(size(neuralRDM,2),size(neuralRDM,2),length(cfg.dynRDMnames));
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
                dRSA(ibin1,:,iRDM) = temp(1,:);
                                
            end% ibin1 loop
            
        end% model loop
        
    end
    
    % combine data from all iterations
    dRSAperIter(iter-(iterbatch-1)*cfg.iterationsPERbatch,:,:,:) = dRSA;
    clear dRSA
    
end

% average over iterations within this batch
dRSAall = squeeze(mean(dRSAperIter));
clear dRSAperIter;

save(fn2save,'dRSAall');

end