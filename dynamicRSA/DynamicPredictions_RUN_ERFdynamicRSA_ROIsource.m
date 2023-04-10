function DynamicPredictions_RUN_ERFdynamicRSA_ROIsource(cfg,iSub,iROI)

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
indir = fullfile(cfg.path,'data','MEG','PreprocessedSensor');% path to preprocessed sensor level MEG data 
indirModel = fullfile(cfg.path,'data','modelRDMs');% path to model RDMs that have been previously computed
indirEye = fullfile(cfg.path,'data','MEG','ELpreprocessed');% path to RDMs based on eyetracker data 

subfilz = dir(fullfile(indir,'preprocessedSensor*'));
load(fullfile(cfg.path,'code','neuralDecoding','ROIdefinitions'),'ROIdefinition');

if cfg.similarity == 0
    corrORglm = 'corr';
elseif cfg.similarity == 1
    corrORglm = ['pcr_' num2str(cfg.nPCAcomps) 'comps'];
end
outdir = fullfile(cfg.path,'data','MEG',['source_' cfg.atlas],'RSA', [corrORglm '_' num2str(cfg.randshuff(1)) 'iterations_' num2str(ceil(cfg.randshuff(2)*1000)) 'msec']);
if ~exist(outdir,'dir')
    mkdir(outdir);
end

% Filename to save
fn2save = sprintf('%s%cdRSA_SUB%02d_%dHz_%s_smMEG%d_smRDMneu%d_smRDMmod%d',outdir, filesep, iSub, cfg.downsample, ROIdefinition.names{iROI}{:}, cfg.smoothingMSec, cfg.smoothNeuralRDM, cfg.smoothModelRDM);

% check if this subject already ran, if so skip
if exist([fn2save '.mat'],'file')
    return
end

% paths to kernel and atlas definitions in Brainstorm database
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

%% time definitions
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
atlasIdx = contains(extractfield(atlas.Atlas,'Name'),cfg.atlas);
atlas = atlas.Atlas(atlasIdx).Scouts;

% do some atlas fixing:
% 1) in Schaefer there is also a background
% 2) in HCP some of the L and R to indicate hemisphere behind the parcel names are gone
idx2remove = contains(extractfield(atlas,'Label'),'Background');
atlas(idx2remove) = [];

for iparcel = 1:length(atlas)
    
    if ~contains(atlas(iparcel).Label,' L') && ~contains(atlas(iparcel).Label,' R')
        
        if strcmp(atlas(iparcel).Region(1),'L')
            atlas(iparcel).Label = [atlas(iparcel).Label ' L'];
        elseif strcmp(atlas(iparcel).Region(1),'R')
            atlas(iparcel).Label = [atlas(iparcel).Label ' R'];
        end
        
    end
    
end

targets = ds.sa.targets;% always gets lost when converting to FT format
ds=cosmo_map2meeg(ds);% convert data to FT format

% select only gradiometers because MNE is only applied to those
ds.label = ds.label(chanIdx);
ds.trial = ds.trial(:,chanIdx,:);

% select ROI
ROIsel = [];
ROIsel.names = ROIdefinition.names(iROI);
ROIsel.parcels = ROIdefinition.parcels(iROI);
new = ROIsel;

% divide in left and right hemisphere
for isource = 1:length(ROIsel.names)
    
    for iparcel = 1:length(ROIsel.parcels{isource})
        
        new.parcels{isource}{iparcel+length(ROIsel.parcels{isource})} = [ROIsel.parcels{isource}{iparcel} ' R'];
        new.parcels{isource}{iparcel} = [ROIsel.parcels{isource}{iparcel} ' L'];
        
    end
    
end
ROIsel = new; clear new;

% select parcels for current ROI, select kernel for those, and apply to data
vertices = [];
for iparcel = 1:length(ROIsel.parcels{isource})
    parcelidx = find(strcmp(extractfield(atlas,'Label'),ROIsel.parcels{isource}{iparcel}),1);
    
    if isempty(parcelidx)
        error(['cannot find matching label in atlas for ' ROIsel.parcels{isource}{iparcel}]);
    end
    
    vertices = [vertices atlas(parcelidx).Vertices(:).'];
    
end

ROIsel.vertices{isource} = vertices;

% Then apply conversion kernel to sensor level data: 
kernelSel = kernel(vertices,:);
sourcedata = zeros(size(ds.trial,1),size(kernelSel,1),size(ds.trial,3));
for itrial = 1:size(ds.trial,1)
    temp = squeeze(ds.trial(itrial,:,:));
    sourcedata(itrial,:,:) = kernelSel*temp;
end
ds.trial = sourcedata;
ds.label = vertices';

ds = cosmo_meeg_dataset(ds);% convert data back to cosmo format
ds.sa.targets = targets;% always gets lost when converting to FT format

%% start iteration loop
% create matrix with times for randomly shuffling sequence onsets over X iterations
% shuffleTime = start times size [iterations X stimuli];

% make sure random numbers are different for each subject, each ROI, and each time this script is ran at a different time and date
rng((round(sum(clock))+iSub*1000+iROI)*100000);

% make sure it's a multiple of 0.01, which is the sampling rate, such that exact same time point is picked for neural and kinematic data
shuffleTime = round((rand(cfg.randshuff(1),nstim)*(cfg.randshuff(2)-1/cfg.downsample))./0.01).*0.01;

% divide jitter iteration loop batches of 50 iterations and average over each batch to save on RAM storage
batchnum = round(cfg.randshuff(1)/50);
dRSAperBatch = zeros(batchnum,newtimelength,newtimelength,length(cfg.models2test));

% loop over iterations in which data are temporally jittered
for iterbatch = 1:batchnum
    
    % iterations in this batch of 50 iterations
    iterations = 1+(iterbatch-1)*50:iterbatch*50;
    
    dRSAperIter = zeros(length(iterations),newtimelength,newtimelength,length(cfg.models2test));
    for iter = iterations
        
        % determine indices for stimulus-specific re-alignment for current iteration
        shuffleIDmodel = dsearchn(tMod',shuffleTime(iter,:)')';
        shuffleIDneural = dsearchn(tNeural',shuffleTime(iter,:)')';     
        
        %% prepare neural RDM for current iteration
        ds_iter = ds;
        
        % re-align neural data to stimulus-specific temporal jitter
        ds_iter=cosmo_map2meeg(ds_iter);% convert data to FT format
        
        temptrial = zeros(nstim,length(ds_iter.label),cfg.downsample*(triallengthsec-cfg.randshuff(2)));
        for istim = 1:nstim
            
            tNewID = shuffleIDneural(istim):shuffleIDneural(istim)+(triallengthsec-cfg.randshuff(2))*cfg.downsample-1;
            
            temptrial(istim,:,:) = squeeze(ds_iter.trial(istim,:,tNewID));
            
        end
        ds_iter.trial = temptrial; clear temptrial
        ds_iter.time = 0:1/cfg.downsample:triallengthsec-cfg.randshuff(2)-1/cfg.downsample;
        ds_iter.trialinfo = 1:nstim;
        
        ds_iter = cosmo_meeg_dataset(ds_iter);% convert data back to cosmo format
        ds_iter.sa.targets = [1:nstim]';% for some reason this field gets lost when coverting back and forth from cosmo to FT format
        
        % define neighborhood over time, here zero
        nbrhood=cosmo_interval_neighborhood(ds_iter,'time','radius',0);
        
        % compute neural RDM
        ds_res=cosmo_searchlight(ds_iter,nbrhood,measure,measure_args);
        
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
        
        neuralRDM = zeros((nstim*nstim-nstim)/2,size(MATxTIME,3));
        for iBin = 1:size(MATxTIME,3)
            matrices = squeeze(MATxTIME(:,:,iBin));            
            matrices = matrices / sum(mean(matrices,2)); % normalize, shouldn't do anything to final result
            
            % simplify matrix by making it symmetrical, which it already is anyway...
            DSM = (matrices + matrices')/2;
            
            DSM(logical(eye(size(DSM)))) = 0;
            neuralRDM(:,iBin) = squareform(DSM);
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
        % scale once across all time points but per model to keep internal model structure of variance over time intact
        for iRDM = 1:size(modelRDM,1)
            modelRDM(iRDM,:,:) = reshape(rescale(reshape(squeeze(modelRDM(iRDM,:,:)),size(modelRDM,2)*size(modelRDM,3),1),0,2),size(modelRDM,2),size(modelRDM,3));
        end        
        neuralRDM = reshape(rescale(reshape(neuralRDM,size(neuralRDM,1)*size(neuralRDM,2),1),0,2),size(neuralRDM,1),size(neuralRDM,2));
        
        % Additionally, RDMs need to be centered per individual time point for PCA and regression
        % This happens below for the model RDM after the first PCA round, i.e., in Xregressout
        neuralRDM = neuralRDM - repmat(nanmean(neuralRDM),size(neuralRDM,1),1);
        
        %% dynamic RSA
        % simple correlation
        if cfg.similarity == 0
            
            dRSA = zeros(size(neuralRDM,2),size(neuralRDM,2),length(cfg.models2test));
            for iRDM = 1:length(cfg.models2test)
                
                % neural - model correlation
                dRSA(:,:,iRDM) = corr(squeeze(modelRDM(cfg.models2test(iRDM),:,:))',neuralRDM);
                
            end
            
        % use regression to regress out principal components of other models
        elseif cfg.similarity == 1
            
            % which models to regress out, i.e., all other models:
            models2regressout = [2:10; 1 3:10; 1 2 4:10; 1:3 5:10; 1:4 6:10; 1:5 7:10; 1:6 8:10; 1:7 9:10; 1:8 10; 1:9];
            models2regressout = models2regressout(cfg.models2test,:);
            
            % to attenuate effects of model autocorrelation, we regress out the model itself, a certain distance away from our time point of interest
            % (i.e., ibin1 in the loops below, i.e., the time point of Xtest). This distancse is the regborder (regression border) variable, and it is
            % previously determined based on the simulations (see methods section of article for details), and loaded in here. 
            load(fullfile(outdir,'..','regressionBorderPerModel_smRDM30msec'),'regborder');
            
            % select subject-variant model for current subject (i.e., eyetracker)
            regborder.subinvarmods(10) = regborder.subvarmods(isub);
            regborder = regborder.subinvarmods;
            regborder = regborder(cfg.models2test);
                        
            dRSA = zeros(size(neuralRDM,2),size(neuralRDM,2),length(cfg.models2test));
            for iRDM = 1:length(cfg.models2test)% models to test
                for ibin1 = 1:size(modelRDM,2)% model time
                    
                    % indices of to-be-regressed out model RDMs for regressor selection: regress out in the -1.5 to 1.5 lag interval, and discarding indices outside of video time
                    regidx = ibin1-cfg.downsample*1.5:ibin1+cfg.downsample*1.5;
                    regidx(logical((regidx<1) + (regidx>size(modelRDM,2)))) = [];
                    
                    % indices at which the to-be-tested model itself will be regressed out, again discarding indices outside of video time
                    test2regidx = [ibin1-cfg.downsample*1.5:ibin1-regborder(iRDM) ibin1+regborder(iRDM):ibin1+cfg.downsample*1.5];
                    test2regidx(logical((test2regidx<1) + (test2regidx>size(modelRDM,2)))) = [];
                    
                    % the to-be-tested model RDM:
                    Xtest = squeeze(modelRDM(cfg.models2test(iRDM),ibin1,:));
                    
                    % the to-be-tested model RDM at indices where it needs to be regressed out to attenuate model autocorrelation 
                    % and first PCA for each to-be-regressed out model to reduce dimensionality
                    Xtest2regressout = squeeze(modelRDM(cfg.models2test(iRDM),test2regidx,:))';
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
        dRSAperIter(iter-(iterbatch-1)*50,:,:,:) = dRSA;
        clear dRSA
        
    end
    
    % average over iterations within this batch
    dRSAperBatch(iterbatch,:,:,:) = squeeze(mean(dRSAperIter));
    clear dRSAperIter;
    
end% iterbatch loop

% average over batches
dRSAall = squeeze(mean(dRSAperBatch));
clear dRSAperBatch;

save(fn2save,'dRSAall','ROIsel');

end