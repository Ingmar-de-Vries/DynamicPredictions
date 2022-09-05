close all; clearvars; clc
addpath('\\cimec-storage5.unitn.it\MORWUR\Projects\INGMAR\ActionPrediction\code\neuralDecoding');

% Parameters for creating RDMs
cfg.cluster=1;% 1 = run on cluster, 2 = run locally
cfg.SubVec = 1:22;%
cfg.radius = 0; % in Time bins
cfg.crossval = 0;%only for correlation/Euclidean, classifiers are always cross-validated
cfg.classifier = 3; % lda, svm, corr, euclidean
cfg.averaging = 0; % number of Trials averaged before classification
cfg.resampling = 0; % number of resampled averages
cfg.chunksize = 0; % number of chunks
cfg.smoothingMSec = 20; % smoothing kernel across time, in msec 
cfg.downsample = 100;
cfg.normalize = 0;
cfg.TestVec = 1; % optionally different classification schemes, for us doesn't make sense to do anything alse than pairwise classification of 14 videos
cfg.TestName = 'multiclassVid';
cfg.type =   'classifyMat';
cfg.chan_type= 'all'; % 'all';% 'meg_axial'; %   'meg_planar'; %
cfg.scale_sensors = 0;
cfg.MNN = 1;% multivariate noise normalisation (Guggenmos et al. 2018 NeuroImage)

% Parameters for dynamic RSA
cfg.smoothNeuralRDM=2;% smoothing of neural RDM. In samples, so with fs = 100 Hz, 100 msec smoothing = 10 samples
cfg.smoothModelRDM=2;% in samples
cfg.randomModelPERM = 0;% 1 = randomize stimuli neural RDM, 2 = randomize stimuli model RDM, 3 = randomize time neural RDM, 4 = randomize time model RDM

cfg.glmRSA=2;% 0 = correlation, 3 = regularized regression, now in combination with PCA, 4 = PCA and non-regularized regression, 5 = partial least-squares regression (PLSR), 6 = principal component regression (PCR) 
cfg.nPCAcomps = 75;% maximum amount of PCA components to regress out
cfg.nstim = 14;% 14 = run the standard analysis on 14 stimuli, 28 = split the stimuli in 2 to create much larger RDMs (i.e., 378 instead of 91)
cfg.randshuff = [250 2];% .94 = largest difference in onset between sequences of any move, 1.4 = average onset second move, 1.86 = max onset second move (so 1st and 2nd moves can even be synchronized)

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOURCE LEVEL ANALYSES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ROI definition
ActionPrediction_defineSourceROIs();

% simple script to check if correct atlases will be loaded for parcel definitions
cfg.atlas = 'HCP';
ActionPrediction_checkAtlases(cfg);

%% ROI-based analysis
cfg.ROIVec = 1:6;%[2 3 7 8 9 10];%V1, V2, V3V4, LOTC, aIPL, PMv, post dors stream, ant dors stream, PMd, pDLPFC
cfg.dynRDMnames = {'graysmooth','opticalFlowFBmag','opticalFlowFB','posabs','posrel','velabs','velrel','accabs','accrel','eyePos'};%,'COMpos','COMvel','COMacc'};
cfg.models2test = [1:7 10];% only test these from the above models, but regress out all 10 others when testing a certain model
cfg.atlas = 'HCP';% HCP, or Schaefer_200, or other 
cfg.splitLR = 0;% separate for left and right hemisphere, if 0 combine hemispheres

script2run = 'ActionPrediction_ERPdynRSA_sourceROIs';
cluster_shell(cfg,script2run);

% statistics
cfg.subnum = length(cfg.SubVec);
cfg.sub4stat = cfg.SubVec;
cfg.SubVec = 1;%only so cluster_shell sends job once rather than for each subject
cfg.plotrandmod = 0;
cfg.fisherz = 1;%fisher transform for individual subject correlation values before stats and plotting
cfg.jackknife = 1;% use jackknife method for estimating peak latency and sustained representation index (SRI)
cfg.timeXtime = 0;% whether to run stats on the full time-by-time dRSA results, or only the line plots after averaging over stimulus time
cfg.pthresh = 0.05;% threshold for cluster size test
cfg.pthreshclust = 0.01;% single sample tests to determine which samples are included in cluster
cfg.ROIVec = 1:6;%V1, V2, V3V4, LOTC, aIPL, PMv, post dors stream, ant dors stream, PMd, pDLPFC
cfg.randshuff(1) = 1000;
script2run = 'ActionPrediction_STATS_ERPdynRSA_sourceROIs';
cluster_shell(cfg,script2run);

% plot ROI analysis
cfg.SubVec = 1:22;
cfg.fisherz = 1;%fisher transform for individual subject correlation values before stats
cfg.normalize = 0;% normalize individual subject dRSA line plots before averaging 
cfg.jackknife = 1;% use jackknife method for estimating peak latency and sustained representation index (SRI)
cfg = rmfield(cfg,'cluster');
cfg.ROIVec = 1:6;%V1, V2, V3V4, LOTC, aIPL, PMv, post dors stream, ant dors stream, PMd, pDLPFC
cfg.randshuff(1) = 1000;

ActionPrediction_PLOT_ERPdynRSA_sourceROIs(cfg);

%% semi-searchlight analysis
cfg.ROIVec = 1;% just call it 1 so cluster_shell sends it for only 1 ROI
cfg.atlas = 'Schaefer2018_400';%'HCP';%'Schaefer2018_400';% use Schaeffer 200 or 400 (because ROIs have roughly same size)
cfg.peaks2test = [1 .11 ; 2 .07 ; 3 -.09 ; 3 .2 ; 4 .18 ; 5 .07 ; 6 -.17 ; 7 -.44 ; 10 .09];% only compute dRSA at one or two lags (in sec) to save time. Lag is determined by finding peak in ROI based analysis
cfg.randshuff = [1000 2];
script2run = 'ActionPrediction_ERPdynRSA_source_semisearchlight';
cluster_shell(cfg,script2run);

% stats and plot semi-searchlight in Brainstorm
cfg.fisherz = 1;%fisher transform for individual subject correlation values before stats
cfg.sub4stat = cfg.SubVec;
cfg.SubVec = 1;%only so cluster_shell sends job once rather than for each subject
cfg.pvalFDR = 0.05;
cfg = rmfield(cfg,'cluster');
ActionPrediction_STATS_ERPdynRSA_source_semisearchlight(cfg);
ActionPrediction_PLOT_ERPdynRSA_source_semisearchlight(cfg);

%% simulations
cfg.lag = 0;%.12;% in sec
cfg.randshuff(1) = 150;
cfg.SubVec = [1:10 99];% these are actually models to implant, just call it SubVec so it works with cluster_shell.m in parallelizing on the cluster
cfg.ROIVec = 7;% only send once to cluster for single randomization, or 1:21 if including eye data so it's send for each subject's eye data (which we can average over afterwards)
cfg.iterationsPERbatch = 30;
cfg.iterbatches = cfg.randshuff(1)./cfg.iterationsPERbatch;
script2run = 'ActionPrediction_ERPdynRSA_sourceROIs_simulation';
cluster_shell_simulations(cfg,script2run);

% simulate combinations of models
cfg.SubVec = 123;% these are actually models to implant, just call it SubVec so it works with cluster_shell.m in parallelizing on the cluster
script2run = 'ActionPrediction_ERPdynRSA_sourceROIs_combisim';
cluster_shell_simulations(cfg,script2run);

% plot main simulations
cfg = rmfield(cfg,'cluster');
ActionPrediction_PLOT_ERPdynRSA_sourceROIs_simulation(cfg);

% plot extra simulations for paper (i.e., 1. combining models, 2. different levels of noise, 3. distributed representation)
cfg = rmfield(cfg,'cluster');
ActionPrediction_PLOT_ERPdynRSA_sourceROIs_simulationExtras(cfg);

