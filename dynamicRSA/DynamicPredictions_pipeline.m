close all; clearvars; clc
addpath('\\cimec-storage5.unitn.it\MORWUR\Projects\INGMAR\ActionPrediction\code\neuralDecoding');

% Parameters for creating neural RDM
cfg.cluster=1;% 1 = run on cluster, 2 = run locally
cfg.SubVec = 1:22;% 
cfg.smoothingMSec = 20; % smoothing kernel across time, in msec. Applied before creating neural RDM.  
cfg.downsample = 100;% downsample neural data before creating neural RDM.  

% Parameters for dynamic RSA
cfg.smoothNeuralRDM=2;% smoothing of neural RDM. In samples, so with fs = 100 Hz, 20 msec smoothing = 2 samples
cfg.smoothModelRDM=2;% in samples
cfg.glmRSA=1;% which measure to use for computing similarity between neural and model RDM: 0 = correlation, 1 = principal component regression (PCR) weights
cfg.nPCAcomps = 75;% maximum amount of PCA components to regress out
cfg.nstim = 14;% number of unique video sequences
cfg.randshuff = [1000 2];% first value indicates amount of iterations, second value indicates maximum start time in sec of resampled sequence, i.e., 2 indicates the original 5 sec trials will be randomly resampled to 3 second trials on each iteration

% ROI definition
% DynamicPredictions_defineSourceROIs();

% check if correct atlases will be loaded for parcel definitions, not necessary to run if you see in the main script below that this is happening correctly
% cfg.atlas = 'HCP';% or Schaefer
% DynamicPredictions_checkAtlases(cfg);

%% ROI-based analysis
cfg.ROIVec = 1:6;% V1, V2, V3V4, LOTC, aIPL, PMv
cfg.dynRDMnames = {'graysmooth','opticalFlowFBmag','opticalFlowFB','posabs','posrel','velabs','velrel','accabs','accrel','eyePos'};% models to include in regression
cfg.models2test = [1:7 10];% only test these from the above models, but regress out all 9 others when testing a certain model
cfg.atlas = 'HCP';% HCP or other. WARNING: CURRENTLY ONLY HCP IS TESTED FOR ROI-BASED ANALYSIS

% send analysis to cluster, where it will run subjects and ROIs in parallel: 
script2run = 'DynamicPredictions_ERFdynamicRSA_ROIsource';
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
script2run = 'DynamicPredictions_STATS_ERFdynamicRSA_ROIsource';
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

