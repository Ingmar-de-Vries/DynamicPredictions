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
% note that this analysis may take up a lot of time! It is worth first running it on a smaller number of iterations, and perhaps see where the
% analysis can be made more time-efficient for your purposes. Tip: On several occasions I ran the analysis e.g., 2 times in parallel with 2 batches of 500 iterations,
% and average over those afterwards. Do make sure to give different names to your output folder or files so you don't overwrite those two batches. 

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

% send analysis to cluster, where it will run subjects and ROIs in parallel
script2run = 'DynamicPredictions_RUN_ERFdynamicRSA_ROIsource';
cluster_shell(cfg,script2run);

% statistics
cfg.subnum = length(cfg.SubVec);
cfg.sub4stat = cfg.SubVec;
cfg.SubVec = 1;%only so cluster_shell sends job once for all subjects, rather than for each subject
cfg.fisherz = 1;%fisher transform for individual subject correlation values before stats and plotting
cfg.jackknife = 1;% use jackknife method for estimating peak latency and representational spread (RS)
cfg.timeXtime = 0;% whether to run stats on the full time-by-time dRSA results, or only the line plots after averaging over stimulus time
cfg.pthresh = 0.05;% or 0.01; threshold for cluster size test
cfg.pthreshclust = 0.01;% or 0.001; single sample tests to determine which samples are included in cluster
cfg.ROIVec = 1:6;%V1, V2, V3V4, LOTC, aIPL, PMv

% send statistical analysis to clsuter, where it will run ROIs in parallel
script2run = 'DynamicPredictions_STATS_ERFdynamicRSA_ROIsource';
cluster_shell(cfg,script2run);

% plot ROI analysis
cfg.SubVec = 1:22;
cfg.fisherz = 1;%fisher transform for individual subject correlation values before stats
cfg.jackknife = 1;% use jackknife method for estimating peak latency and representational spread (RS)
cfg = rmfield(cfg,'cluster');% we do this locally
cfg.ROIVec = 1:6;%V1, V2, V3V4, LOTC, aIPL, PMv

DynamicPredictions_PLOTS_ERFdynamicRSA_ROIsource(cfg);

%% searchlight analysis
cfg.ROIVec = 1;% just call it 1 so cluster_shell sends it only once (at least for ROIs, it will still be send in parallel for each subject)
cfg.atlas = 'Schaefer2018_400';% we used Schaefer2018_400 because ROIs have roughly same size. WARNING: CURRENTLY ONLY Schaefer2018_400 IS TESTED FOR SEARCHLIGHT ANALYSIS
cfg.peaks2test = [1 .11 ; 2 .07 ; 3 -.09 ; 3 .2 ; 4 .18 ; 5 .07 ; 6 -.17 ; 7 -.44 ; 10 .09];% only compute dRSA at peak to save computation time. Lag is determined by finding peak in ROI based analysis
cfg.randshuff = [1000 2];
script2run = 'DynamicPredictions_RUN_ERFdynamicRSA_searchlight';
cluster_shell(cfg,script2run);

% stats and plot searchlight (plotting done using Brainstorm GUI, see script and article for details)
cfg.fisherz = 1;%fisher transform for individual subject correlation values before stats
cfg.sub4stat = cfg.SubVec;
cfg.SubVec = 1;%only so cluster_shell sends job once rather than for each subject
cfg = rmfield(cfg,'cluster');
DynamicPredictions_STATS_ERFdynamicRSA_searchlight(cfg);
DynamicPredictions_PLOTS_ERFdynamicRSA_searchlight(cfg);

%% simulations
cfg.lag = 0;% in sec, at which lag to 'implant' the simulated model RDM
cfg.randshuff(1) = 150;% for simulations fewer iterations are needed
cfg.SimVec = [1:10 99];% models to simulate, 99 = only random neural data, no simulated model RDM
cfg.ROIVec = 7;% only send once to cluster for single randomization
cfg.iterationsPERbatch = 30;% for simulation iterations are send to the cluster in batches in parallel
cfg.iterbatches = cfg.randshuff(1)./cfg.iterationsPERbatch;
cfg.randweight = 0;% weight of random neural RDM (noise) to which simulated model will be added, weight of zero means no random neural RDM
script2run = 'DynamicPredictions_RUN_ERFdynamicRSA_simulations';
cluster_shell_simulations(cfg,script2run);

% plot main simulations
cfg = rmfield(cfg,'cluster');
DynamicPredictions_PLOTS_ERFdynamicRSA_simulations(cfg);

