# Predictive neural representations of naturalistic dynamic input

This repository contains Matlab code accompanying the scientific article, currently available as pre-print, titled: "Predictive neural representations of naturalistic dynamic input", by Ingmar E.J. de Vries and Moritz F. Wurm. 

  -	The article is available at: https://doi.org/10.1101/2022.09.02.506366

This article presents a new dynamic extension to the influential representational similarity analysis (RSA) approach, which allows investigating when naturalistic dynamic stimuli are represented in different parts of the brain, at different hierarchical levels of processing (e.g., low-level visual, body posture or motion of a ballet dancer). The code in this repository allows for exact replication of the dRSA pipeline as presented in the article, but is also meant as inspiration for people interested in using dRSA to answer their own research questions. In the near future I plan to generalize the code a bit more (e.g., less hard coding and more soft coding inside functions, more options, etc.) and create a separate repository dedicated to the dRSA method itself. In principle, dRSA can be applied to many different sensory modalities and contexts (e.g., naturalistic sound scenes, music, language), and on any signal with high enough temporal resolution (M/EEG, ECoG, eyetracking, etc.).  

For details regarding this experiment, stimuli and analysis code, please see methods section of the article. Please contact me for any further questions at i.e.j.de.vries@gmail.com

Note that this custom-written code uses several functions from the Brainstorm, Fieldtrip and EEGLAB toolboxes. 

The code is structured as follows:

  -	Experiment
    - In the “experiment” subdirectory, you will find the Matlab script “DynamicPredictions_MEGexperiment.m” 
    - You need Psychophysics Toolbox Version 3 (PTB-3) to run this experiment. 
    - In the subdirectory “experiment/stimuli”, you will find the 14 unique 5-second-long ballet dancing videos used in the experiment, plus the corresponding and temporally aligned 3D kinematic marker locations at 100 Hz, stored in Matlab matrices. 
    - The experiment script makes use of the following helper scripts also present in the experiment directory:
      - “angle2pix.m” – transform degrees of visual angle to pixels on screen
      - “CreateCatchTrials.m” – create pool of catch trials to pick from.

  -	Pre-processing of MEG and eyetracking data
    - Pre-processing of eyetracking data was done using custom written script "DynamicPredictions_eyetracking_asc2ft.m", which takes raw Eyelink data in asc format as input, and gives pre-processed eyetracking data in Fieldtrip format as output. This is subsequently used to create an eyetracker RDM per individual subject. 
    - Pre-processing of MEG data was done using the Brainstorm toolbox version 3 using GUI operations, which were transformed into Matlab scripts where possible. The main distinction is with preparing the individual-subject anatomy (which was done using the GUI), and preprocessing of the MEG data (which was done with scripts). The following scripts are located in the subdirectory “preprocessing”:
      - "DynamicPredictions_MEGpp1_PSDcheck.m" - initial sanity/quality check of powerspectra. Can be skipped as this will be done after filters anyway. 
      - "DynamicPredictions_MEGpp2_addEvents.m" - read events from trigger channel and give appropriate names.
      - "DynamicPredictions_MEGpp3_addPTBevents.m" - not necessary to use this, check comments in script for details.
      - "DynamicPredictions_MEGpp4_checkVidOnset.m" - only sometimes necessary, i.e., sometimes triggers were erroneously stored double. If that's the case, this script helps finding those duplicates so they can be removed manually in Brainstorm GUI. But only happened in rare cases. 
      - "DynamicPredictions_MEGpp5_filters.m" - notch filter, downsample, and create powerspectra for sanity check.
      - "DynamicPredictions_MEGpp6_ICA.m" - run ICA for ocular and cardiac artifacts, separately for magneto- and gradiometers.
      - "DynamicPredictions_MEGpp7_detectArtifacts.m" - I skipped this step, because I opted for manual artifact detection after epoching. 
      - "DynamicPredictions_MEGpp8_epoch_singletrialDCcorrection.m" - epoch and single-trial baseline correction.
      - "DynamicPredictions_MEGpp9_export2FT.m" - export from Brainstorm to Fieldtrip format.
      - "DynamicPredictions_MEGpp10_realign2photodiode.m" - realign single trials to photodiode. This script is called from MEGpp9. 
      - "DynamicPredictions_MEGppSource1_computeInversionKernel.m" - Apply minimum norm estimation (MNE) and store resulting inversion kernel to transform sensor level data to source level data outside of Brainstorm (which I do in the main dynamic RSA analysis script). 
      - "DynamicPredictions_storeManualBadTrials.m" - store manually detected bad trials, see script for comments

  -	Create model RDMs
    - RDMs based on video data and eyetracker data are created outside of the main "DynamicPredictions_ERFdynRSA_sourceROI.m" script (see section "Run dRSA analysis" below). However, for the 6 RDMs based on kinematic marker data (absolute and relative body posture, motion, and acceleration), these are created inside the resampling iteration in the main script, because relative body posture, motion and acceleration can only be computed after realigning the 3-sec resampled sequences to each other. The pre-created RDMs have size 14x14x250x250 (i.e., stim1 x stim2 x timestim1 x timestim2) at 50 Hz. 
    - In the "modelRDMs" subdirectory, you'll find the following scripts:
      - xxx

  - Run dRSA analysis, statistics, and plotting
    - In the "dynamicRSA" subdirectory, you'll find the following scripts:
      - "DynamicPredictions_pipeline.m" - the main analysis pipeline from which all other functions are called.
      - "cluster_shell.m" - used for sending analysis as parallel jobs to a computing cluster (e.g., with different subjects and ROIs in parallel).
      - "cluster_shell_simulations.m" - same for simulations.
      - "DynamicPredictions_defineSourceROIs.m" - create ROIs based on (combinations of) parcels of HCP atlas.
      - "DynamicPredictions_checkAtlases.m" - just sanity check that correct atlas and inversion kernel will be selected in main analysis
      - "DynamicPredictions_ERFdynamicRSA_ROIsource.m" - main analysis script, which is called from "DynamicPredictions_pipeline.m"
      - "procrustes_constrain_rotationZaxis_IdV.m" - modified version of Matlab's procrustes.m, which now constrains rotation to vertical (Z) axis, because that is how we define viewpoint invariant (i.e., relative) body posture, motion and acceleration. Note that my modified version is correct, but currently very time inefficient, effectively more than doubling the total computation time. This was a later modification and I'm sure I can find a much faster implementation. Feel free to have a look in the script and suggest a faster implementation! 
      - "regressionBorderPerModel_smRDM20msec.mat" - file containing regression borders used to regress out model itself to attenuate effects of model autocorrelation. These borders are determined by the simulations (see methods section in article and explanation in "DynamicPredictions_ERFdynamicRSA_ROIsource.m" for details). 
      - "DynamicPredictions_STATS_ERFdynamicRSA_ROIsource.m" - run statistics on ROI-based dRSA results, and compute peak latency and representational spread (RS) index. This function is called from main script "DynamicPredictions_pipeline.m". 
      - "modelautocorr_slopes.mat" - file containing dRSA curves resulting from PCR on simulated data. This is used to compute the representational spread (RS) index (see methods section in article and explanation in "DynamicPredictions_STATS_ERFdynamicRSA_ROIsource.m" for details).
      - "DynamicPredictions_runFTstats.m" - shell around Fieldtrip functions for running cluster-based permutation tests on 2D dRSA matrix or on averaged dRSA lag-plot. This function is called from "DynamicPredictions_STATS_ERFdynamicRSA_ROIsource.m". See scripts for details. 
      - "DynamicPredictions_PLOTS_ERFdynamicRSA_ROIsource.m" - plot ROI-based results, in article: figure 2a, 3, and S5.
      - "brewermap.m" - creates nice colormaps that are colorblind friendly. Not my code, for all colormaps and source code see: https://colorbrewer2.org/
      - "boundedline.m" - creates nice shading around lines, e.g., with a measure of distribution across subjects (here standard error). Not my code, for source code see https://github.com/kakearney/boundedline-pkg
      - "DynamicPredictions_ERFdynamicRSA_searchlight.m" - searchlight analysis, which is called from "DynamicPredictions_pipeline.m"

- Run dRSA simulations, and plotting
  - In the "simulations" subdirectory, you'll find the following scripts:
    - 'XXX'
