# Predictive neural representations of naturalistic dynamic input

This repository contains Matlab code accompanying the scientific article, currently available as pre-print, titled: "Predictive neural representations of naturalistic dynamic input", by Ingmar E.J. de Vries and Moritz F. Wurm. 

  -	The article is available at: https://doi.org/10.1101/2022.09.02.506366

For details regarding this experiment, stimuli and analysis code, please see methods section of the article. Please contact me for any further questions at i.e.j.de.vries@gmail.com

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
    - RDMs based on video data and eyetracker data are created outside of the main "DynamicPredictions_ERFdynRSA_sourceROI.m" script (see section "Run dRSA analysis" below). However, for the 6 RDMs based on kinematic marker data (absolute and relative body posture, motion, and acceleration), these are created inside the resampling iteration in the main script, because relative body posture, motion and acceleration can only be computed after realigning the 3-sec resampled sequences to each other. The precreated RDMs have size 14x14x250x250 (i.e., stim1 x stim2 x timestim1 x timestim2). 
    - In the "modelRDMs" subdirectory, you'll find the following scripts:
      - xxx

  - Run dRSA analysis
    - In the "dynamicRSA" subdirectory, you'll find the following analysis scripts:
      - "DynamicPredictions_pipeline.m" - the main analysis pipeline from which all other functions are called.
      - "cluster_shell.m" - used for sending analysis as parallel jobs to a computing cluster (e.g., with different subjects and ROIs in parallel).
      - "cluster_shell_simulations.m" - same for simulations.
      - "DynamicPredictions_defineSourceROIs.m" - create ROIs based on (combinations of) parcels of HCP atlas.
      - "DynamicPredictions_checkAtlases.m" - just sanity check that correct atlas and inversion kernel will be selected in main analysis
       
