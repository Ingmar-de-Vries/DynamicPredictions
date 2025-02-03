# Predictive neural representations of naturalistic dynamic input

DOI code: 10.5281/zenodo.7941212

This repository contains Matlab code accompanying the scientific article available at 

De Vries, I.E.J., Wurm, M.F. Predictive neural representations of naturalistic dynamic input. Nat Commun 14, 3858 (2023). https://doi.org/10.1038/s41467-023-39355-y

If you use the data on OSF or the code here, please cite the above article.

This article presents a new dynamic extension to the influential representational similarity analysis (RSA) approach, which allows investigating when naturalistic dynamic stimuli are represented in different parts of the brain, at different hierarchical levels of processing (e.g., low-level visual, body posture or motion of a ballet dancer). The code in this repository allows for exact replication of the dRSA pipeline as presented in the article, but is also meant as inspiration for people interested in using dRSA to answer their own research questions. In the future I will generalize the code a bit more (e.g., less hard coding inside functions, more options, etc.) and create a separate repository dedicated to the dRSA method itself. In principle, dRSA can be applied to many different sensory modalities and contexts (e.g., naturalistic sound scenes, music, language), and on any signal with high enough temporal resolution (M/EEG, ECoG, eyetracking, etc.). Additionally, it should be straightforward to implement different dissimilarity measures for the neural and model RDMs, and a different similarity measure for the dRSA (i.e., here the principal component regression approach). We have found this similarity measure to be most effective for the current experiment, as tested with the simulations, but we have not extensively tested different dissimilarity measures.  

For details regarding this experiment, stimuli and analysis code, please see methods section of the article. Please contact me for any further questions at i.e.j.de.vries@gmail.com

Note that the larger data files belonging to this repository are stored on a public OSF repository (DOI: [10.17605/OSF.IO/ZK42F](https://doi.org/10.17605/OSF.IO/ZK42F); or look for Ingmar de Vries - DynamicPredictions). The OSF repository includes: 
  - MEG data. See more information in the "dataset_readme.txt" file included in this repository. 
  - A zipped folder titled 'Source Data'. This folder includes the data of all 22 subjects for all main results (Fig. 2a and b, Fig. 3a and c, and Fig. S1).
  - The 9 model RDMs used in the reported study.

Note that this custom-written code uses several functions from the Brainstorm (tested version: 3), Fieldtrip (tested version: 20191113) and EEGLAB (tested version: 2019.1) toolboxes, and was written and tested in Matlab 2020a.  

The code in this GitHub repository is structured as follows:

  -	Experiment
    - In the “experiment” subdirectory, you will find the experiment script “DynamicPredictions_MEGexperiment.m” 
    - You need Psychophysics Toolbox Version 3 (PTB-3) to run this experiment. 
    - This experiment can in principle be run as a behaviour-only experiment, include eyetracking, or include eyetracking and MEG. However, I have only tested the latest version of this experiment in the MEG lab at CIMeC, using Matlab 2012b. You might need to make minor adjustments for your setup. 
    - In the subdirectory “experiment/stimuli”, you will find the 14 unique 5-second-long ballet dancing videos used in the experiment, plus the corresponding and temporally aligned 3D kinematic marker locations at 100 Hz, stored in Matlab matrices. 
    - The experiment script makes use of the following helper scripts also present in the experiment directory:
      - “angle2pix.m” – transform degrees of visual angle to pixels on screen
      - “CreateCatchTrials.m” – create pool of catch trials that the experiment script randomly picks from on each run. 
    - In this subdirectory you will additionally find the script "DynamicPredictions_exampleFigure_videoANDstickfigure.m", which plots some example video frames and the respective stick figures based on the kinematic markers. It was used for creating Figure 1a and b in the article. This script also contains information about where each of the 13 kinematic markers were located on the ballet dancer's body. 

  - Analysis and plotting of behavioural results
    - In the "behaviour" subdirectory, you will find the following script:
      - "DynamicPredictions_AnalysisBehaviourCatchTrials.m" - behavioural analysis and plotting of Figure 4c.
      - "plotSpread.m" - can be found in the "plotSpread" subdirectory, together with some helper functions, and is used for the dotcloud figure.

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
      - "DynamicPredictions_MEGppSource1_computeInversionKernel.m" - apply minimum norm estimation (MNE) and store resulting inversion kernel to transform sensor level data to source level data outside of Brainstorm (which I do in the main dynamic RSA analysis script). 
      - "DynamicPredictions_storeManualBadTrials.m" - store manually detected bad trials, see script for comments

  -	Create model RDMs
    - The 9 model RDMs themselves can be found in the OSF repository. 
    - RDMs based on video data, kinematic marker data, and eyetracker data are created outside of the main "DynamicPredictions_ERFdynamicRSA_ROIsource.m" script (see section "Run dRSA analysis" below), and if necessary up- or downsampled to 100 Hz. The pre-created RDMs therefore have size 14x14x500x500 (i.e., stim1 x stim2 x timestim1 x timestim2). This is done to save computation time in the dRSA pipeline, i.e., the RDMs only need to be loaded in, not computed each time the dRSA script is run. 
    - In the "modelRDMs" subdirectory, you'll find the following scripts:
      - "DynamicPredictions_DynamicModelRDMs_eyeTracker.m" - create dynamic RDM of individual subject eyetracking data.
      - "DynamicPredictions_DynamicModelRDMs_pixelwise.m" - create dynamic RDM of smoothed grayscale pixelwise luminance values. 
      - "DynamicPredictions_video2vector.m" - create smoothed grayscale vector representation of videos. Called from "DynamicPredictions_DynamicModelRDMs_pixelwise.m".
      - "DynamicPredictions_DynamicModelRDMs_opticalflow.m" - create dynamic RDM of optical flow vectors.
      - "DynamicPredictions_video2opticalflow.m" - create optical flow vector representation of videos. Called from "DynamicPredictions_DynamicModelRDMs_opticalflow.m".
      - "DynamicPredictions_DynamicModelRDMs_kinematic.m" - create 6 dynamic RDMs of kinematic marker data, i.e., view-dependent and view-invariant posture, motion and acceleration.
      - "procrustes_constrain_rotationZaxis_IdV.m" - modified version of Matlab's procrustes.m, which now constrains rotation to vertical (Z) axis, because that is how we define viewpoint invariant body posture, motion and acceleration. Note that my modified version is correct, but currently very time inefficient, effectively more than doubling the total computation time. This was a later modification and I'm sure I can find a much faster implementation. Feel free to have a look in the script and suggest a faster implementation! Called from "DynamicPredictions_DynamicModelRDMs_kinematic.m".
      - "DynamicPredictions_exampleFigureModels.m" - plots illustrations of the different models for a single frame of 2 videos. It was used for creating Figure 5 in the article. This script also contains information about where each of the 13 kinematic markers were located on the ballet dancer's body.

  - Run dRSA analysis, statistics, and plotting
    - In the "dynamicRSA" subdirectory, you'll find the following scripts:
      - "DynamicPredictions_pipeline.m" - the main analysis pipeline from which all other functions are called.
      - "cluster_shell.m" - used for sending analysis as parallel jobs to a computing cluster (e.g., with different subjects and ROIs in parallel).
      - "DynamicPredictions_defineSourceROIs.m" - create ROIs based on (combinations of) parcels of HCP atlas.
      - "DynamicPredictions_checkAtlases.m" - just sanity check that correct atlas and inversion kernel will be selected in main analysis
      - "DynamicPredictions_RUN_ERFdynamicRSA_ROIsource.m" - main analysis script, which is called from "DynamicPredictions_pipeline.m"
      - "regressionBorderPerModel_smRDM30msec.mat" - file containing regression borders used to regress out model itself to attenuate effects of model autocorrelation. These borders are determined by the simulations (see methods section in article and explanation in "DynamicPredictions_ERFdynamicRSA_ROIsource.m" for details). 
      - "DynamicPredictions_STATS_ERFdynamicRSA_ROIsource.m" - run statistics on ROI-based dRSA results, and compute peak latency and representational spread (RS) index. This function is called from main script "DynamicPredictions_pipeline.m". 
      - "modelautocorr_slopes.mat" - file containing dRSA curves resulting from PCR on simulated data. This is used to compute the representational spread (RS) index (see methods section in article and explanation in "DynamicPredictions_STATS_ERFdynamicRSA_ROIsource.m" for details).
      - "DynamicPredictions_runFTstats.m" - shell around Fieldtrip functions for running cluster-based permutation tests on 2D dRSA matrix or on averaged dRSA lag-plot. This function is called from "DynamicPredictions_STATS_ERFdynamicRSA_ROIsource.m". See scripts for details. 
      - "DynamicPredictions_PLOTS_ERFdynamicRSA_ROIsource.m" - plot ROI-based results, in article: figure 2a, 3, and S1.
      - "brewermap.m" - creates nice colormaps that are colorblind friendly. Not my code, for all colormaps and source code see: https://colorbrewer2.org/
      - "boundedline.m" - creates nice shading around lines, e.g., with a measure of distribution across subjects (here standard error). Not my code, for source code see https://github.com/kakearney/boundedline-pkg
      - "DynamicPredictions_RUN_ERFdynamicRSA_searchlight.m" - searchlight analysis, which is called from "DynamicPredictions_pipeline.m"
      - "DynamicPredictions_STATS_ERFdynamicRSA_searchlight.m" - statistics on searchlight analysis. 
      - "fdr_bh.m" - FDR correction for statistics on searchlight analysis. This function is called in "DynamicPredictions_STATS_ERFdynamicRSA_searchlight.m"
      - "DynamicPredictions_PLOTS_ERFdynamicRSA_searchlight.m" - plot searchlight results, in article: figure 2b. Note that this is done partly using the Brainstorm GUI to create pretty cortical map figures. Go through this script line-by-line and read the comments if you want to create similar figures. 

- Run dRSA simulations, and plotting
  - Note that the simulations are also called from the main "DynamicPredictions_pipeline.m" script. In the "simulations" subdirectory, you'll find the following scripts:
    - "cluster_shell_simulations.m" - used for sending analysis as parallel jobs to a computing cluster (e.g., with different subjects and ROIs in parallel).
    - "DynamicPredictions_RUN_ERFdynamicRSA_simulations.m" - run simulations.
    - "DynamicPredictions_PLOTS_ERFdynamicRSA_simulations.m" - plot simulations. 
