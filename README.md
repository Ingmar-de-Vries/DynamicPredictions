# Predictive neural representations of naturalistic dynamic input

This repository contains Matlab code accompanying the scientific article, currently available as pre-print, titled: "Predictive neural representations of naturalistic dynamic input", by Ingmar E.J. de Vries and Moritz F. Wurm. 

  -	The article is available at: xxx

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
      - "DynamicPredictions_MEGpp2_addEvents.m" - read events from trigger channel and give appropriate names
      - "xxx"

  -	Create model RDMs

  - Run dRSA analysis
    - In the "dynamicRSA" subdirectory, you'll find the following analysis scripts:
      - "cluster_shell.m" - used for sending analysis as parallel jobs to a computing cluster (e.g., with different subjects and ROIs in parallel).
      - "cluster_shell_simulations.m" - same for simulations.
