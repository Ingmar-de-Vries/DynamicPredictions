# Predictive neural representations of naturalistic dynamic input

UNDER CONSTRUCTION

This repository contains Matlab code accompanying the scientific article, currently available as pre-print, titled: "Predictive neural representations of naturalistic dynamic input", by Ingmar de Vries and Moritz Wurm. 

  -	The article is available at: xxx

For details regarding this experiment, stimuli and analysis code, please see methods section of the article. Please contact me for any further questions at i.e.j.de.vries@gmail.com

The code is structured as follows:

  -	Experiment
    o	In the “experiment” subdirectory, you will find the Matlab script “DynamicPredictions_MEGexperiment.m” 
    o	In the subdirectory “experiment/stimuli”, you will find the 14 unique 5-second-long ballet dancing videos used in the experiment, plus the corresponding and temporally aligned 3D kinematic marker locations at 100 Hz, stored in Matlab matrices. 
    o	The experiment script makes use of the following helper scripts or files also present in the experiment directory:
      	“angle2pix.m” – transform degrees of visual angle to pixels on screen
      	“CreateCatchTrials.m” – create pool of catch trials to pick from
    o	Additionally, you’ll need to have Psychophysics Toolbox Version 3 (PTB-3) downloaded and installed to run this experiment. 

  -	Pre-processing of MEG data
    o	Pre-processing of MEG data was done using the Brainstorm toolbox version 3 (201019), using GUI operations, which were transformed into Matlab scripts where possible. These are located in the subdirectory “preprocessing”.

  -	Create model RDMs

-	Run dRSA analysis
