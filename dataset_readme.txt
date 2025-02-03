MEG dataset used for the paper. Please notify me if any problems arise. A few notes about this dataset: 

  - The dataset can be found at DOI: https://doi.org/10.17605/OSF.IO/ZK42F
  
  - This is the real dataset used in the study. Note that for privacy reasons, the individual anatomical scan is not included. However, both the preprocessed data at sensor level and the inversion kernel to transform sensor level data to source level data are included, as well as the individual atlas definitions, which relate each vertex in the inversion kernel to a parcel in the atlas. The inversion kernel is computed using either the individual anatomical scan if available, or the template scan (see which one in the analysis logbook). The combination of the dataset and provided analysis code should allow one to reproduce the figures presented in the article. 
  
  - Data structure:
    - subjectinfo_and_logbook - excel file containing some subject information and a recording and analysis logbook. For privacy purposes non-essential information has been removed.
    - behaviour - Folder containing 1 file per run of output from the experiment script: “DynamicPredictions_MEGexperiment.m”.
    - eyetracker - Folder containing 1 .edf and 1 .asc file per run. The eyelink eyetracker stores data in .edf format, which I have reformatted to .asc. Additionally, the folder contains a single file of preprocessed eyetracker data in which runs are combined. 
    - MEGraw - Folder containing 1 file per run of MEG data in .fif format from the Elekta MEG scanner. This is the raw data, meaning that nothing has been done to it after recording. 
    - MEGmaxfiltered - Folder containing 1 file per run of MEG data in .fif format from which bad channels were manually removed (see which channels in third tab of "subjectinfo_and_logbook" file). Subsequently these files were MaxFiltered and spatially aligned (see article for details). 
    - MEGpreprocessedsensor - MEG dataset (runs combined) in matlab - Fieldtrip format after all preprocessing steps. This dataset is ready to be used in the main dRSA analysis pipeline "DynamicPredictions_pipeline.m". Data is in sensor space. 
    - inversionkernel - matlab file with inversion kernel to transform data from sensor to source space.
    - atlas - matlab file with HCP atlas which defines which vertex in the kernel belongs to which parcel from the HCP atlas. 
