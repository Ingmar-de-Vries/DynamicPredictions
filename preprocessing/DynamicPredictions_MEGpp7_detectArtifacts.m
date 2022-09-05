% Detect additional artifacts: 1-7 Hz high amplitude for movements
% (incl. eye movements), 40-240 Hz high amplitude for muscle activity
clearvars; 

% Input files
rootdir = '\\XXX\ActionPrediction\data\brainstorm_database\ActionPrediction\data\';
cd(rootdir);
sub = 'XXXX';
runfolders = dir([rootdir sub filesep '*resample*']);
sFiles = cell(1,length(runfolders));
for irun=1:length(runfolders)
    runfile = dir([rootdir sub filesep runfolders(irun).name filesep '*raw*']);
    sFiles{1,irun} = [sub filesep runfolders(irun).name filesep runfile.name];
end

% Start a new report
bst_report('Start', sFiles);

% Process: Detect other artifacts
sFiles = bst_process('CallProcess', 'process_evt_detect_badsegment', sFiles, [], ...
    'timewindow',  [], ...
    'sensortypes', 'MEG', ...
    'threshold',   4, ...
    'isLowFreq',   1, ...
    'isHighFreq',  1);

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);

