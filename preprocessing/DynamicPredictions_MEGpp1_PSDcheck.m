% Run initial PSD for first data quality check. For example to discover
% noisy channel not excluded before MaxFilter. If there is a very noisy
% channel, I redid MaxFilter for that run now excluding that channel
clearvars; 

% Input files
rootdir = '\\XXX\ActionPrediction\data\brainstorm_database\ActionPrediction\data\';
cd(rootdir);
sub = 'XXXX';% fill in 4-letter part of subject code
runfolders = dir([rootdir sub filesep '*@raw*']);
sFiles = cell(1,length(runfolders));
for irun=1:length(runfolders)
    runfile = dir([rootdir sub filesep runfolders(irun).name filesep '*raw*']);
    sFiles{1,irun} = [sub filesep runfolders(irun).name filesep runfile.name];
end

% Start a new report
bst_report('Start', sFiles);

% Process: Power spectrum density (Welch)
sFiles = bst_process('CallProcess', 'process_psd', sFiles, [], ...
    'timewindow',  [], ...
    'win_length',  1, ...
    'win_overlap', 50, ...
    'sensortypes', 'MEG', ...
    'win_std',     0, ...
    'edit',        struct(...
         'Comment',         'Avg,Power', ...
         'TimeBands',       [], ...
         'Freqs',           [], ...
         'ClusterFuncTime', 'none', ...
         'Measure',         'power', ...
         'Output',          'all', ...
         'SaveKernel',      0));

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);
