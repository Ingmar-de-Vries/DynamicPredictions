% add events
clearvars

% Input files
rootdir = '\\XXX\ActionPrediction\data\brainstorm_database\ActionPrediction\data\';
cd(rootdir);
sub = 'XXXX';
runfolders = dir([rootdir sub filesep '*@raw*']);
sFiles = cell(1,length(runfolders));
for irun=1:length(runfolders)
    runfile = dir([rootdir sub filesep runfolders(irun).name filesep '*raw*']);
    sFiles{1,irun} = [sub filesep runfolders(irun).name filesep runfile.name];
end

% Start a new report
bst_report('Start', sFiles);

% Process: Read from channel
sFiles = bst_process('CallProcess', 'process_evt_read', sFiles, [], ...
    'stimchan',  'STI102', ...
    'trackmode', 1, ...  % Value: detect the changes of channel value
    'zero',      0);

% Process: Merge events
sFiles = bst_process('CallProcess', 'process_evt_merge', sFiles, [], ...
    'evtnames', '1,2,3,4,5,6,7,8,9,10,11,12,13,14', ...
    'newname',  'VideoOnset');

% Process: Merge events
sFiles = bst_process('CallProcess', 'process_evt_merge', sFiles, [], ...
    'evtnames', '101,102,103,104,105,106,107,108,109,110,111,112,113,114', ...
    'newname',  'CatchVideoOnset');

% Process: Merge events
sFiles = bst_process('CallProcess', 'process_evt_merge', sFiles, [], ...
    'evtnames', '201,202,203,204,205,206,207,208,209,210,211,212,213,214', ...
    'newname',  'CatchTestOnset');

% Process: Merge events
sFiles = bst_process('CallProcess', 'process_evt_merge', sFiles, [], ...
    'evtnames', '21,28', ...
    'newname',  'response');

% Process: Read from channel
sFiles = bst_process('CallProcess', 'process_evt_read', sFiles, [], ...
    'stimchan',  'STI102', ...
    'trackmode', 1, ...  % Value: detect the changes of channel value
    'zero',      0);

% Process: Delete events
sFiles = bst_process('CallProcess', 'process_evt_delete', sFiles, [], ...
    'eventname', '21,28');

% different EOG channels were used for different subjects, as the entrance for the EOG channels in the MSR was falling apart...
if strcmp(sub,'XXXX')   
    EOGchan = 'EEG062';
elseif strcmp(sub,'YYYY') 
    EOGchan = 'EEG063';
elseif strcmp(sub,'ZZZZ')
    EOGchan = 'EEG064';
end

% Process: Detect blink with custom settings
sFiles = bst_process('CallProcess', 'process_evt_detect', sFiles, [], ...
    'eventname',    'blink', ...
    'channelname',  EOGchan, ...
    'timewindow',   [], ...
    'bandpass',     [0.1, 15], ...
    'threshold',    2, ...
    'blanking',     0.5, ...
    'isnoisecheck', 1, ...
    'isclassify',   0);

% Save and display report
ReportFile = bst_report('Save', sFiles);
bst_report('Open', ReportFile);

