%% import data into database and epoch
clearvars;

% Input files
rootdir = '\\XXX\ActionPrediction\data\brainstorm_database\ActionPrediction\data\';
cd(rootdir);
subfilz = dir('*');

for isub = 1:length(subfilz)
    sub = subfilz(isub).name;
    
    runfolders = dir([rootdir sub filesep '*resample']);
    runfolders = runfolders(contains(extractfield(runfolders,'name'),'@raw'));
    sFiles = cell(1,length(runfolders));
    for irun=1:length(runfolders)
        runfile = dir([rootdir sub filesep runfolders(irun).name filesep '*raw*']);
        sFiles{1,irun} = [sub filesep runfolders(irun).name filesep runfile.name];
    end
    
    % Start a new report
    bst_report('Start', sFiles);
    
    % Process: Import MEG/EEG: Events
    sFiles = bst_process('CallProcess', 'process_import_data_event', sFiles, [], ...
        'subjectname', sub, ...
        'condition',   '', ...
        'eventname',   'VideoOnset,CatchVideoOnset', ...%
        'timewindow',  [], ...
        'epochtime',   [-2.5, 7.5], ...
        'createcond',  0, ...
        'ignoreshort', 1, ...
        'usectfcomp',  0, ...
        'usessp',      1, ...
        'freq',        [], ...
        'baseline',    []);
    
    % Save and display report
    ReportFile = bst_report('Save', sFiles);
    bst_report('Open', ReportFile);
end

%% single-trial baseline correction (DC offset)
clearvars;

% Input files
rootdir = '\\XXX\ActionPrediction\data\brainstorm_database\ActionPrediction\data\';
cd(rootdir);
subfilz = dir('*');

for isub = 1:length(subfilz)
    sub = subfilz(isub).name;
    
    runfolders = dir([rootdir sub filesep '*resample']);
    folders2keep = [];
    for irun = 1:length(runfolders)
        if ~contains(runfolders(irun).name,'@raw')
            folders2keep = [folders2keep irun];
        end
    end
    runfolders = runfolders(folders2keep);
    sFiles = cell(1,525);
    for irun=1:length(runfolders)
        filz = dir([rootdir sub filesep runfolders(irun).name filesep '*trial*']);
        
        for ifile = 1:length(filz)
            sFiles{1,(irun-1)*length(filz)+ifile} = [sub filesep runfolders(irun).name filesep filz(ifile).name];
        end
    end
    
    sFiles = sFiles(~cellfun('isempty', sFiles));%if less than the initialized 525 trials, remove empty cells here
    
    % Start a new report
    bst_report('Start', sFiles);
    
    % Process: DC offset correction: [-500ms,-2ms]
    sFiles = bst_process('CallProcess', 'process_baseline', sFiles, [], ...
        'baseline',    [-0.5, -0.002], ...
        'sensortypes', 'MEG', ...
        'method',      'bl', ...  % DC offset correction:    x_std = x - &mu;
        'overwrite',   1);
    
    % Save and display report
    ReportFile = bst_report('Save', sFiles);
    bst_report('Open', ReportFile);
end
