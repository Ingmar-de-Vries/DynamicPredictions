%% ===== CONVERT EPOCHED DATA FROM BRAINSTORM TO FIELDTRIP =====
% this script exports epoched data in fieldtrip ADAM decoding toolbox compatible format

clearvars;

% output directory
outdir = '\\XXX\ActionPrediction\data\MEG\PreprocessedSensor';
addpath('\\XXX\ActionPrediction\toolboxes\fieldtrip-20191113');
addpath('\\XXX\ActionPrediction\code\preprocessing');
ft_defaults

% set some parameters
cfg = [];
cfg.realign2photodiode = 1;
cfg.FTformat = 2;%1 = Fieldtrip data format after ft_preprocessing, 2 = Fieldtrip data format after ft_timelockanalysis

%% ===== CREATE FIELDTRIP STRUCTURE =====
% get the subjs from the actual database
s = bst_get('ProtocolSubjects');

load('\\XXX\ActionPrediction\data\MEG\PreprocessedSensor\allBadTrials');

%% ===== SUBJECTS LOOP =====

for iSubj = [1:numel(s.Subject)]
    
    % Extract Structure of Subj to process
    sFiles = bst_process('CallProcess', 'process_select_files_data', [], [], ...
        'subjectname',  s.Subject(iSubj).Name, ...
        'tag',          'VideoOnset');
    
    %select trials
    sFiles = sFiles(~contains(extractfield(sFiles,'Condition'),'High'));
    sFiles = sFiles(~contains(extractfield(sFiles,'Condition'),'more'));
    
    badTrials = allBadTrials.BadTrials{strcmp(allBadTrials.subject,s.Subject(iSubj).Name)};
    badidx = zeros(1,length(badTrials));
    
    if any(iSubj == [1:3 5:22])
        for itrial = 1:length(badTrials)
            badidx(itrial) = find(contains(extractfield(sFiles,'FileName'),badTrials{itrial}));
        end
    elseif iSubj == 4
        badidx = [];
    end
    badidx(badidx==0) = [];
    
    sFiles(badidx) = [];   
    
    FTdata	= struct;
    
    %% === LOOP THROUGH TRIALS ===
    
    for iTrial = 1:length(sFiles)
        
        % Export the brainstorm data structure of each trial into fieldtrip timelocked format
        [data, dataBST, ChannelMat] = out_fieldtrip_data(sFiles(iTrial).FileName, sFiles(iTrial).ChannelFile, ...
            'MEG,MISC008,STI102,EEG BAD LOC', 0);
        
        % Combining trials
        FTdata.time{1,iTrial}               = cell2mat(data.time);
        FTdata.trial{1,iTrial}				= cell2mat(data.trial);
        
        % Find trigger value
        for ievent = 1:length(dataBST.Events)
            if dataBST.Events(ievent).times == 0
                FTdata.trialinfo(iTrial,1) = str2double(dataBST.Events(ievent).label);
            end
        end

        % Some extra info not necessary for FT but handy nonetheless
        FTdata.BSTevents{iTrial}			= dataBST.Events;
        FTdata.original_trialnum{iTrial}             = dataBST.Comment;
        
    end
    
    % complete the FiledTrip Structure with common fields across trials
    FTdata.label							= data.label;
    FTdata.grad                             = data.grad;
    
    %% realign data to photodiode
    if cfg.realign2photodiode == 1
        FTdata = ActionPrediction_PP10_realign2photodiode(FTdata);
    end
    
    %% transform from FT preprocessed format to FT timelocked format, and remove some unnecessary fields
    if cfg.FTformat == 2
        
        FTdata = rmfield(FTdata,{'BSTevents','original_trialnum'});

        ft_cfg=[];
        ft_cfg.keeptrials = 'yes';
        ft_cfg.removemean = 'no';
        ft_cfg.channel = 'meg';
        FTdata = ft_timelockanalysis(ft_cfg,FTdata);
    end
    
    %% ===== WRITE THE SUBJECT STRUCTURE TO DISK =====
    
    save([outdir filesep 'preprocessedSensor_' s.Subject(iSubj).Name, '.mat'], 'FTdata', '-v7.3');
    
end
