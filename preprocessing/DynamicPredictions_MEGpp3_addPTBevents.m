clearvars; clc
% each run contains 70 'normal' trials, 15 catch trials that will be rejected for analyses,
% and 5 catch trials with the catch onset at 5 sec, which can thus be
% included in the analyses. This piece of code changes the trigger values
% for those 5 trials. Only for first 6 subjects, after that it was
% changed in the experiment itself (i.e. the PTB script) 

% set paths
addpath('\\XXX\ActionPrediction\code\preprocessing');
filedir = '\\XXX\ActionPrediction\data\MEG';

subfilzBST = dir([filedir filesep 'manualEventFixes' filesep '*r1.mrk*']);
subfilzPTB = dir([filedir filesep 'PTBoutput' filesep '*block1*']);

for isub = 1:length(subfilzBST)% Loop over subjects

    % check which runs are present for this subject
    subject = subfilzBST(isub).name(end-9:end-6);
    runfilzBST = dir([filedir filesep 'manualEventFixes' filesep '*trials_' subject '*.mrk']);
    
    subject = subfilzPTB(isub).name(1:end-5);
    runfilzPTB = dir([filedir filesep 'PTBoutput' filesep '*' subject '*']);    
    
    if length(runfilzBST) ~= length(runfilzPTB)
        warning('run files do not match between BST and PTB data');% Just a little safety check to make sure there's 6 run files in both the BST and PTB data 
        break
    end
    
    for irun = 1:length(runfilzBST)% Loop over runs
    
        % Load BST events 
        BSTevents = readBSTevents([filedir filesep 'manualEventFixes' filesep runfilzBST(irun).name]);
        BSTevents(isnan(BSTevents)) = [];
        BSTevents = sort(BSTevents);
        
        % Load PTB data       
        load([filedir filesep 'PTBoutput' filesep runfilzPTB(irun).name]);
        
        %Another safety check:
        clc
        disp(['You have just loaded in the following two files: ' runfilzBST(irun).name ' and ' runfilzPTB(irun).name ', check if that is correct']);
        
        PTBevents = extractfield(output,'trigger');
        catchOnset = extractfield(output,'time_catch'); 
        
        PTBevents2keep = PTBevents(catchOnset == 5);
        
        %% Create new event file
        % one for each direction
        vidOnset = BSTevents(catchOnset == 5);
        
        if length(PTBevents)~=length(BSTevents)
            error('Oops! Trigger lists do not match, something is wrong...');
        end
        
        %And one more safety check:
        disp(['Files will be saved as:' 10 subfilzBST(isub).name(1:end-6) '_r' num2str(irun) '_vidOnset' 10 ]);
        
        save([filedir filesep 'manualEventFixes' filesep subfilzBST(isub).name(1:end-6) '_r' num2str(irun) '_vidOnset'],'vidOnset');

    end
end
%% Some info
% just as a reminder here subject list:
% subject number in PTB     subject name in BST
% 01                            XXXX
% 02                            YYYY 
% 03                            ZZZZ
% 04                            ....


