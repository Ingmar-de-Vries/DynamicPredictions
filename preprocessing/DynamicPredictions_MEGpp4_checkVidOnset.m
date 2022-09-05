% In some rare cases, we find more than 70 video onset events in a run, because of erroneously stored duplicate triggers.
% In this piece of code we search for those duplicates. 
% To double check you can compare it to the trigger values in the PTB output.
% After identifying them, remove them manually in the BST gui.
% Ingmar de Vries October 2020
clearvars;

% set paths
addpath('\\XXX\ActionPrediction\code\preprocessing');
filedir = '\\XXX\ActionPrediction\data\MEG\manualEventFixes\';
file = dir(fullfile(filedir, '*mrk*'));

file = file(24)

% Load BST events
BSTevents = readBSTevents([filedir filesep file.name]);

for itrial = 2:length(BSTevents)
    if BSTevents(itrial) <= BSTevents(itrial-1)+4
        disp(['Dulicate is trial ' num2str(itrial) ' at t = ' num2str(BSTevents(itrial))]);
    end
end