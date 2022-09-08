close all; clearvars; clc

% load psychtoolbox output files
addpath(genpath('\\XXX\ActionPrediction\code\behaviouralAnalysis'));
indir = '\\XXX\ActionPrediction\data\MEG\PTBoutput';

% type of catch trial can be traced back from catch trials pool, which was created with "CreateCatchTrials.m" function
% so you can rerun the function here, or load the previously stored info:
load('\\XXX\ActionPrediction\code\behaviouralAnalysis\CatchTrialPool');

% subject files
subfilz = dir(fullfile(indir,'*block1*'));

data = struct;
RTpersub = zeros(length(subfilz),1);
ACCpersub = zeros(length(subfilz),1);
RTpercon = zeros(length(subfilz),2);
ACCpercon = zeros(length(subfilz),2);
for isub = 1:length(subfilz)
    
    % run files for current subject
    runfilz = dir(fullfile(indir,['*' subfilz(isub).name(end-12:end-11) '*']));
    
    % for subject 3 there was a problem in both run 3 and 6 and they were stopped halfway through
    if isub == 3
        runfilz([3 6]) = [];
    end
    
    % combine data from all runs
    for irun = 1:length(runfilz)
    
        load(fullfile(indir,runfilz(irun).name));
    
        if irun == 1
        	data = output;
        else
            data = [data output];
        end
        
    end
    
    % select only catch trials
    catchtrials = logical(extractfield(data,'catch'));
    data = data(catchtrials);
    
    if length(data) ~= length(runfilz)*20
        error('Not right catch trial amount...');
    end
    
    % extract reaction time and response
    RT = extractfield(data,'RT');
    response = extractfield(data,'correct_response');
    
    % if no response, RT is negative. correct this
    RT(RT < 0 | RT > 5) = NaN;
    
    % condition-average per subject
    RTpersub(isub) = nanmean(RT);
    ACCpersub(isub) = nanmean(response);
    
    % sort catch trials per type
    allCatchTrials = [CatchTrialsKeep CatchTrialsKinematics CatchTrialsSemantics];
    CatchKin = [1 8 9 10 12 13 15:28];
    CatchSem = [2:7 11 14 29:42];
    AllKinTrials = false(length(data),1);
    AllSemTrials = false(length(data),1);
    AllCrossTrials = false(length(data),1);
    for itrial = 1:length(data)
        currentCatchTrial = find(((cell2mat(allCatchTrials(1,:)) == data(itrial).sequence) + (cell2mat(allCatchTrials(2,:)) == data(itrial).time_catch)) == 2);
        if isempty(currentCatchTrial)
            AllCrossTrials(itrial) = true;
        elseif any(currentCatchTrial == CatchKin)
            AllKinTrials(itrial) = true;
        elseif any(currentCatchTrial == CatchSem)
            AllSemTrials(itrial) = true;
        end
            
    end

    if min(AllCrossTrials + AllSemTrials + AllKinTrials) ~= 1 || max(AllCrossTrials + AllSemTrials + AllKinTrials) ~= 1
        error('Error in defining catch trial types');
    end
    
    % combine catch trial types in single condition variable:
    conID = logical([AllCrossTrials AllKinTrials+AllSemTrials]);

    % Now determine performance per catch trial type
    for icon = 1:size(conID,2)
        RTpercon(isub,icon) = nanmean(RT(conID(:,icon)));
        ACCpercon(isub,icon) = nanmean(response(conID(:,icon)));
    end

end

%% Figures
% per condition
close(figure(1));
figure
set(gcf,'color','w');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [5 5 10 8]);

%accuracy
subplot(211)
hold on
plotSpread(ACCpercon,'spreadWidth',1,'distributionColors',[.7 .7 .7 ; .4 .4 .4],'distributionMarkers','.','showMM',2);
hold off
set(gca,'ylim',[0.5 1]);
set(gca,'xlim',[0.5 2.5]);
set(gca,'XTick',[1 2])
set(gca,'XTickLabel','')
set(gca,'YTick',[0.5 .75 1])
set(gca,'YTickLabel',{'50','75', '100'})
set(gca, 'FontSize', 8,'Fontname','Arial');
ylabel('Performance [%]','FontSize',8);

%RT
subplot(212)
hold on
plotSpread(RTpercon,'spreadWidth',1,'distributionColors',[.7 .7 .7 ; .4 .4 .4],'distributionMarkers','.','showMM',2);
hold off
set(gca,'ylim',[0.4 2]);
set(gca,'xlim',[0.5 2.5]);
set(gca,'XTick',[1 2])
set(gca,'XTickLabel',{'fixation cross','motion dancer'})
set(gca,'YTick',[.5 1 1.5 2])
set(gca,'YTickLabel',{'0.5', '1.0','1.5', '2.0'})
set(gca, 'FontSize', 8,'Fontname','Arial');
ylabel('RT [sec]','FontSize',8);
xlabel('catch trial type','FontSize',8);

% cd('XXX\ActionPrediction\figures');
% print -depsc -r600 ActionPrediction_FigureS1c_behaviour.eps

% % set(gcf,'renderer','Painters')
% print -depsc -tiff -r600 -painters ActionPrediction_FigureS1c_behaviour.eps
