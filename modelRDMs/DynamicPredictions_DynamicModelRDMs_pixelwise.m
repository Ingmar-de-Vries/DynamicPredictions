close all; clearvars; clc

% input and output directories
StimDir = '\\XXX\ActionPrediction\experiment\stimuli\videos';
DataDir = '\\XXX\ActionPrediction\data\modelRDMs';
addpath('\\XXX\ActionPrediction\Code\modelRDMs');

vidNames = dir(fullfile(StimDir, '*.mp4'));
    
% Loop over each combination of videos to fill in dynamic RDM matrix
frames = 250;
RDMgraysmooth = zeros(length(vidNames),length(vidNames),frames,frames);
for iMovie1=1:length(vidNames)
     
    %Call function to extract smoothed grayscale pixelwise luminance values as vector representations for first video
    cfg = [];
    cfg.videoName = fullfile(StimDir,vidNames(iMovie1).name);
    [~,~,~, vecrepGRAYSMOOTH1] = ProjectAction_video2vector(cfg);
    
    for iMovie2=1:length(vidNames)
        
        %Call function to extract smoothed grayscale pixelwise luminance values as vector representations for second video
        cfg = [];
        cfg.videoName = fullfile(StimDir,vidNames(iMovie2).name);
        [~,~,~, vecrepGRAYSMOOTH2] = DynamicPredictions_video2vector(cfg);

        clc;
        disp(['Correlating video ' num2str(iMovie1) ' with ' num2str(iMovie2)]);
        
        % compute dynamic RDM 
        RDMgraysmooth(iMovie1,iMovie2,:,:) = 1 - corr(vecrepGRAYSMOOTH1,vecrepGRAYSMOOTH2);
        
    end
end

% save correlation matrices
save([DataDir filesep 'DynamicPredictions_dynRDM_graysmooth'],'RDMgraysmooth');
