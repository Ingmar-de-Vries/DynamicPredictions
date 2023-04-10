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
    vecrepGRAYSMOOTH1 = DynamicPredictions_video2vector(cfg);
    
    for iMovie2=1:length(vidNames)
        
        %Call function to extract smoothed grayscale pixelwise luminance values as vector representations for second video
        cfg = [];
        cfg.videoName = fullfile(StimDir,vidNames(iMovie2).name);
        vecrepGRAYSMOOTH2 = DynamicPredictions_video2vector(cfg);

        clc;
        disp(['Correlating video ' num2str(iMovie1) ' with ' num2str(iMovie2)]);
        
        % compute dynamic RDM 
        RDMgraysmooth(iMovie1,iMovie2,:,:) = 1 - corr(vecrepGRAYSMOOTH1,vecrepGRAYSMOOTH2);
        
    end
end

% interpolate to 100Hz to match neural and kinematic RDMs in main script
newframes = 500;
fsVid = 50;
fsNew = 100;
tVid = 0:1/fsVid:5-1/fsVid;
tNew = 0:1/fsNew:5-1/fsNew;
tempall = zeros(length(vidNames),length(vidNames),newframes,newframes);
for istim1 = 1:length(vidNames)
    for istim2 = 1:length(vidNames)
        
        temp = squeeze(RDMgraysmooth(istim1,istim2,:,:));
        
        % bilinear interpolation:
        temp = interp1(tVid,temp,tNew,'pchip','extrap');
        temp = interp1(tVid,temp',tNew,'pchip','extrap')';
        tempall(istim1,istim2,:,:) = temp;
        
    end% second stimuli loop
end% first stimuli loop

RDMgraysmooth = tempall;

% save correlation matrices
save([DataDir filesep 'DynamicPredictions_dynRDM_graysmooth'],'RDMgraysmooth');
