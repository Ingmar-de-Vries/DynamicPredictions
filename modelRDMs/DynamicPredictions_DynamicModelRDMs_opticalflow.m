function DynamicPredictions_DynamicModelRDMs_opticalflow()

% set input and output directories    
addpath('//XXX/ActionPrediction/code/neuralDecoding');
rootdir = '//XXX/ActionPrediction';
StimDir = fullfile(rootdir,'experiment','stimuli','videos');
DataDir = fullfile(rootdir,'data','modelRDMs');

vidNames = dir(fullfile(StimDir, '*.mp4'));

%% Next load vector representations and compute RDMs 
frames = 250;
RDMoptflow_dir = zeros(length(vidNames),length(vidNames),frames,frames);
RDMoptflow_mag = zeros(length(vidNames),length(vidNames),frames,frames);
for imov1 = 1:length(vidNames)
    
    % create optical flow vector representations of video 1
    cfg = [];
    cfg.videoName = fullfile(StimDir,vidNames(imov1).name);
    [vecrep_dir1, vecrep_mag1] = DynamicPredictions_video2opticalflow(cfg);
    
    for imov2 = 1:length(vidNames)
        
        clc;
        disp(['Correlation optical flow representation video ' num2str(imov1) ' with video ' num2str(imov2)]);
        
        % create optical flow vector representations of video 2
        cfg = [];
        cfg.videoName = fullfile(StimDir,vidNames(imov2).name);
        [vecrep_dir2, vecrep_mag2] = DynamicPredictions_video2opticalflow(cfg);
        
        % create dynamic RDMs
        RDMoptflow_dir(imov1,imov2,:,:) = 1 - corr(vecrep_dir1',vecrep_dir2');
        RDMoptflow_mag(imov1,imov2,:,:) = 1 - corr(vecrep_mag1',vecrep_mag2');
        
    end
    
end

% interpolate to 100Hz to match neural and kinematic RDMs in main script
newframes = 500;
fsVid = 50;
fsNew = 100;
tVid = 0:1/fsVid:5-1/fsVid;
tNew = 0:1/fsNew:5-1/fsNew;
tempall1 = zeros(length(vidNames),length(vidNames),newframes,newframes);
tempall2 = tempall1;
for istim1 = 1:length(vidNames)
    for istim2 = 1:length(vidNames)
        
        temp = squeeze(RDMoptflow_dir(istim1,istim2,:,:));
        % bilinear interpolation:
        temp = interp1(tVid,temp,tNew,'pchip','extrap');
        temp = interp1(tVid,temp',tNew,'pchip','extrap')';
        tempall1(istim1,istim2,:,:) = temp;
        
        temp = squeeze(RDMoptflow_mag(istim1,istim2,:,:));
        % bilinear interpolation:
        temp = interp1(tVid,temp,tNew,'pchip','extrap');
        temp = interp1(tVid,temp',tNew,'pchip','extrap')';
        tempall2(istim1,istim2,:,:) = temp;
        
    end% second stimuli loop
end% first stimuli loop

RDMoptflow_dir = tempall1;
RDMoptflow_mag = tempall2;

save([DataDir filesep 'ActionPrediction_dynRDM_OFdir'],'RDMoptflow_dir');
save([DataDir filesep 'ActionPrediction_dynRDM_OFmag'],'RDMoptflow_mag');

