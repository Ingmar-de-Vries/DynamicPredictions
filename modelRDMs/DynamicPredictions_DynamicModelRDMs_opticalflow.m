function DynamicPredictions_DynamicModelRDMs_opticalflow()

% set input and output directories    
addpath('//XXX/ActionPrediction/code/neuralDecoding');
rootdir = '//XXX/ActionPrediction';
StimDir = fullfile(rootdir,'experiment','stimuli','videos');
DataDir = fullfile(rootdir,'data','modelRDMs');

vidNames = dir(fullfile(StimDir, '*.mp4'));

%% Next load vector representations and compute RDMs 
frames = 250;
RDMopticalflow_dir = zeros(length(vidNames),length(vidNames),frames,frames);
RDMopticalflow_mag = zeros(length(vidNames),length(vidNames),frames,frames);
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
        RDMopticalflow_dir(imov1,imov2,:,:) = 1 - corr(vecrep_dir1',vecrep_dir2');
        RDMopticalflow_mag(imov1,imov2,:,:) = 1 - corr(vecrep_mag1',vecrep_mag2');
        
    end
    
end

save([DataDir filesep 'DynamicPredictions_dynRDM_opticalFlow'],'RDMopticalflow_dir','RDMopticalflow_mag');


