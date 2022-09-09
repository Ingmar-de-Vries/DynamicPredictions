function [vecrep_dir, vecrep_mag] = DynamicPredictions_video2opticalflow(cfg)
%% Transform video frames to an optical flow vector representation for each frame
%INPUT
% video = full name (incl path) of video

% estimate optical flow field of each frame using the Farneback implementation
opticFlow = opticalFlowFarneback;

% load video
vidReader = VideoReader(cfg.videoName);

% initialize variables
matrep_dir = zeros(vidReader.NumFrames,2,vidReader.Height,vidReader.Width);
matrep_mag = zeros(vidReader.NumFrames,vidReader.Height,vidReader.Width);
clear flowFB
iframe = 0;
while hasFrame(vidReader)
    
    iframe = iframe+1;
    
    frameRGB = readFrame(vidReader);
    frameGray = rgb2gray(frameRGB);
    
    % estimate optical flow field
    flowFB = estimateFlow(opticFlow,frameGray);
    matrep_dir(iframe,1,:,:) = flowFB.Vx;
    matrep_dir(iframe,2,:,:) = flowFB.Vy;
    matrep_mag(iframe,:,:) = flowFB.Magnitude;
    
end

% optical flow first frame undefined (because determined based on previous frame), so for simplicity we say frame 1 = frame 2
matrep_dir(1,:,:,:) = matrep_dir(2,:,:,:);
matrep_mag(1,:,:) = matrep_mag(2,:,:);

%% last steps before creating vector, following Efros et al. 2003 
%smooth with Gaussian filter
sigma = 5;% ~ as Kriegeskorte 2008 for pixelwise similarity
for iframe = 1:size(matrep_dir,1)
    
    % x or y direction
    for idir = 1:2       
        matrep_dir(iframe,idir,:,:) = imgaussfilt(squeeze(matrep_dir(iframe,idir,:,:)),sigma);
    end
    
    % magnitude
    matrep_mag(iframe,:,:) = imgaussfilt(squeeze(matrep_mag(iframe,:,:)),sigma);
    
end

% store optical flow vector representation of each video
vecrep_dir = reshape(matrep_dir,size(matrep_dir,1),numel(squeeze(matrep_dir(1,:,:,:))));
vecrep_mag = reshape(matrep_mag,size(matrep_mag,1),numel(squeeze(matrep_mag(1,:,:,:))));

end