function [vecrepGRAYSMOOTH] = DynamicPredictions_video2vector(cfg)
%% Transform video frames to a low level vector representation for each frame
%INPUT
% video = full name (incl path) of video

videoHeader = VideoReader(cfg.videoName);
video = read(videoHeader);

% transform uint8 to double
video = im2double(video);

% Vector representation of CIELAB values for each frame
videoLAB = rgb2lab(video);

% Vector representation of grayscale / luminance values approximated by the L (first) dimension in the CIELAB colorspace 
videoGRAY = squeeze(videoLAB(:,:,1,:));

% Vector representation of grayscale smoothed with Gaussian kernel
% Kriegeskorte et al. 2008 FSN uses Gaussian kernel with FWHM = 11.75 pixels
% Conversion to sigma is roughly: sigma = FWHM / 2.355, i.e.: 11.75/2.355 = 5ish
videoGRAYSMOOTH = imgaussfilt(videoGRAY,5);
vecrepGRAYSMOOTH = reshape(videoGRAYSMOOTH,numel(videoGRAYSMOOTH(:,:,1)),size(videoGRAYSMOOTH,3));

end