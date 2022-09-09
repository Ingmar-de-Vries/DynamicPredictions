% This script plots an illustration of all models for a single frame of 2 example videos
% It was used for creating Figure S2 in the article 

close all; clearvars; clc;
cmap = brewermap([],'*Greys');

%% set paths
vidir = '\\XXX\ActionPrediction\experiment\stimulil\videos';

vidName = dir(fullfile(vidir, '*.mp4'));

% example frame selected to illustrate what 'relative' means, i.e., the body posture should be very dissimilar for absolute, but very similar for
% relative
frame2show1 = 145;
frame2show2 = 165;

%% load videos
videoHeader = VideoReader(fullfile(vidir,vidName(2).name));
video1 = read(videoHeader);
video1 = video1(21:380,51:326,:,:);
frames(:,:,:,1) = video1(:,:,:,frame2show1);

videoHeader = VideoReader(fullfile(vidir,vidName(7).name));
video2 = read(videoHeader);
video2 = video2(21:380,51:326,:,:);
frames(:,:,:,2) = video2(:,:,:,frame2show2);

close(figure(1));
figure(1)
set(gcf,'color','w');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [2 2 16 6]);

%% show example frame
subplot_tight(2,8,1)
imshow(squeeze(frames(:,:,:,1)));

subplot_tight(2,8,9)
imshow(squeeze(frames(:,:,:,2)));

%% pixelwise smoothed grayscale
framesdouble = im2double(frames);% correctly transform uint8 to double

% create CIELAB values for each frame
framesLAB = rgb2lab(framesdouble);

% grayscale / luminance values approximated by the L (first) dimension in the CIELAB colorspace 
framesGRAY = squeeze(framesLAB(:,:,1,:));

% grayscale smoothed with Gaussian kernel
% Kriegeskorte et al. 2008 FSN uses Gaussian kernel with FWHM = 11.75 pixels
% Conversion to sigma is roughly: sigma = FWHM / 2.355, i.e.: 11.75/2.355 ~= 5
framesGRAYSMOOTH = imgaussfilt(framesGRAY,5);

climit = squeeze(max(max(framesGRAYSMOOTH)));

subplot_tight(2,8,2)
imshow(squeeze(framesGRAYSMOOTH(:,:,1)),[0 climit(1)]);

subplot_tight(2,8,10)
imshow(squeeze(framesGRAYSMOOTH(:,:,2)),[0 climit(2)]);

%% optical flow vector magnitude and direction
% estimate optical flow of each frame
opticFlowFB = opticalFlowFarneback;

% load video 1
vidReader = VideoReader(fullfile(vidir,vidName(2).name));

% initialize variables
clear flowFB
iframe = 0;
while hasFrame(vidReader)
    
    iframe = iframe+1;
    
    frameRGB = readFrame(vidReader);
    frameGray = rgb2gray(frameRGB);
    
    flowFB = estimateFlow(opticFlowFB,frameGray);
    
    if iframe == frame2show1
       vectors = flowFB;
       break
    end
    
end

subplot_tight(2,8,3)
imagesc(flowFB.Magnitude);
colormap(cmap);
set(gca,'xticklabel','');
set(gca,'yticklabel','');

subplot_tight(2,8,4)
imshow(frameRGB)
hold on
plot(flowFB,'DecimationFactor',[15 15],'ScaleFactor',10);
hold off

% load video 2
vidReader = VideoReader(fullfile(vidir,vidName(7).name));

% initialize variables
clear flowFB
iframe = 0;
while hasFrame(vidReader)
    
    iframe = iframe+1;
    
    frameRGB = readFrame(vidReader);
    frameGray = rgb2gray(frameRGB);
    
    flowFB = estimateFlow(opticFlowFB,frameGray);
    
    if iframe == frame2show2
       vectors = flowFB;
       break
    end
    
end

subplot_tight(2,8,11)
imagesc(flowFB.Magnitude);
colormap(cmap);
set(gca,'xticklabel','');
set(gca,'yticklabel','');

subplot_tight(2,8,12)
imshow(frameRGB)
hold on
plot(flowFB,'DecimationFactor',[15 15],'ScaleFactor',10);
hold off

%% marker positions
kindir = '\\XXX\ActionPrediction\experiment\stimuli\kinematics50Hz';
kinName = dir(fullfile(kindir, '*kin*'));

load(fullfile(kindir,kinName(2).name));
kin(:,:,1) = squeeze(seqkin(:,:,frame2show1));

load(fullfile(kindir,kinName(7).name));
kin(:,:,2) = squeeze(seqkin(:,:,frame2show2));

%define bones
%legs
IDbones(1,:) = [12 10];%left lower leg
IDbones(2,:) = [10 8];%left upper leg
IDbones(3,:) = [13 11];%right lower leg
IDbones(4,:) = [11 9];%right upper leg
%trunk
IDbones(5,:) = [8 2];%left trunk
IDbones(6,:) = [9 3];%right trunk
IDbones(7,:) = [8 9];%hips
IDbones(8,:) = [2 3];%shoulders
%arms
IDbones(9,:) = [6 4];%left lower arm
IDbones(10,:) = [4 2];%left upper arm
IDbones(11,:) = [7 5];%right lower arm
IDbones(12,:) = [5 3];%right upper arm
IDhead = 1;

subplot_tight(2,8,5)
hold on

%head
plot(kin(IDhead,1,1),kin(IDhead,3,1),'.k','MarkerSize',10);

%bones
for i = 1:size(IDbones,1)
    plot([kin(IDbones(i,1),1,1) kin(IDbones(i,2),1,1)],[kin(IDbones(i,1),3,1) kin(IDbones(i,2),3,1)],'k','LineWidth',1,'MarkerSize',10);
end

hold off
set(gca,'xlim',[-900 900],'ylim',[0 2000]);
set(gca,'xtick',1);
set(gca,'xticklabel','');
set(gca,'yticklabel','');
box on
axis square

subplot_tight(2,8,13)
hold on

%head
plot(kin(IDhead,1,2),kin(IDhead,3,2),'.k','MarkerSize',10);

%bones
for i = 1:size(IDbones,1)
    plot([kin(IDbones(i,1),1,2) kin(IDbones(i,2),1,2)],[kin(IDbones(i,1),3,2) kin(IDbones(i,2),3,2)],'k','LineWidth',1,'MarkerSize',10);
end

hold off
set(gca,'xlim',[-900 900],'ylim',[0 2000]);
set(gca,'xtick',1);
set(gca,'xticklabel','');
set(gca,'yticklabel','');
box on
axis square

%% relative positions
subplot_tight(2,8,6)
hold on

%head
plot(kin(IDhead,1,1),kin(IDhead,3,1),'.k','MarkerSize',10);

%bones
for i = 1:size(IDbones,1)
    plot([kin(IDbones(i,1),1,1) kin(IDbones(i,2),1,1)],[kin(IDbones(i,1),3,1) kin(IDbones(i,2),3,1)],'k','LineWidth',1,'MarkerSize',10);
end

hold off
set(gca,'xlim',[-900 900],'ylim',[0 2000]);
set(gca,'xtick',1);
set(gca,'xticklabel','');
set(gca,'yticklabel','');
box on
axis square

% now align frame 2 to frame 1 using procrustes
[~, kin(:,:,3)] = procrustes(squeeze(kin(:,:,1)),squeeze(kin(:,:,2)));

subplot_tight(2,8,14)
hold on

%head
plot(kin(IDhead,1,3),kin(IDhead,3,3),'.k','MarkerSize',10);

%bones
for i = 1:size(IDbones,1)
    plot([kin(IDbones(i,1),1,3) kin(IDbones(i,2),1,3)],[kin(IDbones(i,1),3,3) kin(IDbones(i,2),3,3)],'k','LineWidth',1,'MarkerSize',10);
end

hold off
set(gca,'xlim',[-900 900],'ylim',[0 2000]);
set(gca,'xtick',1);
set(gca,'xticklabel','');
set(gca,'yticklabel','');
box on
axis square

%% motion vectors
load(fullfile(kindir,kinName(2).name));
kinnext(:,:,1) = squeeze(seqkin(:,:,frame2show1+5));
load(fullfile(kindir,kinName(7).name));
kinnext(:,:,2) = squeeze(seqkin(:,:,frame2show2+5));
[~, kinnext(:,:,3)] = procrustes(squeeze(kinnext(:,:,1)),squeeze(kinnext(:,:,2)));

subplot_tight(2,8,7)
hold on

%head
plot(kin(IDhead,1,1),kin(IDhead,3,1),'.k','MarkerSize',10);
plot(kinnext(IDhead,1,1),kinnext(IDhead,3,1),'.','Color',[.4 .4 .4],'MarkerSize',10);

%bones
for i = 1:size(IDbones,1)
    plot([kin(IDbones(i,1),1,1) kin(IDbones(i,2),1,1)],[kin(IDbones(i,1),3,1) kin(IDbones(i,2),3,1)],'k','LineWidth',1,'MarkerSize',10);
end
%bones
for i = 1:size(IDbones,1)
    plot([kinnext(IDbones(i,1),1,1) kinnext(IDbones(i,2),1,1)],[kinnext(IDbones(i,1),3,1) kinnext(IDbones(i,2),3,1)],'Color',[.4 .4 .4],'LineWidth',1,'MarkerSize',10);
end

%motion vectors from each marker of one frame to the same marker on the next frame
for i = 1:size(kin,1)
    quiver(kin(i,1,1),kin(i,3,1),kinnext(i,1,1)-kin(i,1,1),kinnext(i,3,1)-kin(i,3,1),0,'b','lineWidth',2);
end
        
hold off
set(gca,'xlim',[-900 900],'ylim',[0 2000]);
set(gca,'xtick',1);
set(gca,'xticklabel','');
set(gca,'yticklabel','');
box on
axis square

subplot_tight(2,8,15)
hold on

%head
plot(kin(IDhead,1,2),kin(IDhead,3,2),'.k','MarkerSize',10);
plot(kinnext(IDhead,1,2),kinnext(IDhead,3,2),'.','Color',[.4 .4 .4],'MarkerSize',10);

%bones
for i = 1:size(IDbones,1)
    plot([kin(IDbones(i,1),1,2) kin(IDbones(i,2),1,2)],[kin(IDbones(i,1),3,2) kin(IDbones(i,2),3,2)],'k','LineWidth',1,'MarkerSize',10);
end
%bones
for i = 1:size(IDbones,1)
    plot([kinnext(IDbones(i,1),1,2) kinnext(IDbones(i,2),1,2)],[kinnext(IDbones(i,1),3,2) kinnext(IDbones(i,2),3,2)],'Color',[.4 .4 .4],'LineWidth',1,'MarkerSize',10);
end

%motion vectors from each marker of one frame to the same marker on the next frame
for i = 1:size(kin,1)
    quiver(kin(i,1,2),kin(i,3,2),kinnext(i,1,2)-kin(i,1,2),kinnext(i,3,2)-kin(i,3,2),0,'b','lineWidth',2);
end
        
hold off
set(gca,'xlim',[-900 900],'ylim',[0 2000]);
set(gca,'xtick',1);
set(gca,'xticklabel','');
set(gca,'yticklabel','');
box on
axis square

%% relative
subplot_tight(2,8,8)
hold on

%head
plot(kin(IDhead,1,1),kin(IDhead,3,1),'.k','MarkerSize',10);
plot(kinnext(IDhead,1,1),kinnext(IDhead,3,1),'.','Color',[.4 .4 .4],'MarkerSize',10);

%bones
for i = 1:size(IDbones,1)
    plot([kin(IDbones(i,1),1,1) kin(IDbones(i,2),1,1)],[kin(IDbones(i,1),3,1) kin(IDbones(i,2),3,1)],'k','LineWidth',1,'MarkerSize',10);
end
%bones
for i = 1:size(IDbones,1)
    plot([kinnext(IDbones(i,1),1,1) kinnext(IDbones(i,2),1,1)],[kinnext(IDbones(i,1),3,1) kinnext(IDbones(i,2),3,1)],'Color',[.4 .4 .4],'LineWidth',1,'MarkerSize',10);
end

%motion vectors from each marker of one frame to the same marker on the next frame
for i = 1:size(kin,1)
    quiver(kin(i,1,1),kin(i,3,1),kinnext(i,1,1)-kin(i,1,1),kinnext(i,3,1)-kin(i,3,1),0,'b','lineWidth',2);
end
        
hold off
set(gca,'xlim',[-900 900],'ylim',[0 2000]);
set(gca,'xtick',1);
set(gca,'xticklabel','');
set(gca,'yticklabel','');
box on
axis square

subplot_tight(2,8,16)
hold on

%head
plot(kin(IDhead,1,3),kin(IDhead,3,3),'.k','MarkerSize',10);
plot(kinnext(IDhead,1,3),kinnext(IDhead,3,3),'.','Color',[.4 .4 .4],'MarkerSize',10);

%bones
for i = 1:size(IDbones,1)
    plot([kin(IDbones(i,1),1,3) kin(IDbones(i,2),1,3)],[kin(IDbones(i,1),3,3) kin(IDbones(i,2),3,3)],'k','LineWidth',1,'MarkerSize',10);
end
%bones
for i = 1:size(IDbones,1)
    plot([kinnext(IDbones(i,1),1,3) kinnext(IDbones(i,2),1,3)],[kinnext(IDbones(i,1),3,3) kinnext(IDbones(i,2),3,3)],'Color',[.4 .4 .4],'LineWidth',1,'MarkerSize',10);
end

%motion vectors from each marker of one frame to the same marker on the next frame
for i = 1:size(kin,1)
    quiver(kin(i,1,3),kin(i,3,3),kinnext(i,1,3)-kin(i,1,3),kinnext(i,3,3)-kin(i,3,3),0,'b','lineWidth',2);
end
        
hold off
set(gca,'xlim',[-900 900],'ylim',[0 2000]);
set(gca,'xtick',1);
set(gca,'xticklabel','');
set(gca,'yticklabel','');
box on
axis square

% cd('XXX\ActionPrediction\figures');
% print -depsc -r600 ActionPrediction_figureS2_stimulusmodels.eps

