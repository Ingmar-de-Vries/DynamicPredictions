%% plot example video frames and stick figures, as used e.g., for Fig. 1a and b in article
% this script also contains info about which kinematic marker was placed on which location on the ballet dancer's body
close all; clearvars; clc;

%% set paths
vidir = '\\XXX\ActionPrediction\experiment\stimuli\videos';

vidName = dir(fullfile(vidir, '*.mp4'));

% pick which frames to show:
frames2show1 = 73:12:145;
frames2show2 = 93:12:165;

%% load variables
%video
videoHeader = VideoReader(fullfile(vidir,vidName(2).name));
video1 = read(videoHeader);
video1 = video1(21:380,51:326,:,frames2show1);

videoHeader = VideoReader(fullfile(vidir,vidName(7).name));
video2 = read(videoHeader);
video2 = video2(21:380,51:326,:,frames2show2);

% plot video frames
close(figure(1));
figure(1)
set(gcf,'color','w');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [2 2 10 5]);

for frame = 1:size(video1,4)
    
    subplot_tight(2,size(video1,4),frame)
    imshow(squeeze(video1(:,:,:,frame)));
    
    subplot_tight(2,size(video1,4),frame+size(video1,4))
    imshow(squeeze(video2(:,:,:,frame)));
    
end

%% example illustration of markers
kindir = '\\XXX\ActionPrediction\experiment\stimuli\kinematics50Hz';
kinName = dir(fullfile(kindir, '*kin*'));

load(fullfile(kindir,kinName(2).name));
kin(:,:,:,1) = squeeze(seqkin(:,:,frames2show1));

load(fullfile(kindir,kinName(7).name));
kin(:,:,:,2) = squeeze(seqkin(:,:,frames2show2));

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

close(figure(2));
figure(2)
set(gcf,'color','w');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [2 2 10 3]);

for ikin = 1:size(kin,3)
    subplot_tight(2,size(kin,3),ikin)
    hold on

    %head
    plot(kin(IDhead,1,ikin,1),kin(IDhead,3,ikin,1),'.k','MarkerSize',10);

    %bones
    for i = 1:size(IDbones,1)
        plot([kin(IDbones(i,1),1,ikin,1) kin(IDbones(i,2),1,ikin,1)],[kin(IDbones(i,1),3,ikin,1) kin(IDbones(i,2),3,ikin,1)],'k','LineWidth',1,'MarkerSize',10);
    end

    hold off
    set(gca,'xlim',[-900 900],'ylim',[0 2000]);
    set(gca,'xticklabel',[]);
    set(gca,'yticklabel',[]);
    box on
    
    subplot_tight(2,size(kin,3),ikin+size(kin,3))
    hold on

    %head
    plot(kin(IDhead,1,ikin,2),kin(IDhead,3,ikin,2),'.k','MarkerSize',10);

    %bones
    for i = 1:size(IDbones,1)
        plot([kin(IDbones(i,1),1,ikin,2) kin(IDbones(i,2),1,ikin,2)],[kin(IDbones(i,1),3,ikin,2) kin(IDbones(i,2),3,ikin,2)],'k','LineWidth',1,'MarkerSize',10);
    end

    hold off
    set(gca,'xlim',[-900 900],'ylim',[0 1800]);
    set(gca,'xticklabel',[]);
    set(gca,'yticklabel',[]);
    box on

end

% cd('XXX\ActionPrediction\figures');
% print -depsc -r600 ActionPrediction_figure1a_videoframes.eps
% print -depsc -r600 ActionPrediction_figure1b_stickfigures.eps

