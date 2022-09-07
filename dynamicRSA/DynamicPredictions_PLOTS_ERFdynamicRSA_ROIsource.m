function DynamicPredictions_PLOTS_ERFdynamicRSA_ROIsource(cfg)

addpath('\\cimec-storage5.unitn.it\MORWUR\Projects\INGMAR\ActionPrediction\toolboxes\fieldtrip-20191113');
addpath('\\cimec-storage5.unitn.it\MORWUR\Projects\INGMAR\ActionPrediction\code\neuralDecoding');

% addpath('\\XXX\ActionPrediction\toolboxes\fieldtrip-20191113');
% addpath(genpath('\\XXX\ActionPrediction\code'));
% addpath(genpath('\\XXX\matlab_toolboxes\CoSMoMVPA-master'));
% cfg.path = '\\XXX\ActionPrediction';

ft_defaults

%% some parameters for plotting:
models2plot = [1:7 10];
maxLag = 1;% maximum lag to plot
tRange=maxLag*cfg.downsample;% maximum lag in samples
TimeVec = 0:1/cfg.downsample:5-cfg.randshuff(2)-1/cfg.downsample;% dRSA time vector
% binVec = dsearchn(TimeVec',[0:ceil(max(TimeVec))]');

TimeVec2 = -maxLag:1/cfg.downsample:maxLag;% time vector for dRSA lag plot (i.e., [-maxLag maxLag])
binVec2 = [1  tRange+1  tRange*2+1];% for X-ticks in figure

% colors
cmap = brewermap([],'*RdBu');% colormap for time-time plots
lineColors = brewermap(7,'YlGnBu');% colormap for lineplots 
lineColors(1,:) = [];% first light yellow is too bright...
lineColors(1,:) = lineColors(1,:) - .1;% second still too bright so we darken it a little...

% fontsize
fs = 7;

%% input and output folders
if cfg.glmRSA == 0
    corrORglm = 'corr';
elseif cfg.glmRSA == 1
    corrORglm = ['pcaANDpcr_' num2str(cfg.nPCAcomps) 'comps'];
end
indirSTAT = fullfile(cfg.path,'data','MEG',['source_' cfg.atlas],'RSA','statistics',[corrORglm '_' num2str(cfg.randshuff(1)) 'iterations_' num2str(ceil(cfg.randshuff(2)*1000)) 'msec']);

% ROI names:
load(fullfile(cfg.path,'code','neuralDecoding','ROIdefinitions'),'ROIdefinition');
ROIdefinition = ROIdefinition.names;
ROInames = vertcat(ROIdefinition{:});
ROInames = ROInames(cfg.ROIVec);

rdmLabels = {'pixelwise grayscale','optical flow magnitude','optical flow direction','absolute position','relative position','absolute motion','relative motion','eye position'};

% initialize variables
peaksAll = zeros(length(cfg.SubVec),length(ROInames),length(rdmLabels),2);
omegaROIs = [];
% rstackROIs = [];
rstackLineSEM = [];
rstackLineROIs = [];
signLineROIs = [];
signLineStrictROIs = [];
signRS_ROIs = [];
RS_SEM = [];
RS_ROIs = [];
for iROI = 1:length(ROInames)
    
    %% load data resulting from STATS script
    ROIname = ROInames{iROI};
    fn = sprintf('%s%cSTATS_p%s_p%s_fisherz%d_%dHz_%s_smMEG%d_smRDMneu%d_smRDMmod%d',indirSTAT, filesep, '01', '05', cfg.fisherz, cfg.downsample, ROIname, cfg.smoothingMSec, cfg.smoothNeuralRDM, cfg.smoothModelRDM);
    load(fn,'omegaAll','rstackLine','signLine','RS','peakLatency');
    
    % stricter significance threshold for additional significance bars in line plots and for posthoc representational spread (RS) analysis
    fn = sprintf('%s%cSTATS_p%s_p%s_fisherz%d_%dHz_%s_smMEG%d_smRDMneu%d_smRDMmod%d',indirSTAT, filesep, '001', '01', cfg.fisherz, cfg.downsample, ROIname, cfg.smoothingMSec, cfg.smoothNeuralRDM, cfg.smoothModelRDM);
    strictSTATS = load(fn,'signLine','signRS');
    
    % combine data from all ROIs
    omegaROIs(iROI,:,:,:) = squeeze(mean(omegaAll));
    rstackLineSEM(iROI,:,:) = squeeze(stdErrMean(rstackLine,1));
    rstackLineROIs(iROI,:,:) = squeeze(mean(rstackLine));
    RS_SEM(iROI,:,:,:) = squeeze(stdErrMean(RS,1));
    RS_ROIs(iROI,:,:,:) = squeeze(mean(RS));
    
    signLineROIs(iROI,:,:) = signLine;
    signLineStrictROIs(iROI,:,:) = strictSTATS.signLine;
    signRS_ROIs(iROI,:,:,:) = strictSTATS.signSRI;
    peaksAll(:,iROI,:,:) = peakLatency;
    
end% iROI loop

% change zeros to NaN so it's easier to plot significance lines while leaving gaps for non-significant values:
signLineROIs(signLineROIs==0) = NaN;
signRS_ROIs(signRS_ROIs==0) = NaN;
signLineStrictROIs(signLineStrictROIs==0) = NaN;

%% dRSA lag-plots
close(figure(1));

figure(1);
set(gcf,'color','w');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [1 5 22 10]);

ylims4lines = [.25 .13 .04 .04 .04 .04 .04 .05];
for iRDM = 1:length(models2plot)
    
    subplot(2,4,iRDM);
    
    ylimit = ylims4lines(iRDM);
    rstack2plot = squeeze(rstackLineROIs(:,iRDM,:));
    SEM2plot = permute(repmat(squeeze(rstackLineSEM(:,iRDM,:)),[1 1 2]),[2 3 1]);% this is the format that "boundedline.m" wants
    
    hold on;
    clear h
    
    % plot lines with standard error shading around it
    boundedline(TimeVec2,rstack2plot, SEM2plot , 'alpha','cmap', lineColors(1:6,:));
    
    for isubROI = 1:size(rstack2plot,1)/(cfg.splitLR+1)% loop over ROIs
        
        % plot lines again but now wider and on top of the shadings
        h(isubROI) = plot(TimeVec2,rstack2plot(isubROI,:),'color',lineColors(isubROI,:),'lineWidth',1.5);
        
        % plot horizontal signifance bars 
        plot(TimeVec2,(-ylimit*0.2-(ylimit/18*(isubROI-1)))*squeeze(signLineROIs(isubROI,iRDM,:)),'lineWidth',3,'Color',lineColors(isubROI,:));
        if isubROI<5
            plot(TimeVec2,(-ylimit*0.2-(ylimit/18*(isubROI-1)))*squeeze(signLineStrictROIs(isubROI,iRDM,:)),'lineWidth',1.5,'Color','k');
        elseif isubROI>4
            plot(TimeVec2,(-ylimit*0.2-(ylimit/18*(isubROI-1)))*squeeze(signLineStrictROIs(isubROI,iRDM,:)),'lineWidth',1.5,'Color','w');
        end
        
    end
    
    % some settings for plotting
    plot([TimeVec2(1) TimeVec2(end)],[0 0],'--k');
    plot([0 0],[-1 1],'--k');
    set(gca,'xlim',[TimeVec2(1) TimeVec2(end)]);
    set(gca,'ylim',[-ylimit*.5 ylimit]);
    set(gca,'xtick',TimeVec2(binVec2));
    set(gca,'Fontsize',6,'FontName','Helvetica');
    hold off
    
    if iRDM == 1
        xlabel('lag [sec]','Fontsize',fs,'FontName','Helvetica');
    end
    if iRDM == 1
        ylabel('dynRSA [Rho]','Fontsize',fs,'FontName','Helvetica');
    end
    title(rdmLabels{iRDM},'Fontsize',fs,'FontName','Helvetica','FontWeight','normal');

end% iRDM

% cd('XXX\ActionPrediction\figures');
% print -depsc -r600 ActionPrediction_Figure2a_ROIlineplots.eps

% or try like this, if figures not stored correctly:  
% set(gcf,'renderer','Painters')
% print -depsc -tiff -r600 -painters ActionPrediction_Figure2a_ROIlineplots.eps

%% time-time 2D dRSA plot
% downsample because otherwise gigantic file size, and such a high resolution is not necessary anyway
omegaROIs = omegaROIs(:,1:5:end,1:5:end,:);
TimeVec = TimeVec(1:5:end);
binVec = [1 20 40 60];

close(figure(2));

figure(2);
set(gcf,'color','w');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [0 0 19 14]);
climits = squeeze(max(max(max(abs(omegaROIs)),[],2),[],3));
count = 0;
for iROI = 1:length(cfg.ROIVec)
    for iRDM = 1:length(models2plot)
        count = count + 1;
        subplot(6,length(models2plot),count)
        hold on
        contourf(TimeVec,TimeVec,squeeze(omegaROIs(iROI,:,:,iRDM)),60,'lineColor','none');
        plot(TimeVec,TimeVec,'--k','lineWidth',1)
        
        climit = climits(iRDM);
        caxis([-climit climit]);
        
        set(gca,'YDir','normal');
        set(gca,'xtick',TimeVec(binVec));
        set(gca,'xticklabel','');
        set(gca,'ytick',TimeVec(binVec));
        set(gca,'yticklabel','');
        
        if iRDM == 1
            ylabel(ROInames{iROI},'FontSize',fs);
        end
        
        set(gca,'tickLength',[0 0]);
        
        colormap(cmap);
        axis square
        box on
        hold off
                
    end
end

% cd('XXX\ActionPrediction\figures');
% set(gcf,'renderer','Painters')
% print -depsc -tiff -r600 -painters ActionPrediction_FigureS5_ROItimeXtime.eps

%% representational spread (RS) plot
% only plot ROI x model combi if significant in main dRSA analysis
close(figure(3));

figure(3);
set(gcf,'color','w');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [1 5 22 8]);

ylims4lines = [1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.2];
for iRDM = 1:length(models2plot)
    
    if iRDM == 3% for optical flow vector direction there are 2 peaks
        
        intervals = [1 116 ; 91 length(TimeVec2)];% only plot until through in Fig. 2a, because after that the other peak starts
        for ipeak = 1:2
            subplot(2,8,4+ipeak);
            
            ylimit = ylims4lines(iRDM);
            SRI2plot = squeeze(RS_ROIs(:,iRDM,ipeak,:));
            SEM2plot = permute(repmat(squeeze(RS_SEM(:,iRDM,2,:)),[1 1 2]),[2 3 1]);
            
            hold on;
            clear h
            
            for isubROI = 1:size(SRI2plot,1)
                if max(signLineROIs(isubROI,iRDM,:)) == 1
                    id2plot = ~isnan(SRI2plot(isubROI,:));% make sure not to plot NaNs coz boundedline doesn't handle them well
                    boundedline(TimeVec2(id2plot),SRI2plot(isubROI,id2plot), SEM2plot(id2plot,:,isubROI), 'alpha','cmap', lineColors(isubROI,:));
                end
            end
            for isubROI = 1:size(SRI2plot,1)
                if max(signLineROIs(isubROI,iRDM,:)) == 1
                    h(isubROI) = plot(TimeVec2,SRI2plot(isubROI,:),'color',lineColors(isubROI,:),'lineWidth',1.5);
                    plot(TimeVec2,(-ylimit*0.2-(ylimit/18*(isubROI-1)))*squeeze(signRS_ROIs(isubROI,iRDM,:,2)),'lineWidth',3,'Color',lineColors(isubROI,:));
                end
            end
            
            plot([TimeVec2(1) TimeVec2(end)],[0 0],'--k');
            plot([0 0],[-1 1],'--k');
            set(gca,'xlim',[TimeVec2(intervals(ipeak,1)) TimeVec2(intervals(ipeak,2))]);
            set(gca,'ylim',[-ylimit*.5 ylimit]);
            set(gca,'xtick',TimeVec2(binVec2));
            set(gca,'Fontsize',6,'FontName','Helvetica');
            hold off
        end
    else% for the other models there is only a single peak
        % averaged over video time
        subplot(2,4,iRDM);
        
        ylimit = ylims4lines(iRDM);
        SRI2plot = squeeze(RS_ROIs(:,iRDM,2,:));
        SEM2plot = permute(repmat(squeeze(RS_SEM(:,iRDM,2,:)),[1 1 2]),[2 3 1]);
        
        hold on;
        clear h
        
        for isubROI = 1:size(SRI2plot,1)
            if max(signLineROIs(isubROI,iRDM,:)) == 1
                id2plot = ~isnan(SRI2plot(isubROI,:));% make sure not to plot NaNs coz boundedline doesn't handle them well                    
                boundedline(TimeVec2(id2plot),SRI2plot(isubROI,id2plot), SEM2plot(id2plot,:,isubROI), 'alpha','cmap', lineColors(isubROI,:));
            end
        end
        for isubROI = 1:size(SRI2plot,1)
            if max(signLineROIs(isubROI,iRDM,:)) == 1
                h(isubROI) = plot(TimeVec2,SRI2plot(isubROI,:),'color',lineColors(isubROI,:),'lineWidth',1.5);
                plot(TimeVec2,(-ylimit*0.2-(ylimit/18*(isubROI-1)))*squeeze(signRS_ROIs(isubROI,iRDM,:,2)),'lineWidth',3,'Color',lineColors(isubROI,:));  
            end
        end
        
        plot([TimeVec2(1) TimeVec2(end)],[0 0],'--k');
        plot([0 0],[-1 1],'--k');
        set(gca,'xlim',[TimeVec2(1) TimeVec2(end)]);
        set(gca,'ylim',[-ylimit*.5 ylimit]);
        set(gca,'xtick',TimeVec2(binVec2));
        set(gca,'Fontsize',6,'FontName','Helvetica');
        hold off
        
    end
    
    if iRDM == 1
        xlabel('lag [sec]','Fontsize',fs,'FontName','Helvetica');
    end
    if iRDM == 1
        ylabel('SRI','Fontsize',fs,'FontName','Helvetica');
    end
    title(rdmLabels{iRDM},'Fontsize',fs,'FontName','Helvetica','FontWeight','normal');

end% iRDM

% cd('XXX\ActionPrediction\figures');
% print -depsc -r600 ActionPrediction_Figure3b_SR.eps

% set(gcf,'renderer','Painters')
% print -depsc -tiff -r600 -painters ActionPrediction_Figure3b_SR.eps

%% barplots of peak latency
close(figure(4));

figure(4)
set(gcf,'color','w');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [1 5 22 4]);

limits = [0 0.3 ; 0 0.3 ; -.7 0 ; 0 0.3 ; 0 0.3 ; 0 0.3 ; -.7 0 ; -.7 0];

subplotcount = 0;
for iRDM = 1:length(models2plot)-1
    
    if iRDM == 3
        for ipeak = 1:2
            subplotcount = subplotcount+1;
            subplot(1,8,subplotcount)
            bar2plot = squeeze(mean(peaksAll(:,end:-1:1,iRDM,ipeak)));
            sem2plot = squeeze(std(peaksAll(:,end:-1:1,iRDM,ipeak)))/sqrt(size(peaksAll,1));
            
            hold on
            for ibar = 1:size(bar2plot,2)
                barh(ibar,bar2plot(ibar),'FaceColor',lineColors(end+1-ibar,:));
                errorbar(bar2plot(ibar),ibar,sem2plot(ibar),'k.','horizontal');
            end
            hold off
            set(gca,'YTick',1:6);
            set(gca,'YTickLabel','');
            set(gca,'xlim',limits(subplotcount,:));
        end
    else
        subplotcount = subplotcount+1;
        subplot(1,8,subplotcount);
        bar2plot = squeeze(mean(peaksAll(:,end:-1:1,iRDM,2)));
        sem2plot = squeeze(std(peaksAll(:,end:-1:1,iRDM,2)))/sqrt(size(peaksAll,1));
        
        hold on
        for ibar = 1:size(bar2plot,2)
            barh(ibar,bar2plot(ibar),'FaceColor',lineColors(end+1-ibar,:));
            errorbar(bar2plot(ibar),ibar,sem2plot(ibar),'k.','horizontal');
        end
        hold off
        set(gca,'YTick',1:6);
        set(gca,'YTickLabel','');
        set(gca,'xlim',limits(subplotcount,:));
    end

end

% cd('XXX\ActionPrediction\figures');
% print -depsc -r600 ActionPrediction_Figure3a_peakLatency.eps

% % set(gcf,'renderer','Painters')
% % print -depsc -tiff -r600 -painters ActionPrediction_figure3a_peakLatency.eps
%
%% stats on peak latency (ANOVAs)
% peak timing for pixelwise vs optical flow magnitude, i.e. 2 models x 4 ROIs
Y = [];
MODEL = [];
ROI = [];
SUB = [];
for iSub = 1:length(cfg.SubVec)
    for iRDM = 1:2% all models (i.e., motion) with a predictive peak
        for iROI = 1:length(cfg.ROIVec)-2
            Y = [Y peaksAll(iSub, iROI ,iRDM, 2)];
            MODEL = [MODEL iRDM];
            ROI = [ROI iROI];
            SUB = [SUB iSub];
        end
    end
end

FACTORS={MODEL ROI SUB};
[P,table,STATS,TERMS] = anovan(Y, FACTORS, 'random',3,'model','full','display','on','varnames', {'MODEL' 'ROI' 'SUBJECT'});

% peak timing for predictive peaks, i.e. 3 models x 4 ROIs
peaks2compare = [3 1; 6 2; 7 2];% 3 x 2 matrix with model number in first column, and (predictive) peak in second column
Y = [];
MODEL = [];
ROI = [];
SUB = [];
for iSub = 1:length(cfg.SubVec)
    RDMcount = 0;
    for iRDM = peaks2compare(:,1)'% all models (i.e., motion) with a predictive peak
        RDMcount = RDMcount+1;
        for iROI = 1:length(cfg.ROIVec)-2
            Y = [Y peaksAll(iSub, iROI ,iRDM, peaks2compare(RDMcount,2))];
            MODEL = [MODEL RDMcount];
            ROI = [ROI iROI];
            SUB = [SUB iSub];
        end
    end
end

FACTORS={MODEL ROI SUB};
[P,table,STATS,TERMS] = anovan(Y, FACTORS, 'random',3,'model','full','display','on','varnames', {'MODEL' 'ROI' 'SUBJECT'});


end