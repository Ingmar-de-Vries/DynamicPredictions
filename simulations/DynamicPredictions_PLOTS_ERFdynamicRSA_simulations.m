function DynamicPredictions_PLOTS_ERFdynamicRSA_simulations(cfg)

addpath('\\XXX\ActionPrediction\toolboxes\fieldtrip-20191113');
addpath(genpath('\\XXX\ActionPrediction\code'));
cfg.path = '\\XXX\ActionPrediction';
ft_defaults

addpath('\\cimec-storage5.unitn.it\MORWUR\Projects\INGMAR\ActionPrediction\toolboxes\fieldtrip-20191113');
addpath('\\cimec-storage5.unitn.it\MORWUR\Projects\INGMAR\ActionPrediction\code\neuralDecoding');

ft_defaults
cfg = ActionPrediction_config(cfg);

%% some parameters for plotting:
tRange=1*cfg.downsample;% maximum lag to plot
maxLag = 1;% maximum lag to plot
TimeVec2 = -maxLag:1/cfg.downsample:maxLag;% time vector for dRSA lag plot (i.e., [-maxLag maxLag])
binVec2 = [1  tRange+1  tRange*2+1];% for X-ticks in figure

% colors
colors = brewermap(9,'*RdYlBu');
colors(1:5,:) = colors(1:5,:) - .1;% first 5 colors are a bit too bright
colors(10,:) = [0 0 0];% last color is eye position which can be black

%% input and output folders
if cfg.glmRSA == 0
    corrORglm = 'corr';
elseif cfg.glmRSA == 1
    corrORglm = ['pcr_' num2str(cfg.nPCAcomps) 'comps'];
end
indirSIM = fullfile(cfg.path,'data','MEG',['source_' cfg.atlas],'RSA','simulations',[corrORglm '_' num2str(cfg.lag*1000) 'lag_' num2str(cfg.randshuff(1)) 'iterations_' num2str(ceil(cfg.randshuff(2)*1000)) 'msec']);

rdmLabels = cfg.dynRDMnames;

%% load simulation results
omega = [];
for ibatch = 1:cfg.iterbatches
    for implant = 1:length(rdmLabels)
        fn = sprintf('%s%comega_SUB%02d_implant%d_batch%d_%dHz_smMEG%d_smRDMneu%d_smRDMmod%d',indirSIM, filesep, 7, implant, ibatch, cfg.downsample, cfg.smoothingMSec, cfg.smoothNeuralRDM, cfg.smoothModelRDM);
        load(fn,'omegaAll');
        omega(ibatch,implant,:,:,:) = omegaAll;
    end
    
    % randomized data with no models implanted
    fn = sprintf('%s%comega_SUB%02d_implant%d_batch%d_%dHz_smMEG%d_smRDMneu%d_smRDMmod%d',indirSIM, filesep, 7, 99, ibatch, cfg.downsample, cfg.smoothingMSec, cfg.smoothNeuralRDM, cfg.smoothModelRDM);
    load(fn,'omegaAll');
    omega(ibatch,length(rdmLabels)+1,:,:,:) = omegaAll;
end
omega = squeeze(mean(omega));

%% compute rstackLine for lag-plots (see main script "DynamicPredictions_STATS_ERFdynamicRSA_ROIsource.m" for more info on this step)
rstack = zeros(size(omega,1),size(omega,4),size(omega,2),length(TimeVec2));
for implant = 1:size(omega,1)
    for iRDM = 1:size(omega,4)
        for iModelTime = 1:size(omega,2)
            
            timeidx = iModelTime - tRange:iModelTime + tRange;
            NotInVid = logical((timeidx < 1)+(timeidx > size(omega,2)));
            timeidx(NotInVid) = 1;%remove indices that fall before or after video
            
            slice = squeeze(omega(implant,iModelTime,timeidx,iRDM));
            slice(NotInVid) = NaN;%remove indices that fall before or after video
            rstack(implant,iRDM,iModelTime,:) = slice;
            
        end
    end
end

% Average over video time
rstackLine = squeeze(nanmean(rstack,3));

%% compute and save slopes of simulated autocorrelation (or in this case PCAandPCR dRSA lag-plots resulting from simulated data)
% NOTE: this should be run on the results of the simulations using PCR
% for itime = 1:size(rstackLine,3)
%     modelautocorr_slopes(:,itime) = diag(squeeze(rstackLine(1:end-1,:,itime)));
% end
% save(fullfile(indirSIM,'..','..','modelautocorr_slopes'),'modelautocorr_slopes');

%% compute and save regression borders, i.e., at which lag the model itself should be regressed out as well to attenuate model autocorrelation effects
% NOTE: this should be run on the results of the simulations using simple correlation
% zeroID = dsearchn(TimeVec2',0);
% regborder = zeros(size(rstackLine,2),2);
% for iRDM = 1:size(rstackLine,2)
%     
%     peak = max(rstackLine(iRDM,iRDM,:));
%     
%     % based on 50%, 25%, 32% or 71% of peak of corr = 1 (i.e., 25%, 6.25%, 10% or 50% of shared variance)
%     regborder(iRDM,1) = ceil(nanmean([find(flip((squeeze(rstackLine(iRDM,iRDM,1:zeroID)))) < .5*peak,1) find(squeeze(rstackLine(iRDM,iRDM,zeroID:end)) < .5*peak,1)]));
%     regborder(iRDM,2) = ceil(nanmean([find(flip((squeeze(rstackLine(iRDM,iRDM,1:zeroID)))) < .25*peak,1) find(squeeze(rstackLine(iRDM,iRDM,zeroID:end)) < .25*peak,1)]));
%     regborder(iRDM,3) = ceil(nanmean([find(flip((squeeze(rstackLine(iRDM,iRDM,1:zeroID)))) < .3162*peak,1) find(squeeze(rstackLine(iRDM,iRDM,zeroID:end)) < .3162*peak,1)]));
%     regborder(iRDM,4) = ceil(nanmean([find(flip((squeeze(rstackLine(iRDM,iRDM,1:zeroID)))) < .71*peak,1) find(squeeze(rstackLine(iRDM,iRDM,zeroID:end)) < .71*peak,1)]));  
%     
% end
% regborder(isnan(regborder)) = 100; 

% save(fullfile(indirSIM,'..','..','regressionBorderPerModel_smRDM20msec'),'regborder');

%% figures
close(figure(1));

figure(1);
set(gcf,'color','w');
set(gcf, 'Units', 'centimeters');
set(gcf, 'Position', [1 5 24 3]);

for iRDM = 1:length(cfg.models2test)
    subplot(1,8,iRDM);
    
    h = [];
    
    hold on;
    h(1) = plot([TimeVec2(1) TimeVec2(end)],[0 0],'--k');
    for implant = 1:size(rstackLine,1)-1
        h(implant+1) = plot(TimeVec2,squeeze(rstackLine(implant,cfg.models2test(iRDM),:)),'Color',colors(implant,:),'lineWidth',1.5);
    end
    h(end+1) = plot(TimeVec2,squeeze(rstackLine(end,cfg.models2test(iRDM),:)),'--','Color',[.2 .2 .2],'lineWidth',1.5);
        
    plot([0 0],[-1 1.5],'--k');
    plot([cfg.lag cfg.lag],[-1 1],'--k');
    set(gca,'xlim',[TimeVec2(1) TimeVec2(end)]);
    set(gca,'ylim',[-0.1 1]);
    set(gca,'xtick',TimeVec2(binVec2));
    set(gca,'ytick',[0 .5 1 1.5]);
    if iRDM ~= 1
        set(gca,'yticklabel','');
    end
    set(gca,'Fontsize',6,'FontName','Helvetica');
    hold off
    
    if iRDM == 1
        if cfg.glmRSA == 0
            ylabel('dRSA [corr]'); 
        elseif cfg.glmRSA == 1
            ylabel('dRSA [beta]'); 
        end
    end
    
end% iRDM

% cd('XXX\ActionPrediction\figures');
% print -depsc -r600 ActionPrediction_FigureS4a_simulations.eps
% print -depsc -r600 ActionPrediction_FigureS4b_simulations.eps

% % set(gcf,'renderer','Painters')
% print -depsc -tiff -r600 -painters ActionPrediction_FigureS4a_simulations.eps
% print -depsc -tiff -r600 -painters ActionPrediction_FigureS4b_simulations.eps

end