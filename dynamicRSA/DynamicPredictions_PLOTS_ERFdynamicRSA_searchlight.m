%% Prepare wholebrain source level dRSA values for cortical surface plot in Brainstorm
addpath(genpath('\\XXX\ActionPrediction\code'));
addpath(genpath('\\XXX\matlab_toolboxes\CoSMoMVPA-master'));
cfg.path = '\\XXX\ActionPrediction';

if cfg.glmRSA == 0
    corrORglm = 'corr';
elseif cfg.glmRSA == 1
    corrORglm = 'pcaANDpcr_75comps_';
end

% load measure you want to display
indir = fullfile(cfg.path,'data','MEG',['searchlight_' cfg.atlas],'RSA','statistics',[corrORglm num2str(cfg.randshuff(1)) 'iterations_' num2str(ceil(cfg.randshuff(2)*1000)) 'msec']);
fn = sprintf('%s%cSTATS_onesided_fisherz%d_%dHz_smMEG%d_smRDMneu%d_smRDMmod%d',outdir, filesep, cfg.fisherz, cfg.downsample, cfg.smoothingMSec, cfg.smoothNeuralRDM, cfg.smoothModelRDM);
load(fn,'atlasALL','omegaMean','pvals','FDR05','FDR01');

% this can be any colormap of your preference
% we save it here, and manually import it into the Brainstorm GUI as colormap:
% cmap = brewermap(128,'*RdYlBu');% 'YlGnBu' or 'RdPu'
% cmap = cmap(65:128,:);
% cmap = brewermap(128,'*RdBu');
% cmap = cmap(19:128,:);
% save('\\XXX\ActionPrediction\code\neuralDecoding\cmap4searchlight_BuRdX','cmap');
% scoutnum = length(atlas);

% load standard atlas of template subject of which we will change the colors
atlasBST = load('\\XXX\ActionPrediction\data\brainstorm_database\ActionPrediction\anat\@default_subject\tess_cortex_pial_low');
atlasID = contains(extractfield(atlasBST.Atlas,'Name'),cfg.atlas);
atlasBST = extractfield(atlasBST.Atlas(atlasID),'Scouts');
atlasBST = atlasBST{:};

idx2remove = contains(extractfield(atlasBST,'Label'),'Background');% in case of Schaeffer atlas there's a background...
atlasBST(idx2remove) = [];

% give each vertex value according to it's location in atlas
FDR = 1;
valuePERvertex = zeros(15002,size(omegaMean,2));
for ipeak = 1:size(omegaMean,2)
        
    if FDR
        % for FDR corrected, only color parcels that survive a certain threshold:
        roi2color = find(FDR05(:,ipeak))';% only those surviving
    else
        % for uncorrected maps, color all parcels:
        roi2color = 1:400;
    end

    for iroi = roi2color
        vertPERroi = extractfield(atlasBST(iroi),'Vertices');
        valuePERvertex(vertPERroi,ipeak) = omegaMean(iroi,ipeak);
    end
end

% Now a manual part starts:
% start Brainstorm GUI, and first export to your Matlab workspace from Brainstorm GUI a source map from the 'group analysis'...
% here I called this new variable 'temporary', see below.
% then replace values with seachlight values: 
temporary.ImageGridAmp = valuePERvertex;
temporary.Time = 1:size(valuePERvertex,2);
% and load this temporary back into Brainstorm GUI for creating pretty figures
% additionally in the Brainstorm GUI I manually applied spatial smoothing with 2mm FWHM
% I change the background to white, and the non-colored parcels to light gray
% and manually save the figures as images, which I combined in illustrator to create Figure 2b and S6