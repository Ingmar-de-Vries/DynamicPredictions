function DynamicPredictions_STATS_ERFdynamicRSA_searchlight(cfg,~,~)

% You can run this function locally (by looping over subjects), or send to a cluster (which I did using the accompanying "cluser_shell.m" function)
% set some directories depending on running this locally or on a cluster:
if isfield(cfg,'cluster')
    addpath('//XXX/ActionPrediction/toolboxes/fieldtrip-20191113');
    addpath(genpath('//XXX/ActionPrediction/code'));
    addpath(genpath('//XXX/matlab_toolboxes/CoSMoMVPA-master'));
    cfg.path = '//XXX/ActionPrediction'; 
else
    addpath('\\XXX\ActionPrediction\toolboxes\fieldtrip-20191113');
    addpath(genpath('\\XXX\ActionPrediction\code'));
    addpath(genpath('\\XXX\matlab_toolboxes\CoSMoMVPA-master'));
    cfg.path = '\\XXX\ActionPrediction';
end
ft_defaults

if cfg.glmRSA == 0
    corrORglm = 'corr';
elseif cfg.glmRSA == 1
    corrORglm = ['pcaANDpcr_' num2str(cfg.nPCAcomps) 'comps'];
end

%% input and output folders
indir = fullfile(cfg.path,'data','MEG',['searchlight_' cfg.atlas],'RSA',[corrORglm num2str(cfg.randshuff(1)) 'iterations_' num2str(ceil(cfg.randshuff(2)*1000)) 'msec']);
outdir = fullfile(cfg.path,'data','MEG',['searchlight_' cfg.atlas],'RSA','statistics',[corrORglm num2str(cfg.randshuff(1)) 'iterations_' num2str(ceil(cfg.randshuff(2)*1000)) 'msec']);
if ~exist(outdir,'dir')
    mkdir(outdir);
end

%% load data, combine subjects
if strcmp(cfg.atlas,'HCP')
    numsource = 360;
elseif contains(cfg.atlas,'Schaefer2018_400')
    numsource=400;
elseif contains(cfg.atlas,'Schaefer2018_200')
    numsource=200;
end

atlasALL = cell(numsource,length(cfg.sub4stat));
omegaGA = zeros(length(cfg.sub4stat),numsource,length(cfg.peaks2test));
for isub = cfg.sub4stat
    
    fn = sprintf('%s%cdRSA_SUB%02d_%dHz_smMEG%d_smRDMneu%d_smRDMmod%d',indir, filesep, isub, cfg.downsample, cfg.smoothingMSec, cfg.smoothNeuralRDM, cfg.smoothModelRDM);
    fprintf('loading %s..\n',fn);
    
    load(fn,'omegaAll','atlas');
    
    omegaGA(isub,:,:) = omegaAll;
    atlasALL(:,isub) = extractfield(atlas,'Label'); 
    
    % used to check other subjects against to ensure atlases align across subjects
    if isub == 1
        atlas2check = extractfield(atlas,'Label');
    end
    
    % fix subject 12 atlas such that it matches the one from all other subject
    if isub == 12
        sub12 = extractfield(atlas,'Label')';

        % first order according to hemisphere (left -- right) 
        lefties = false(size(sub12));
        for iroi = 1:length(atlas)
            if strcmp(sub12{iroi}(end),'L')
                lefties(iroi) = true;
            end
        end
        atlas = [atlas(lefties) atlas(~lefties)];% and the actual atlas
%         sub12 = [sub12(lefties) ; sub12(~lefties)];% temporary atlas labels
        
%         % then manually change the mismatching names:
%         new2oldSchaefer200 = [1:5 9 6:8 10:19 24:26 20:23 27:42 44:48 43 49:53 58 59 54:57 60 61 64 62 63 65:67 69 70 68 71:80 84 85 81 86 82 ...
%             83 87:90 94 91:93 95:97 100 98 99 101:110 112:114 111 115:120 125 121:124 126:132 134:138 133 139:143 148:151 144:147 152 157 158 ...
%             153 155 156 154 159:161 165 164 166 162 163 167:179 181:184 180 185:190 194 191:193 195 196 200 197:199];
%         load('\\XXX\ActionPrediction\data\MEG\searchlight_Schaefer2018_200\new2oldSchaefer200','new2oldSchaefer200');

        % then manually change the mismatching names:
        new2oldSchaefer400 = [1:30 40:46 31:39 47:100 108:112 101:107 113:115 118:121 116 117 122:127 131:135 128:130 136:155 166 167 159:163 ...
            156:158 164 165 168:170 171:177 187 178 179 180:184 188 185 186 189 196 197 198 190:193 199 200 194 195 201:235 245:249 236:244 ...
            250:266 269:280 267 268 281:292 299:304 293:298 312:315 309:311 305 316 317 306 318 319 307 308 320:324 330 327:329 331 332 325 ...
            326 333:353 365:367 356:358 354 359:361 355 362:364 368:378 388 379 380:385 389 386 387 390 397 398 391:394 399 400 395 396];
% %         load('\\XXX\ActionPrediction\data\MEG\searchlight_Schaefer2018_400\new2oldSchaefer400','new2oldSchaefer400');
        
        % reorder data, first left-right, then remaining mismatch 
        omegaAll = [omegaAll(lefties,:) ; omegaAll(~lefties,:)];
        omegaAll = omegaAll(new2oldSchaefer400,:); 
        omegaGA(isub,:,:) = omegaAll;
        atlas = atlas(new2oldSchaefer400);% change atlas names 
        atlasALL(:,isub) = extractfield(atlas,'Label'); 
    end
    
    % check whether atlases match between subjects, with atlas of subject 1 being the reference
    if any(~strcmp(atlas2check,extractfield(atlas,'Label'))) && isub ~= 12% except for 12 because we know it doesn't match, but we fixed it manually above
        error(['The atlas of subject ' num2str(isub) ' does not match the atlas of subject 1'])
    end
    
end% subject loop

if cfg.fisherz
    omegaGA = atanh(omegaGA);
end

%% Statistics
fnSTATS = sprintf('%s%cSTATS_onesided_fisherz%d_%dHz_smMEG%d_smRDMneu%d_smRDMmod%d',outdir, filesep, cfg.fisherz, cfg.downsample, cfg.smoothingMSec, cfg.smoothNeuralRDM, cfg.smoothModelRDM);

pvals = zeros(size(omegaGA,2),size(omegaGA,3));
tvals = zeros(size(omegaGA,2),size(omegaGA,3));
FDR05 = pvals;
FDR01 = pvals;
for ipeak = 1:size(omegaGA,3)
    
    data2test = squeeze(omegaGA(:,:,ipeak));
    
    if isnan(min(min(data2test))) && isnan(max(max(data2test)))
        
        pvals(:,ipeak) = NaN(size(omegaGA,2),1);
        FDR05(:,ipeak) = NaN(size(omegaGA,2),1);
        FDR01(:,ipeak) = NaN(size(omegaGA,2),1);
    else
        
        [~,p,~,stats] = ttest(data2test,0,'tail','right');
        tvals(:,ipeak) = stats.tstat;
        pvals(:,ipeak) = p;
        [FDR05(:,ipeak),~,~,~] = fdr_bh(p,.05);
        [FDR01(:,ipeak),~,~,~] = fdr_bh(p,.01);
    end% if loop
    
end% peak loop

omegaMean = squeeze(mean(omegaGA));
save(fnSTATS,'omegaMean','atlasALL','tvals','pvals','FDR05','FDR01');

end
