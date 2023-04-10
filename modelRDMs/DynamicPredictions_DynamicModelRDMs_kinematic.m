function DynamicPredictions_DynamicModelRDMs_kinematic(cfg,istim1,istim2)
% kinematic marker RDMs take a very long time to compute, but this only needs to be done once
% my approach here was to compute each pairwise comparison (14x14 = 196) on a different node on the cluster, and combine them afterwards into the RDMs
% outside of the cluster. 

% input and output directories
if isfield(cfg,'cluster')
    indirKin = '\\XXX\ActionPrediction\experiment\kinematics';
    DataDir = '\\XXX\ActionPrediction\data\modelRDMs';
else
    indirKin = '/XXX/ActionPrediction/experiment/kinematics';
    DataDir = '/XXX/ActionPrediction/data/modelRDMs';
end

% tNew and t4deriv are only used for interpolation of derivatives of kinematic models
fs = 100;
tNew = 0:1/fs:5-1/fs;% after jittering
t4deriv = tNew+(1/fs)/2;% timesteps in between for any derivative (e.g., position --> motion, or motion --> acceleration)
t4deriv(end) = [];

%load the kinematic data
kinNames = dir(fullfile(indirKin, '*.mat'));

seqs = zeros(length(kinNames),13,3,500);
for iseq=1:length(kinNames)
    
    load(fullfile(indirKin, kinNames(iseq).name),'Trajs');
    seqs(iseq,:,:,:) = Trajs;
    
end

% extract kinematic data sequences for current two stimuli
seq1 = squeeze(seqs(istim1,:,:,:));
seq2 = squeeze(seqs(istim2,:,:,:));

% posture vectors
pos1 = reshape(seq1,numel(seq1(:,:,1)),size(seq1,3));
pos2 = reshape(seq2,numel(seq2(:,:,1)),size(seq2,3));

% compute dynamic RDM of view dependent body posture
RDMposture_dep = 1 - corr(pos1,pos2);

% motion vectors
vel1 = diff(pos1,[],2)';
vel2 = diff(pos2,[],2)';

% interpolate because motion vectors now placed in between time points
vel1 = interp1(t4deriv,vel1,tNew,'pchip','extrap')';
vel2 = interp1(t4deriv,vel2,tNew,'pchip','extrap')';

% compute dynamic RDM of view dependent body motion
RDMmotion_dep = 1 - corr(vel1,vel2);

% acceleration vectors
acc1 = diff(vel1,[],2)';
acc2 = diff(vel2,[],2)';

% interpolate because acceleration vectors now placed in between time points
acc1 = interp1(t4deriv,acc1,tNew,'pchip','extrap')';
acc2 = interp1(t4deriv,acc2,tNew,'pchip','extrap')';

% compute dynamic RDM of view dependent body motion
RDMacc_dep = 1 - corr(acc1,acc2);

%% view invariant
% for each comparison between seq1 and seq2 at time1 and time2, we need to align first for view invariance, then compute motion and
% acceleration, then compute all RDMs. To do that we need to iteratively compute snippets of data of 9 samples long.
tNewShort = tNew(1:9);
t4derivShort = t4deriv(1:8);

% first add padding to sequence 1 and 2
seq1 = permute(seq1,[3 1 2]);
seq1 = interp1(1:500,seq1,-3:504,'linear','extrap');
seq1 = permute(seq1,[2 3 1]);
seq2 = permute(seq2,[3 1 2]);
seq2 = interp1(1:500,seq2,-3:504,'linear','extrap');
seq2 = permute(seq2,[2 3 1]);

% initialize
RDMposture_invar = zeros(size(RDMposture_dep));
RDMmotion_invar = RDMposture_invar;
RDMacc_invar = RDMposture_invar;
for itime1 = 5:size(seq1,3)-4
    
    for itime2 = 5:size(seq2,3)-4
        
        % 9-sample snippets of sequence for current (=C) iteration
        Cseq1 = seq1(:,:,itime1-4:itime1+4);
        Cseq2 = seq2(:,:,itime2-4:itime2+4);
        
        % Align sequence 2 to sequence 1 using constrained procrustes
        Cseq2trans = zeros(size(Cseq2));
        for itime2align = 1:size(Cseq1,3)
            [~, Cseq2trans(:,:,itime2align)] = procrustes_constrain_rotationZaxis_IdV(squeeze(Cseq1(:,:,itime2align)),squeeze(Cseq2(:,:,itime2align)),'Reflection',false,'Scaling',false);
        end
        
        % Position vectors
        Cpos1 = reshape(Cseq1,numel(Cseq1(:,:,1)),size(Cseq1,3));
        Cpos2trans = reshape(Cseq2trans,numel(Cseq2trans(:,:,1)),size(Cseq2trans,3));
        clear Cseq1 Cseq2trans
        
        % compute dynamic RDM of view invariant body posture
        RDMposture_invar(itime1-4,itime2-4) = 1 - corr(Cpos1(:,end/2+1/2),Cpos2trans(:,end/2+1/2));
        
        % motion vectors
        Cvel1 = diff(Cpos1,[],2)';
        Cvel2trans = diff(Cpos2trans,[],2)';
        clear Cpos1 Cpos2trans
        
        Cvel1 = interp1(t4derivShort,Cvel1,tNewShort,'pchip','extrap')';
        Cvel2trans = interp1(t4derivShort,Cvel2trans,tNewShort,'pchip','extrap')';
        
        % compute dynamic RDM of view invariant body motion
        RDMmotion_invar(itime1-4,itime2-4) = 1 - corr(Cvel1(:,end/2+1/2),Cvel2trans(:,end/2+1/2));
        
        % acceleration vectors
        Cacc1 = diff(Cvel1,[],2)';
        Cacc2trans = diff(Cvel2trans,[],2)';
        clear Cvel1 Cvel2trans
        
        Cacc1 = interp1(t4derivShort,Cacc1,tNewShort,'pchip','extrap')';
        Cacc2trans = interp1(t4derivShort,Cacc2trans,tNewShort,'pchip','extrap')';
        
        % compute dynamic RDM of view invariant body acceleration
        RDMacc_invar(itime1-4,itime2-4) = 1 - corr(Cacc1(:,end/2+1/2),Cacc2trans(:,end/2+1/2));
        clear Cacc1 Cacc2trans
    end% time of stimulus 2
    
end% time of stimulus 1

% save correlation matrices
save([DataDir filesep 'ActionPrediction_dynRDM_posturedep_stim' num2str(istim1) '_stim' num2str(istim2)],'RDMposture_dep');
save([DataDir filesep 'ActionPrediction_dynRDM_postureinvar_stim' num2str(istim1) '_stim' num2str(istim2)],'RDMposture_invar');
save([DataDir filesep 'ActionPrediction_dynRDM_motiondep_stim' num2str(istim1) '_stim' num2str(istim2)],'RDMmotion_dep');
save([DataDir filesep 'ActionPrediction_dynRDM_motioninvar_stim' num2str(istim1) '_stim' num2str(istim2)],'RDMmotion_invar');
save([DataDir filesep 'ActionPrediction_dynRDM_accdep_stim' num2str(istim1) '_stim' num2str(istim2)],'RDMacc_dep');
save([DataDir filesep 'ActionPrediction_dynRDM_accinvar_stim' num2str(istim1) '_stim' num2str(istim2)],'RDMacc_invar');

% now combine each pairwise comparison into the single RDM outside of the cluster, i.e., locally
if ~isfield(cfg,'cluster')

    RDMposture_depAll = zeros(size(seqs,1),size(seqs,1),size(seqs,4),size(seqs,4));
    RDMposture_invarAll = RDMposture_depAll;
    RDMmotion_depAll = RDMposture_depAll;
    RDMmotion_invarAll = RDMposture_depAll;
    RDMacc_depAll = RDMposture_depAll;
    RDMacc_invarAll = RDMposture_depAll;
    for istim1 = 1:size(seqs,1)
        for istim2 = 1:size(seqs,1)
                        
            load([DataDir filesep 'ActionPrediction_dynRDM_posturedep_stim' num2str(istim1) '_stim' num2str(istim2)],'RDMposture_dep');
            load([DataDir filesep 'ActionPrediction_dynRDM_postureinvar_stim' num2str(istim1) '_stim' num2str(istim2)],'RDMposture_invar');
            load([DataDir filesep 'ActionPrediction_dynRDM_motiondep_stim' num2str(istim1) '_stim' num2str(istim2)],'RDMmotion_dep');
            load([DataDir filesep 'ActionPrediction_dynRDM_motioninvar_stim' num2str(istim1) '_stim' num2str(istim2)],'RDMmotion_invar');
            load([DataDir filesep 'ActionPrediction_dynRDM_accdep_stim' num2str(istim1) '_stim' num2str(istim2)],'RDMacc_dep');
            load([DataDir filesep 'ActionPrediction_dynRDM_accinvar_stim' num2str(istim1) '_stim' num2str(istim2)],'RDMacc_invar');
            
            RDMposture_depAll(istim1,istim2,:,:) = RDMposture_dep;
            RDMposture_invarAll(istim1,istim2,:,:) = RDMposture_invar;
            RDMmotion_depAll(istim1,istim2,:,:) = RDMmotion_dep;
            RDMmotion_invarAll(istim1,istim2,:,:) = RDMmotion_invar;
            RDMacc_depAll(istim1,istim2,:,:) = RDMacc_dep;
            RDMacc_invarAll(istim1,istim2,:,:) = RDMacc_invar;
        end
    end
    
    RDMposture_dep = RDMposture_depAll;
    RDMposture_invar = RDMposture_invarAll;
    RDMmotion_dep = RDMmotion_depAll;
    RDMmotion_invar = RDMmotion_invarAll;
    RDMacc_dep = RDMacc_depAll;
    RDMacc_invar = RDMacc_invarAll;
    
    save([DataDir filesep 'ActionPrediction_dynRDM_posturedep'],'RDMposture_dep');
    save([DataDir filesep 'ActionPrediction_dynRDM_postureinvar'],'RDMposture_invar');
    save([DataDir filesep 'ActionPrediction_dynRDM_motiondep'],'RDMmotion_dep');
    save([DataDir filesep 'ActionPrediction_dynRDM_motioninvar'],'RDMmotion_invar');
    save([DataDir filesep 'ActionPrediction_dynRDM_accdep'],'RDMacc_dep');
    save([DataDir filesep 'ActionPrediction_dynRDM_accinvar'],'RDMacc_invar');
        
end


end