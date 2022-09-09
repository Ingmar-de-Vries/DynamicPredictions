function DynamicPredictions_DynamicModelRDMs_eyeTracker(~,eyeSub,~)

% path to preprocessed eyetracker data
indirEye = '//XXX/ActionPrediction/data/MEG/ELpreprocessed';

%% load and prepare data
fn= sprintf('%s%cSUB%02d',indirEye, filesep, eyeSub);
load(fn,'data');

nstim = 14;
eyePos = data;
clear data

% only select eye positions
eyePos.trial = eyePos.trial(:,[1 2 4 5],:);

% select video time
toi = dsearchn(eyePos.time',[0 5]');
toi = toi(1):toi(end);
eyePos.trial = eyePos.trial(:,:,toi);
eyePos.time = eyePos.time(toi);

%% interpolate missing data due to blinks
temp = eyePos.trial;
ntrials = size(temp,1);
ntime = size(temp,3);

% concatenate all trials to do gap search in one go
temp = reshape(permute(temp,[2 3 1]),size(temp,2),ntrials*ntime);

% find data points with value zero or nan
missingdata = find(logical(sum(temp==0 + isnan(temp))));

% find start and end of gaps in data
gapends = [find(diff(missingdata)~=1) length(missingdata)];
gapstarts = [1 gapends(1:end-1)+1];

% add some padding to the gaps, because e.g., for a blink the signal might already reflect the blink at the gaps edges
padding = 5;% in samples
paddedgaps = [];
for igap = 1:length(gapstarts)  
    paddedgaps = [paddedgaps missingdata(gapstarts(igap))-padding:missingdata(gapends(igap))+padding];
end% gap loop
paddedgaps(paddedgaps<1)=[];

% before interpolation, remove the gaps from the data (temp) and from the time indices (idx2interp)
idx = 1:size(temp,2);
idx2interp = idx;
idx2interp(paddedgaps) = [];
temp(:,paddedgaps) = [];

% now interpolate to new time indices without gaps
temp = interp1(idx2interp,temp',idx,'pchip','extrap')';

% separate trials again
eyePos.trial = permute(reshape(temp,size(temp,1),ntime,ntrials),[3 1 2]);

% average over repetitions of same video
temp = zeros(nstim,4,length(toi));
for istim = 1:nstim
    temp(istim,:,:) = squeeze(mean(eyePos.trial(eyePos.trialinfo == istim,:,:)));
end
eyePos = temp;

%% create dynamic RDM
RDMeyePOS = zeros(nstim,nstim,size(eyePos,3),size(eyePos,3));
for istim1 = 1:nstim
    for istim2 = 1:nstim
        for itime1 = 1:length(toi)
            for itime2 = 1:length(toi)
                
                % Euclidean distance using matlab's norm function
                RDMeyePOS(istim1,istim2,itime1,itime2) = norm(squeeze(eyePos(istim1,:,itime1))-squeeze(eyePos(istim2,:,itime2)));

            end
        end% time 1 loop
    end% stim 2 loop
end% stim 1 loop

% save correlation matrices
save([fn '_dynRDM_eyeTracker'],'RDMeyePOS');

end
