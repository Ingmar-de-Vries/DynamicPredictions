function FTdata = DynamicPredictions_MEGpp10_realign2photodiode(FTdata)
%% ===== load FT sensor data and realign to photodiode =====

% find photodiode delay
photochan = strcmp(FTdata.label,'MISC008');
sampledelay = zeros(length(FTdata.trial),1);
for itrial = 1:length(FTdata.trial)
    tzero = dsearchn(FTdata.time{itrial}',0);
    
    for isample = tzero:length(FTdata.time{itrial})
        if FTdata.trial{itrial}(photochan,isample) > 0.05
            sampledelay(itrial) = isample-tzero;
            break
        end
    end
end

% now realign such that photodiode comes 2 samples after zero, as is
% the case in most trials
% By far most trials have a delay of 2 samples in the onset of the photodiode.
% We take that as the zero point, because it actually takes 2 samples for the photodiode to reach our threshold of >0.05 (which can be seen by plotting the photodiode signal).
% Also, this way as little trials as possible will actually be realigned / changed, which is preferrable
sampledelay = sampledelay - 2;

% Shift all data forward according to the sample delay. Of course this
% results in some useless samples at the end of the trial, but that
% doesn't matter as they're only their to prevent edge artifacts during
% TF analyses anyway.
for itrial = 1:length(FTdata.trial)
    if sampledelay(itrial) > 0
        FTdata.trial{itrial}(:,1:end-sampledelay(itrial)+1) = FTdata.trial{itrial}(:,sampledelay(itrial):end);
    end
end


end