function [stat,signpos,signneg] = DynamicPredictions_runFTstats(powspctrm,addinfo,template,testtype,pthresh,pthreshclust,baseline)
% shell around Fieldtrip functions for statistical testing, where parameters are set for our specific case
% this function is called from "DynamicPredictions_STATS_ERFdynamicRSA_ROIsource.m", which in turn is called from "DynamicPredictions_pipeline.m"

subnum = size(powspctrm,1);

if contains(testtype,'uni')
    
    % we need some channel file
    ft_cfg = [];
    ft_cfg.template = template;
    ft_cfg.channel = addinfo.label;
    ft_cfg.method = 'template';
    neighbours = ft_prepare_neighbours(ft_cfg);
    
elseif contains(testtype,'multi')
    neighbours = [];
    
    %duplicate as if multiple sensors to trick Fieldtrip into thinking it's the correct format
    addinfo.label = {'1','2'};
    if contains(testtype,'timefreq')
        powspctrm = repmat(powspctrm,[1 1 1 2]);
        powspctrm = permute(powspctrm,[1 2 4 3]);
    elseif contains(testtype,'timeseries')
        powspctrm = repmat(powspctrm,[1 1 2]);
        powspctrm = permute(powspctrm,[1 3 2]); 
    end
end

% set parameters that are the same for all testtypes
% see Fieldtrip documentation for details on these settings
ft_cfg = [];
ft_cfg.neighbours = neighbours;
ft_cfg.method = 'montecarlo';
ft_cfg.statistic = 'depsamplesT';
ft_cfg.correctm = 'cluster';
ft_cfg.clusteralpha = pthreshclust;
ft_cfg.clusterstatistic = 'maxsum';
ft_cfg.minnbchan = 1;
ft_cfg.tail = addinfo.tail;
ft_cfg.clustertail = addinfo.tail;
ft_cfg.alpha = pthresh;

% the more the better, but for the time-time plots anything more than 5000 takes a very long time
if contains(testtype,'timefreq')
    ft_cfg.numrandomization = 5000;
elseif contains(testtype,'timeseries')
    ft_cfg.numrandomization = 25000;
end

design = zeros(2,2*subnum);
for i = 1:subnum
    design(1,i) = i;
end
for i = 1:subnum
    design(1,subnum+i) = i;
end
design(2,1:subnum)        = 1;
design(2,subnum+1:2*subnum) = 2;

ft_cfg.design = design;
ft_cfg.uvar  = 1;
ft_cfg.ivar  = 2;

%now set testtype specific parameters, and run statistical test
if contains(testtype,'timeseries')% let fieldtrip average over the sensors, so test is only run over the dimension time
    
    ft_cfg.avgoverchan = 'yes';
    
    %change data to fieldtrip stats format
    ss.dimord = 'chan_time';
    ss.label = addinfo.label;
    ss.time = addinfo.time;
    ss.fsample = 1/((ss.time(end)-ss.time(1))/(length(ss.time)-1));

    data4stat = cell(1,subnum);
    for isub=1:subnum
        data4stat{isub} = ss;
        data4stat{isub}.avg = squeeze(powspctrm(isub,:,:));
    end
    
    %create datastruct filled with baseline value to test observed values against (zero in case of dRSA)
    datazero = data4stat;
    for isub=1:subnum
        datazero{isub}.avg = ones(size(datazero{isub}.avg)).*baseline; 
    end
     
    %run statistical test
    stat = ft_timelockstatistics(ft_cfg, data4stat{:},datazero{:});
    
    %initialize
    signpos = false(1,length(addinfo.time));
    signneg = signpos;
    
elseif contains(testtype,'timefreq')% let fieldtrip average over the sensors, so test is only run over the dimensions time and frequency (or time and time for 2D dRSA)
    
    ft_cfg.avgoverchan = 'yes';
    ft_cfg.frequency        = 'all';
    
    %change data to fieldtrip stats format
    data4stat.dimord = 'subj_chan_freq_time';
    data4stat.label = addinfo.label;
    data4stat.freq = addinfo.freq;
    data4stat.time = addinfo.time;
    data4stat.powspctrm = permute(powspctrm,[1 3 2 4]);

    %create datastruct filled with baseline value to test observed values against (zero in case of dRSA)
    datazero = data4stat;
    datazero.powspctrm = ones(size(datazero.powspctrm)).*baseline;  
    
    %run statistical test
    stat = ft_freqstatistics(ft_cfg, data4stat,datazero);
    
    %initialize
    signpos = false(length(addinfo.freq),length(addinfo.time));
    signneg = signpos;

end
    
% get relevant (significant) values, cluster-corrected
if isfield(stat,'posclusters')
    if isfield(stat.posclusters,'prob')
        pos_cluster_pvals = [stat.posclusters(:).prob];
        pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);% find positive clusters that survive cluster-correction
        signpos = squeeze(ismember(stat.posclusterslabelmat, pos_signif_clust));
    end
end
if isfield(stat,'negclusters')
    if isfield(stat.negclusters,'prob')
        neg_cluster_pvals = [stat.negclusters(:).prob];
        neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);% find negative clusters that survive cluster-correction
        signneg = squeeze(ismember(stat.negclusterslabelmat, neg_signif_clust));
    end
end

end