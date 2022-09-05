function cluster_shell_simulations(cfg,script2run)
% determine whether to send job to cluster, or run script locally

if isfield(cfg,'cluster')% send to cluster
    
    matlab_version = sscanf(version('-release'),'%d%s');
    if matlab_version(1) ~= 2020
        error('The package for sending jobs to the cluster only works on matlab 2020, please use that');
    end

    cluster = parcluster;
    job = createJob(cluster);
    
    % Assign tasks using your function
    for iterbatch = 1:cfg.iterbatches% to speed up process, run different iterations in parallel on different nodes on the cluster, and average over those afterwards
        for iroi=cfg.ROIVec
            for isub = cfg.SubVec
                createTask(job,str2func(script2run), 0,{cfg, isub, iroi, iterbatch});%
            end
        end
    end
    % run task
    submit(job);
    
else% run locally
    
    feval(script2run,cfg,0,0);%if run locally, looping over subjects and ROIs happens within function (saves time as you only have to load the data of a subject once for all ROIs)
    
end