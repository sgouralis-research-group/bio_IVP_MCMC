function sample = sampler_update( sample, params )

%% update counter
sample.i = sample.i + 1;

%% No random sweep since only 2 cases
  
tag = 1;
[sample.h,sample.g,sample.rec(:,tag)] ...
    = update_hg(sample.g,sample.rec(:,tag),sample.y,params);
    
tag = 2;
[sample.h,sample.y,sample.rec(:,tag)] ...
    = update_hy(sample.y,sample.rec(:,tag),sample.g,params);
    

%% book-keeping
sample.P = get_log_probs(sample.g,sample.h,sample.y,params);

    
        