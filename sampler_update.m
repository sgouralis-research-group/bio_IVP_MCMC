function sample = sampler_update( sample, params )

%% update counter
sample.i = sample.i + 1;

%% No random sweep since only 2 cases
  
[sample.h,sample.g,sample.rec(:,1)] = update_hg( ...
          sample.g,sample.rec(:,1), sample.y,params);
    
[sample.h,sample.y,sample.rec(:,2)] = update_hy( ...
          sample.y,sample.rec(:,2), sample.g,params);
    

%% book-keeping
sample.P = get_log_probs(sample.g,sample.h,sample.y,params);