clear
format compact
clc

% Specify dir and file ending 
save_file = ['Runs/',mfilename,'_1.mat'];

%% get data
load('synthetic_data.mat', 't', 'z', 'v', 'K', 't_min', 't_max', 'units_time', 'units_conc', 'ground');

%% prepare input
opts.units_time = units_time;
opts.units_conc = units_conc;
opts.t_min = t_min;
opts.t_max = t_max; 
opts.t = t;
opts.z = z;
opts.v = v;
opts.K = K;

if exist('ground','var')
    opts.ground = ground;
end


%% init chain
chain = chainer_main([],0,opts,true,[]);
clear opts % clear leftovers
chain_temp = chain; % store a temporary copy for debugging 
save(save_file,'chain','save_file') % save chain for future processing
disp(['SAVED: ', save_file])

%% expand chain
run_expander


