clear
format compact
clc

% specify a filename 
save_file = 'MCMC_results.mat';


%% get data
load('example_data.mat', 't', 'm', 's', 'K', 't_min', 't_max', 'units_time', 'units_conc', 'ground');


%% prepare input
opts.units_time = units_time;
opts.units_conc = units_conc;
opts.t_min = t_min;
opts.t_max = t_max; 
opts.t = t;
opts.K = K;

% opts.z = m; % means (sets linear constraints)
opts.z = [m,s]; % means and stds (sets linear and quadratic constraints)

if exist('ground','var')
    opts.ground = ground;
end


%% init chain
chain = chainer_main([],0,opts,true,[]);
clear opts % clear leftovers

save(save_file,'chain','save_file') % save chain for future processing
disp(['SAVED: ', save_file])


%% expand chain in batches
mcmc_batch_size = 150;

while true
    
    % stop expanding, if chain is too long
    if chain.i(end) >= 1200
        break
    end

    % reset acceptance record
    chain = chainer_main(chain,[],[],true,[]);
    
    % expand chain
    chain = chainer_main(chain,mcmc_batch_size,[],true,true); % with active visuals (too slow)
    % chain = chainer_main(chain,mcmc_batch_size,[],true,false); % without active visuals (much faster)
    
    % thin chain, if chain is too big
    if chain.sizeGB > 1.0
        chain = chainer_main(chain,-fix(chain.length/2),[],true,[],[]); 
    end
    
    % save chain for future processing
    save(save_file,'chain','save_file') 
    disp(['SAVED: ', save_file])

end    


%% get results
clf

% remove burn-in
idx = find(chain.i>0.2*chain.i(end));

% show posterior
plotmatrix(chain.g(:,idx)')

% get map
[~,i] = max(chain.P(1,:));
g_map = chain.g(:,i);

