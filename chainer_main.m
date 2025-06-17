function chain = chainer_main(chain_init,d_length,opts,flag_status,flag_visual)
% to init:
% chain = chainer_main([]   ,  0, opts, true, []  );
% to expand:
% chain = chainer_main(chain,+25, []  , true, true);
% to reduce:
% chain = chainer_main(chain,-10, []  , true, []  );
% to reset record:
% chain = chainer_main(chain, [], []  , true, []  );

tic_id = tic;


% initialize the seed or use seed for last expansion
if isempty(chain_init)
    rng('shuffle');
else
    rng(chain_init.random);
end


% init chain --------------------------------------------------------------
if d_length == 0

    % MCMC
    chain.params = chainer_init_params(opts);
    % chain.params = SA_chainer_init_params(opts);
    chain.length = 1;
    chain.stride = 1;
    chain.ledger = nan(0,2); % total wall time without Initialization
    chain.sizeGB = nan;      % current chain memory size
    chain.record = [];       % acceptance rates
    chain.sample = [];

    chain.sample = chainer_init_sample(chain.params,opts);
    
    % history
    chain.i = cast( chain.sample.i, 'uint64');
    chain.P = cast( chain.sample.P, 'double');
    chain.g = cast( chain.sample.g, 'single');
    chain.h = cast( chain.sample.h, 'single');
    chain.y = cast( chain.sample.y, 'single');
    
    if flag_status
        disp('CHAINER: chain initiated')
    end

    

% expand chain ------------------------------------------------------------
elseif d_length > 0

    chain.params = chain_init.params;
    chain.length = chain_init.length + d_length;
    chain.stride = chain_init.stride;
    chain.ledger = chain_init.ledger;
    chain.sizeGB = nan;
    chain.record = chain_init.record;
    chain.sample = chain_init.sample;
    
    chain.i = cat( 2 , chain_init.i, zeros( 1                             , d_length , 'like',chain_init.i ) );
    chain.P = cat( 2 , chain_init.P,   nan( get_log_probs                 , d_length , 'like',chain_init.P ) );
    chain.g = cat( 2 , chain_init.g,   nan( 4                             , d_length , 'like',chain_init.g ) );
    chain.h = cat( 2 , chain_init.h,   nan( 1                             , d_length , 'like',chain_init.h ) );
    chain.y = cat( 3 , chain_init.y,   nan( chain.params.N, chain.params.K, d_length , 'like',chain_init.y ) );
    
    
    if flag_visual
        Gim = chainer_visualize([],chain);
    end
    
    %---------------------------- expand chain
    r = chain_init.length+1;
    while r <= chain.length
        
        chain.sample = sampler_update(chain.sample,chain.params);
        
        if mod(chain.sample.i,chain.stride) == 0
            
            chain.i(    r) = chain.sample.i;
            chain.P(:,  r) = chain.sample.P;
            chain.g(:,  r) = chain.sample.g;
            chain.h(    r) = chain.sample.h;
            chain.y(:,:,r) = chain.sample.y;
            
            if flag_visual
                chainer_visualize(Gim,chain);
            end
            
            if flag_status
                disp([  'i = ', num2str(chain.sample.i,'%d'), ...
                     ' - acc = ', ...
                                num2str( chain.sample.rec(1,:)./chain.sample.rec(2,:)  * 100 ,'%#6.2f') , ' %', ...
                     ])
            end
            
            r = r+1;
        end
    end    

    if flag_status
        disp('CHAINER: chain expanded')
    end


% reduce chain ------------------------------------------------------------
elseif d_length < 0

    d_length = min(-d_length,chain_init.length);

    chain.params = chain_init.params;
    chain.length = d_length;
    chain.stride = nan;
    chain.ledger = chain_init.ledger;
    chain.sizeGB = nan;
    chain.record = chain_init.record;
    chain.sample = chain_init.sample;
    
    ind = mod(chain_init.length,d_length)+(floor(chain_init.length/d_length)*(1:d_length));

    chain.i = chain_init.i(    ind);
    chain.P = chain_init.P(:  ,ind);
    chain.g = chain_init.g(:  ,ind);
    chain.h = chain_init.h(    ind);
    chain.y = chain_init.y(:,:,ind);
    
    chain.stride = double(chain.i(2)-chain.i(1));
    
    if flag_status
        disp('CHAINER: chain reduced')
    end
    
    
% reset chain -------------------------------------------------------------
elseif isempty(d_length)

    chain = chain_init;
    
    chain.record = [chain.record; [chain.sample.i,chain.sample.rec(1,:)./chain.sample.rec(2,:)] ];

    chain.sample.rec(1,:) = 0;
    chain.sample.rec(2,:) = realmin;
    chain.sample.rec(3,:) = realmin;

    if flag_status
        disp('CHAINER: chain reset')
    end
    
end


% store the seed for future expansion
chain.random = rng();



%% book-keeping
chain.sizeGB = get_sizeGB(chain);               % mem size

% ledger
wall_time = toc(tic_id);
chain.ledger = [chain.ledger; double(chain.i(end)), wall_time];

if flag_status
    disp(['( wall time = ',num2str(wall_time),' s, overall wall time = ',num2str(sum(chain.ledger(:,2))),' s )'])
end



end





%% auxiliary functions

function sizeGB = get_sizeGB(chain)
    sizeGB = whos( inputname(1) );
    sizeGB = sizeGB.bytes/1024^3;
end

