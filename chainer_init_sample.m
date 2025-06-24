function sample = chainer_init_sample(params,opts)


%% counter
sample.i = 0;


%% model variables
sample.g = params.g_prior_psi./params.g_prior_phi .* randg(params.g_prior_phi);
sample.h = params.h_prior_psi /params.h_prior_phi  * randg(params.h_prior_phi);
sample.y = get_initial_y(sample,params);


%% book-keeping
sample.P = get_log_probs(sample.g,sample.h,sample.y,params);

sample.rec = repmat([0;realmin],1,2);
end

%% optimizer
function y = get_initial_y(sample,params)
    x = get_x(sample.g,params.t,params.t_min,params.t_max);
    y = x(:,2) .* exp( randn(params.N,params.K) / sqrt(sample.h) );
    y = reshape(y',params.N*params.K,1);
    y = abs(fsolve(@(y_c) get_c(abs(y_c),params), y, params.fsolve_options));
    y = reshape(y,params.K,params.N)';
    for i = 1:50
        [~,y,~] = update_hy(y, nan(3,1), sample.g, params);
    end
end


