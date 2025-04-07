function sample = chainer_init_sample(params,opts)


% % start from ground truth
% sample.g = params.ground.g;
% sample.h = params.ground.h;
% sample.y = params.ground.y;
% %


%% counter
sample.i = 0;


%% model variables
sample.g = params.g_prior_psi./params.g_prior_phi .* randg(params.g_prior_phi);
sample.h = params.h_prior_psi /params.h_prior_phi  * randg(params.h_prior_phi);
sample.y = get_initial_y(sample,params);
init_y_count = 1;
disp('Initializing')
while max(abs(get_c(reshape(sample.y',params.N*params.K,1),params)),[],'all') > params.initial_constraint_tolerance ...
    || min(sample.y,[],'all') < 0
    sample.g = params.g_prior_psi./params.g_prior_phi .* randg(params.g_prior_phi);
    sample.h = params.h_prior_psi /params.h_prior_phi  * randg(params.h_prior_phi);
    sample.y = get_initial_y(sample,params);
    % y_error = max(abs(get_c(reshape(sample.y',params.N*params.K,1),params)),[],'all')
    init_y_count = init_y_count + 1;
    if init_y_count > 1000
        error('Initial sample cannot satisfy constraints')
    end
end

%% book-keeping
sample.P = get_log_probs(sample.g,sample.h,sample.y,params);

sample.rec = repmat([0;realmin],1,2);
end

function y = get_initial_y(sample,params)
    x = get_x(sample.g,params.t,params.t_min,params.t_max);
    y = x(:,2) .* exp( randn(params.N,params.K) / sqrt(sample.h) );
    % log_y = log(y);
    % log_y = reshape(log_y',params.N*params.K,1);
    % log_y = fsolve(@(log_y_c) get_c(exp(log_y_c),params), log_y, params.fsolve_options);
    % y = reshape(exp(log_y),params.K,params.N)';
    y = reshape(y',params.N*params.K,1);
    y = fsolve(@(y_c) get_c(y_c,params), y, params.fsolve_options);
    y = reshape(y,params.K,params.N)';
end


