function log_P = get_log_probs(g,h,y,params)
% log_P = log( post,like,priors )

if nargin<1
    log_P = 5;
    return
end

log_P = nan(get_log_probs,1);


%% Likelihood
if max(abs(get_c(y',params))) < 100 * eps( max(params.Z(:,1)) )
    log_P(2) = 0;  % constraints satisfied within tolerance
else
    log_P(2) = -inf; % constraints not satisfied
end


%% Priors
% (Q,P,m,a)
log_P(3) = sum( (params.g_prior_phi-1).*log(g./params.g_prior_psi) ...
               - params.g_prior_phi.*g./params.g_prior_psi );

% (y's)
x = get_x(g,params.t,params.t_min,params.t_max);
log_P(4) = -sum( log(y), "all" ) + ( (params.N * params.K) / 2 ) * log(h/(2*pi)) - (h/2) * ( sum( log( y ./ x(:,2) ).^2 , "all" ) );

% (h)
log_P(5) = (params.h_prior_phi-1)*log(h/params.h_prior_psi) ...
         - h*params.h_prior_phi/params.h_prior_psi;


%% Posterior
log_P(1) = sum(log_P(2:end));


