function log_P = get_log_probs(g,h,y,params)
% log_P = log( post,like,priors )

if nargin<1
    log_P = 5; % posterior, likelihood, priors: Q,P,m,a,y,h
    return
end

% compute the log-posterior

log_P = nan(get_log_probs,1);

%% Likelihood
if params.HMC_RATTLE_constraint_error < params.HMC_RATTLE_constraint_tolerance
    log_P(2) = 0;  % constraint satisfied within tolerance
else
    log_P(2) = -Inf; % constraint not satisfied
end


%% Priors
% (Q,P,m,a)
log_P(3) = sum((params.g_prior_phi-1) .* log(g)) - sum(g .* params.g_prior_phi ./ params.g_prior_psi);

% (y's)
x = get_x(g,params.t,params.t_min,params.t_max);
rho = x(:,2);
log_P(4) = -sum( log(y), "all" ) + ( (params.N * params.K) / 2 ) * log(h/(2*pi)) - (h/2) * ( sum( log( y ./ rho ).^2 , "all" ) );
% (h)
log_P(5) = (params.h_prior_phi  - 1)* log(h) - ...
    h * params.h_prior_phi / params.h_prior_psi;


%% Posterior
log_P(1) = sum(log_P(2:end));


