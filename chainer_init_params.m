function params = chainer_init_params(opts)


%% setup
params.units.time = opts.units_time; % units of time
params.units.conc = opts.units_conc; % units of concentration


%% data
params.t_min = opts.t_min; % [time]
params.t_max = opts.t_max; % [time]

params.t = opts.t; % [time]
params.z = opts.z; % [conc]
params.v = opts.v; % [conc]^2

params.N = length(params.t);
params.K = opts.K;


%% sanity checks

if ~iscolumn(params.t)
    error('t must be column')
end

if any( diff([params.t_min;params.t;params.t_max])<=0 )
    error('t must be ordered')
end

if (size(params.z,1)~=params.N)
    error('z and t must have consistent dimensions')
end


%% priors

% dynamical parameters
params.g_prior_phi = nan(4,1);
params.g_prior_psi = nan(4,1);

% Q
params.g_prior_phi(1) = 2;
params.g_prior_psi(1) = 1.3e-1;   % [conc]
% P
params.g_prior_phi(2) = 2;
params.g_prior_psi(2) = 0.3e-3; % [conc]
% m
params.g_prior_phi(3) = 2;
params.g_prior_psi(3) = .5;    % 1/[time]
% a
params.g_prior_phi(4) = 2;
params.g_prior_psi(4) = 10; % 1/([conc]*[time])

% observational parameters
params.h_prior_phi = 250;
params.h_prior_psi = 25; % [1]

%% MESS sampling variance of (Q,P,m,a)
params.T_sig = 0.05*ones(4,1);
params.T_rep = 3; % number of internal repetitions

%% HMC parameters
params.constraint_type = opts.constraint_type; 
params.fsolve_options = optimset('Display','off', ...
                                 'TolFun', 1e-10, ...
                                 'Algorithm', 'levenberg-marquardt');
params.initial_constraint_tolerance = 1e-6;
switch params.constraint_type
    case "linear"
        params.num_constraints = params.N;
        params.HMC_step_size = 1e-9; % RATTLE step size
    case "quadratic"
        params.num_constraints = 2*params.N;
        params.HMC_step_size = 2e-9; % RATTLE step size
end
params.HMC_L = 20; % number of RATTLE steps
params.HMC_M    = 1e-6*kron(1./params.z,ones(params.K,1)); % mass matrix
params.HMC_Minv = 1e+6*kron(   params.z,ones(params.K,1)); % inverse mass matrix
params.HMC_RATTLE_constraint_tolerance = 1e-8;
params.HMC_RATTLE_constraint_error = 0;
params.HMC_RATTLE_tangency_error = 0;




%% ground
if isfield(opts,'ground')
    params.ground = opts.ground;
end

