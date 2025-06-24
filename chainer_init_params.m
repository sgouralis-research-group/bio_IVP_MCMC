function params = chainer_init_params(opts)


%% setup
params.units.time = opts.units_time; % units of time
params.units.conc = opts.units_conc; % units of concentration


%% data
params.t_min = opts.t_min; % [time]
params.t_max = opts.t_max; % [time]

params.t = opts.t; % [time]
params.z = opts.z; % [conc] stores either (mean) or (mean,std)

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


%% pre-processing
switch size(params.z,2)
    case 1
        % 1st raw moment
        params.Z = params.z;
    case 2
        % 1st and 2nd raw moments
        params.Z = [ params.z(:,1) , ...
                     params.z(:,1).^2 + (1-1/params.K)*params.z(:,2).^2 ];
    otherwise
        error('z must contain 1 or 2 statistics only')
end


%% priors

% dynamical parameters
params.g_prior_phi = nan(4,1);
params.g_prior_psi = nan(4,1);

% Q
params.g_prior_phi(1) = 1;
params.g_prior_psi(1) = params.z(end,1); % [conc]
% P
params.g_prior_phi(2) = 1;
params.g_prior_psi(2) = params.z(  1,1); % [conc]
% m
params.g_prior_phi(3) = 1;
params.g_prior_psi(3) = 10/(params.t_max-params.t_min); % 1/[time]
% a
params.g_prior_phi(4) = 1;
params.g_prior_psi(4) = params.g_prior_psi(3)/params.g_prior_psi(1); % 1/([conc]*[time])

% observational parameters
params.h_prior_phi = 1;
params.h_prior_psi = 10; % [1]


%% mESS settings for (Q,P,m,a)
params.T_sig = 0.05;
params.T_rep = 5;


%% cHMC settings for y
switch size(params.Z,2)
    case 1
        params.constraint_type = "linear";
        params.num_constraints = params.N;
        params.HMC_step_size   = 1e-9; % RATTLE step size
    case 2
        params.constraint_type = "quadratic";
        params.num_constraints = 2*params.N;
        params.HMC_step_size   = 2e-9; % RATTLE step size
end
params.HMC_L    = 20; % number of RATTLE steps
params.HMC_M    = 1e-6*kron(1./params.Z(:,1),ones(params.K,1)); % mass matrix
params.HMC_Minv = 1e+6*kron(   params.Z(:,1),ones(params.K,1)); % inverse mass matrix

params.fsolve_options = optimset('Display','off', ...
                                 'TolFun', 1e-10, ...
                                 'Algorithm', 'levenberg-marquardt');
params.initial_constraint_tolerance = 1e-6;


%% ground
if isfield(opts,'ground')
    params.ground = opts.ground;
end

