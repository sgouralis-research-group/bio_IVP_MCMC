function [h,y,rec] = update_hy(y, rec, g, params)

x = get_x(g,params.t,params.t_min,params.t_max);
dyn_params = params;
dyn_params.rho = x(:,2);

%% sample y with constrained HMC

% get initial position and momentum
sample_position = reshape(y',params.N*params.K,1);
grad_c = get_grad_c(sample_position,params);
sample_momentum = fsolve( @(momentum) grad_c'*momentum, ...
                          randn(params.N*params.K,1), ...
                          params.fsolve_options ...
                        );

% integrate Hamiltonian
% sample random step size (exponential distribution) FIXME
[propos_position, propos_momentum] = integrate_Hamiltonian(sample_position,sample_momentum,dyn_params);

% % check reversibility
% [revers_position, revers_momentum] = integrate_Hamiltonian(propos_position,-propos_momentum,dyn_params);
% revpos_relerr = norm(sample_position-revers_position)/norm(sample_position);
% revmom_relerr = norm(sample_momentum+revers_momentum)/norm(sample_momentum);
% disp(revpos_relerr);
% disp(revmom_relerr);

% accept or reject
sample_energy = get_Hamiltonian_energy(sample_position,sample_momentum,dyn_params);
propos_energy = get_Hamiltonian_energy(propos_position,propos_momentum,dyn_params);
logR = sample_energy - propos_energy;

u = rand;

if log(u) <= logR && isequal(propos_position,abs(propos_position))
    
    % % check constraint satisfaction
    % params.HMC_RATTLE_constraint_error = max(abs(get_c(propos_position,params)));
    % params.HMC_RATTLE_tangency_error = max((normc(get_grad_c(propos_position,params)))' * normc(params.HMC_Minv * propos_momentum),[],'all');
    % disp(params.HMC_RATTLE_constraint_error)
    % disp(params.HMC_RATTLE_tangency_error);

    y = reshape(propos_position,params.K,params.N)'; % accept
    rec(1) = rec(1) + 1;
end
rec(2) = rec(2) + 1;

%% sample h with direct sampling

A =  ( params.N * params.K ) /2 + params.h_prior_phi;
B = get_B(y,dyn_params.rho,params);
h = gamrnd(A,B);
% h = params.ground.h;


end

%%

function [q,p] = integrate_Hamiltonian(q,p,params)
% Inputs
% q = position column vector
% p = momentum column vector
% params = struct of parameters

lambda_q = zeros(params.num_constraints,1);
p = params.HMC_Minv * p;
for l=1:params.HMC_L
    [q,p,lambda_q] = one_step_RATTLE(q,p,lambda_q,params);
end
p = params.HMC_M * p;


end

%%

function H = get_Hamiltonian_energy(y,p,params)
% Inputs
% q = position
% p = momentum
% params = struct of fixed parameters

potential_energy = get_potential_energy(y, params); % potential energy
kinetic_energy = p' * p / 2; % kinetic energy
H = potential_energy + kinetic_energy; % total energy of unconstrained system
end

%%

function V = get_potential_energy(y, params)
% This needs to be changed or a switch added depending on the problem

% recover params
y_mat = reshape(y,params.K,params.N)';
A = ( params.N * params.K ) / 2 + params.h_prior_phi;
B = get_B(y_mat,params.rho,params);

V = A * log(1/B) + sum( log(y) );
end

%%

function [q,p,lambda_q] = one_step_RATTLE(q,p,lambda_q,params)
% one step RATTLE
%
% q = position
% p = momentum
% lambda_q = previous iterate's position constraint multiplier

h = params.HMC_step_size; % integration step size

% (half momentum update)
grad_V = get_grad_V(q,params);
grad_c = get_grad_c(q,params);

switch params.constraint_type

    case "linear"
        % get lambda_q for just sample mean
        grad_c_dag = linsolve(grad_c'*grad_c,grad_c');
        f = params.K*grad_c_dag'*params.z;
        lambda_q = grad_c_dag * ( (2/h^2)*(q-f) + (2/h)*p - grad_V );
    
    case "quadratic"
        % get lambda_q with mean and variance constraints with one step Newton-Raphson
        f = q + h*p - (h^2/2)*(grad_V + grad_c*lambda_q);
        J = -(h^2/2)*get_grad_c(f,params)'*grad_c;
        correction = linsolve(J,get_c(f,params));
        lambda_q = lambda_q - correction;

end

% (half momentum update)
p = p - (h/2) * (grad_V + grad_c * lambda_q);

% (full position update)
q = q + h * p;
% q = abs(q + h * p);

% (half momentum update)
grad_V = get_grad_V(q,params);
grad_c = get_grad_c(q,params);
lambda_p = linsolve(grad_c' * grad_c, grad_c' * (p - (h/2) * grad_V));
% lambda_p_residual = grad_c' * grad_c * lambda_p - grad_c' * (p - (h/2) * grad_V);
% disp(lambda_p_residual);
p = p - (h/2) * (grad_V + grad_c * lambda_p);

end

%%

function grad_V = get_grad_V(y,dyn_params)
  % This needs to be changed or a switch added depending on the problem
  
  y_mat = reshape(y,dyn_params.K,dyn_params.N)';
  A = ( dyn_params.N*dyn_params.K ) / 2 + dyn_params.h_prior_phi;
  B = get_B(y_mat, dyn_params.rho, dyn_params);
  rep_rho = kron(dyn_params.rho,ones(dyn_params.K,1));
  grad_V = (A*B*rep_rho./y) .* log(y./rep_rho) + 1./y;

end
    
    