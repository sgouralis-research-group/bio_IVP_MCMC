function [h,y,rec] = update_hy(y, rec, g, params)
% constraint Hamiltonian Monte Carlo (cHMC)

% prepare
x = get_x(g,params.t,params.t_min,params.t_max);

% get initial position and momentum
sample_position = reshape(y',params.N*params.K,1);
grad_c = get_grad_c(sample_position,params);
sample_momentum = fsolve( @(momentum) grad_c'*momentum, ...
                          randn(params.N*params.K,1), ...
                          params.fsolve_options );

% integrate Hamiltonian
[propos_position, propos_momentum] = integrate_Hamiltonian(sample_position,sample_momentum,x(:,2),params);

% carry acc test
[sample_H,sample_B] = get_Hamiltonian_energy(sample_position,sample_momentum,x(:,2),params);
[propos_H,propos_B] = get_Hamiltonian_energy(propos_position,propos_momentum,x(:,2),params);
if all(propos_position>0,'all') && rand < exp(sample_H-propos_H)
    y = reshape(propos_position,params.K,params.N)';
    sample_B = propos_B;
    rec(1) = rec(1) + 1; % count acceptances
end
    rec(2) = rec(2) + 1; % count proposals


%% sample h with direct sampling
h = sample_B * randg( 0.5*(params.N*params.K)+params.h_prior_phi );


end % function




%%
function [q,p] = integrate_Hamiltonian(q,p,rho,params)

lambda_q = zeros(params.num_constraints,1);

p = params.HMC_Minv .* p;
for l=1:params.HMC_L
    [q,p,lambda_q] = one_step_RATTLE(q,p,lambda_q,rho,params);
end
p = params.HMC_M .* p;

end


%%
function [H,B] = get_Hamiltonian_energy(y,p,rho,params)
[V,B] = get_potential_energy(y,rho,params);
H = V + 0.5*(p'*p);
end


%%
function [V,B] = get_potential_energy(y,rho,params)
y = reshape(y,params.K,params.N)';
B = get_B(y,rho,params);
V = sum(log(y),'all') ...
  - ( 0.5*(params.N*params.K)+params.h_prior_phi )*log(B);
end


%%
function [q,p,lambda_q] = one_step_RATTLE(q,p,lambda_q,rho,params)
h = params.HMC_step_size; % integration step size

% (half momentum update)
grad_V = get_grad_V(q,rho,params);
grad_c = get_grad_c(q,    params);

switch params.constraint_type

    case "linear"
        % get lambda_q for just sample mean
        grad_c_dag = linsolve(grad_c'*grad_c,grad_c');
        f = params.K*grad_c_dag'*params.Z(:,1);
        lambda_q = grad_c_dag * ( (2/h^2)*(q-f) + (2/h)*p - grad_V );
    
    case "quadratic"
        % get lambda_q with mean and variance constraints with one step Newton-Raphson
        f = q + h*p - (h^2/2)*(grad_V + grad_c*lambda_q);
        J = -(h^2/2)*get_grad_c(f,params)'*grad_c;
        lambda_q = lambda_q - linsolve(J,get_c(f,params));

end

% (half momentum update)
p = p - (h/2) * (grad_V + grad_c * lambda_q);

% (full position update)
q = q + h*p;

% (half momentum update)
grad_V = get_grad_V(q,rho,params);
grad_c = get_grad_c(q,    params);
lambda_p = linsolve(grad_c' * grad_c, grad_c' * (p - (h/2) * grad_V));
p = p - (h/2) * (grad_V + grad_c * lambda_p);

end


%%
function grad_V = get_grad_V(y,rho,params)
A = 0.5*(params.N*params.K) + params.h_prior_phi;
B = get_B( reshape(y,params.K,params.N)' ,rho,params);
rep_rho = kron(rho,ones(params.K,1));
grad_V = A*B*rep_rho./y .* log(y./rep_rho) + 1./y;
end