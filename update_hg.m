function [h,sample_g,rec] = update_hg(sample_g,rec, y, params)

% Multiplicative elliptical slice sampler (MESS)

M=length(params.T_sig);

sample_P = get_log_target(sample_g,y,params);

for rep = 1:params.T_rep

    % prepare sampler M=4
    sample_T = params.T_sig .* randn(M,1);
    tempor_T = params.T_sig .* randn(M,1);

    % pick slice
    U_prop = log(rand);
    
    % pick interval
    T = 2*pi*rand;
    T_min = T - 2*pi;
    T_max = T;
    
    while true 
        
        rec(2) = rec(2) + 1; % count proposals

        % generate proposal
        propos_T = cos(T)*sample_T + sin(T)*tempor_T;
        propos_t = sample_g.*exp(sample_T-propos_T);
        
        propos_P = get_log_target(propos_t,y,params);

        % first sum is jacobian terms
        log_a = sum(sample_T-propos_T)+(propos_P-sample_P);

        % carry acc test
        if U_prop < log_a
            sample_P = propos_P;
            sample_g = propos_t;

            rec(1) = rec(1) + 1; % count acceptances
            break % while true
        else
            % update intervals
            if T<0 
                T_min=T;
            
            else 
                T_max=T; 
            
            end
            % prepare next proposal
            T = T_min + (T_max-T_min)*rand;
        end % acc

    end % while true

end % rep

%%%  RECOVER H
x = get_x(sample_g,params.t,params.t_min,params.t_max);
h = gamrnd((params.N * params.K) /2 + params.h_prior_phi,get_B(y, x(:,2),params) );

end % ends update g function 

%% LOG TARGET FUNCTION
function out = get_log_target(g,y,params)
out = 0; 
% eval ODE
x = get_x(g,params.t,params.t_min,params.t_max); 
inner = sum(sum((log(y ./ x(:,2))).^2));

% y's 
out = out - ((params.N * params.K)/2 + params.h_prior_phi)...
* log(inner/2 + params.h_prior_phi / params.h_prior_psi );
%  (Q,P,m,a)
out = out + sum((params.g_prior_phi-1) .* log(g)) - sum(g .* params.g_prior_phi ./ params.g_prior_psi);

end % ends log target function

