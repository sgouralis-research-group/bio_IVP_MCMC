function [h,sample_g,rec] = update_hg(sample_g,rec,y,params)
% multiplicative elliptical slice sampler (mESS)

[sample_P,sample_B] = get_log_target(sample_g,y,params);

for rep = 1:params.T_rep

    % prepare sampler
    sample_T = params.T_sig*randn(4,1);
    tempor_T = params.T_sig*randn(4,1);

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
        propos_g = sample_g.*exp(sample_T-propos_T);

        [propos_P,propos_B] = get_log_target(propos_g,y,params);

        % first sum is jacobian terms
        log_a = sum(sample_T-propos_T) ...
                 + (propos_P-sample_P);

        % carry acc test
        if U_prop < log_a
            sample_P = propos_P;
            sample_B = propos_B;
            sample_g = propos_g;

            rec(1) = rec(1) + 1; % count acceptances
            break % while true
        else
            % update intervals
            if T<0; T_min=T; else; T_max=T; end
            % prepare next proposal
            T = T_min + (T_max-T_min)*rand;
        end % acc

    end % while true

end % rep


%%  recover h
h = sample_B * randg( 0.5*(params.N*params.K)+params.h_prior_phi );


end % function 



%% log target
function [P,B] = get_log_target(g,y,params)

x = get_x(g,params.t,params.t_min,params.t_max); 
B = get_B(y,x(:,2),params);

P = ( 0.5*(params.N*params.K)+params.h_prior_phi ) * log(B) ...
  + sum( (params.g_prior_phi-1).*log(g./params.g_prior_psi) - params.g_prior_phi.*g./params.g_prior_psi );

end

