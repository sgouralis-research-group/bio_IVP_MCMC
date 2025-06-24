function c = get_c(y,params)

y = reshape(y,params.K,params.N)';

switch params.constraint_type

    case "linear"
        % sample mean only
        c =   sum(y,2)-params.K*params.Z(:,1) ;
    
    case "quadratic"
        % sample mean and second moment
        c = [ sum(y,2)-params.K*params.Z(:,1) ; ...
              0.5*params.Z(:,1).*(sum(y.^2./params.Z(:,2),2)-params.K) ];

end
