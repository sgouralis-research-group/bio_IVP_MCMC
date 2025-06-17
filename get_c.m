function c = get_c(y,params)
    % units: [c(y)] = [y]
    y = reshape(y,params.K,params.N)';
    switch params.constraint_type
        case "linear"
            % sample mean only
            c = sum(y,2)-params.K*params.z;
        case "quadratic"
            % sample mean and second moment
            c = [ sum(y,2)-params.K*params.z     ; ...
                params.z .*(sum((y.^2)./params.v,2)-params.K)/2 ];
    end
end