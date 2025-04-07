function grad_c = get_grad_c(y,params)

grad_c1 = kron(eye(params.N),repmat([1],params.K,1));

switch params.constraint_type

    case "linear"
        grad_c = grad_c1;

    case "quadratic"
        scale_y = reshape(y,params.K,params.N)'.* (params.z ./ params.v);
        scale_y = reshape(scale_y',params.N*params.K,1);
        rep_y = repmat(scale_y,1, params.N);
        % rep_y = repmat(y,1, params.N);
        grad_c2 = 0*rep_y;
        idx = logical(grad_c1);
        grad_c2(idx) = rep_y(idx);
        grad_c = [ grad_c1, grad_c2 ];
end

end