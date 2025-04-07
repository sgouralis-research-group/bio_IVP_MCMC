function B = get_B(y, rho, params)
  B = 1 / ( sum( log( y ./ rho ).^2 , "all" ) / 2 + params.h_prior_phi / params.h_prior_psi );
end