function B = get_B(y, rho, params)
  B = 1 / ( 0.5 * sum( log( y./rho ).^2 , "all" ) ...
          + params.h_prior_phi/params.h_prior_psi );
end