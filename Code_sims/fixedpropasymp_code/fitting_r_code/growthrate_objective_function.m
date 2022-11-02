function f = growthrate_objective_function(x,params)

params.beta_a = x;
params.beta_s = params.k_relR0*(params.beta_a/params.gamma_a)*params.gamma_s;

curr_r_fixedpropasymp = get_r_SEIR_twodiseases_fixedpropasymp(params);

fixed_r = params.fixed_r;
f = (curr_r_fixedpropasymp-fixed_r)^2;

end