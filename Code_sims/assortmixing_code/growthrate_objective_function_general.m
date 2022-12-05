function f = growthrate_objective_function_general(x,params)

params.beta_a = x;
params.beta_s = params.relR0*(params.beta_a/params.gamma_a)*params.gamma_s;

curr_r_assortmixing = get_r_SEIR_twodiseases_assortmixing(params);

fixed_r = params.fixed_r;
f = (curr_r_assortmixing-fixed_r)^2;

end