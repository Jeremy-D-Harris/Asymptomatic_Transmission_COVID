function f = growthrate_objective_function_samebetas(x,params)

params.beta_a = x;
params.beta_s = params.beta_a;

curr_r_fixedpropasymp = get_r_SEIR_twodiseases_fixedpropasymp(params);

fixed_r = params.fixed_r;
f = (curr_r_fixedpropasymp-fixed_r)^2;

end