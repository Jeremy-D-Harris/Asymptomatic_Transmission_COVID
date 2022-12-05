function f = growthrate_objective_function_propasymp(x,params)

params.p_aa = x;
params.p_as = params.p_aa;

curr_r_fixedpropasymp = get_r_SEIR_twodiseases_assortmixing(params);

fixed_r = params.fixed_r;
f = (curr_r_fixedpropasymp-fixed_r)^2;

end