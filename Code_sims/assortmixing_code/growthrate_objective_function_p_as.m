function f = growthrate_objective_function_p_as(x,params)

params.p_as = x;

curr_r_fixedpropasymp = get_r_SEIR_twodiseases_assortmixing(params);

fixed_r = params.fixed_r;
f = (curr_r_fixedpropasymp-fixed_r)^2;

end