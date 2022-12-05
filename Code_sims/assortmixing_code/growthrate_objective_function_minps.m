function f = growthrate_objective_function_minps(x,params)

params.p_aa = x(1);
params.p_as = x(2);

curr_r_fixedpropasymp = get_r_SEIR_twodiseases_assortmixing(params);

% need to get eigen proportion direction
eigen_direction_assortmixing = get_eigendirection_SEIR_twodiseases_assortmixing(params);

init_total_incidence = params.beta_a*eigen_direction_assortmixing(4)+params.beta_s*eigen_direction_assortmixing(5);
init_asymp_incidence = params.beta_a*eigen_direction_assortmixing(4);
curr_init_prop_asymp = init_asymp_incidence/init_total_incidence;


fixed_r = params.fixed_r;
fixed_p = params.fixed_p;


fixed_p_aa = params.fixed_p_aa; fixed_p_as = params.fixed_p_as;
% fixed_ps = [fixed_p_aa;fixed_p_as]; 
f = (curr_r_fixedpropasymp-fixed_r)^2*(1+(curr_init_prop_asymp-fixed_p)^2+(curr_r_fixedpropasymp-fixed_r)^2*((x(1)-fixed_p_aa).^2+(x(2)-fixed_p_as).^2));
% f = (curr_r_fixedpropasymp-fixed_r)^2+(curr_init_prop_asymp-fixed_p)^2;
% f = (curr_r_fixedpropasymp-fixed_r)^2;

end