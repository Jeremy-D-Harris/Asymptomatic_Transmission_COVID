function f = mitigation_objective_function(x,params)

params.mitigation_level = x;
% params.beta_s = params.k_relR0*(params.beta_a/params.gamma_a)*params.gamma_s;

% need to get eigen proportion direction
eigen_direction_fixedpropasymp = get_eigendirection_SEIR_twodiseases_fixedpropasymp(params);

% R0_fixedpropasymp = get_R0_SEIR_twodiseases_fixedpropasymp(params);
% fprintf('Basic reproductive number \n');
% fprintf('R_0 =  %2.4f \n\n',R0_fixedpropasymp);

% r_fixedpropasymp = get_r_SEIR_twodiseases_fixedpropasymp(params);
% fprintf('Exponential growth rate \n');
% fprintf('r =  %2.4f \n\n',r_fixedpropasymp);

perturb = 1e-11;
if eigen_direction_fixedpropasymp(1)<0
    init_conds = [1;0;0;0;0;0;0] + perturb*eigen_direction_fixedpropasymp;
else
    init_conds = [1;0;0;0;0;0;0] - perturb*eigen_direction_fixedpropasymp;
end

t_start = 0; t_end = 70;

params.t_span = t_start:0.01:t_end;

options = odeset('RelTol',1e-10,'AbsTol',1e-12);

[t,y_traj_burnin] = ode45(@(t,y)simulate_SEIR_twodiseases_fixedpropasymp(t,y,params), params.t_span, init_conds,options);

t_start = 0; t_end = 250;

params.t_span = t_start:0.01:t_end;

init_conds = transpose(y_traj_burnin(end,:));

[t,y_traj] = ode45(@(t,y)simulate_SEIR_twodiseases_fixedpropasymp(t,y,params), params.t_span, init_conds,options);

Rt_final_nomitigation_traj = get_Rt_SEIR_twodiseases_fixedpropasymp(params,y_traj);
Rt_final_nomitigation = Rt_final_nomitigation_traj(end);

[t,y_traj] = ode45(@(t,y)simulate_SEIR_twodiseases_fixedpropasymp_mitigation(t,y,params), params.t_span, init_conds,options);

this_Rt_fixedpropasymp_traj = get_Rt_SEIR_twodiseases_fixedpropasymp(params,y_traj);

curr_Rt_final = this_Rt_fixedpropasymp_traj(end);

f = (curr_Rt_final-Rt_final_nomitigation)^2;

end