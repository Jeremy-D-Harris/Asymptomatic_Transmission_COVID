
%% simulate SEIR model with fixed proportion of asymptomatic incidence, p
% same reproduction numbers, R0_a = R0_s

clear all; close all; clc;


%% want to save?
save_ans = 0;
% 0: don't save
% 1: save


%% which set of time scales?
% which_timescales = 3; % 1,2,3
% 1: same time scales: Ta=Ts=5 days
% 2: longer time scales of asymptomatic transmission: Ta=6,Ts=5 days
% 3: even longer time scales of asymptomatic transmission: Ta=8,Ts=5 days


filename = 'SEIR_fixedpropasymp_twodiseases_Ta8_varyrelR0_103022.mat';


%% set up colors and parameters
cbf_colors_db = [15,32,128]/255; % dark blue - same time scales
cbf_colors_v = [169,90,161]/255; % violet - longer time scale of asymptomatic
cbf_colors_lb = [133,192,249]/255; % light blue - even longer time scale of asymptomatc

cbf_colors_vector = [cbf_colors_db;cbf_colors_v;cbf_colors_lb];

cbf_colors = cbf_colors_vector(3,:);

%% vary infectious period, decay rates, days^-1
gamma_a=1/8; gamma_s=1/5;

% inital conditions to findint betas such that r=0.14
beta_a_init = 0.3275; beta_s_init = (beta_a_init/gamma_a)*gamma_s; 

% burnin time depends on parameters
t_end_burnin = 71.94;

% parameters
gamma_e=1/3; % 3 day exposure period

% params.beta_a = beta_a_init;
% params.beta_s =beta_s_init;
params.gamma_a = gamma_a;
params.gamma_s = gamma_s;
params.gamma_e = gamma_e;

% p is the proportion of asymptomatic incidence
proportion_asymp = 0.4;
params.p = proportion_asymp;
fixed_r = 0.14;
params.fixed_r = fixed_r;


%% next fit beta_a to little r
x0=beta_a_init;

fprintf('finding minimum wrt transmission rates... \n\n');

[x_soln,f_val] = fminsearch(@(x)objective_function(x,params),x0);

beta_a = x_soln(1);
beta_s = (beta_a/gamma_a)*gamma_s; 

params.beta_a = beta_a;
params.beta_s =beta_s;

fprintf('beta_a =  %2.4f \n\n',beta_a);
fprintf('beta_s =  %2.4f \n\n',beta_s);

% bestfit_rsquared = 1-(x_soln-fixed_r)^2;
% fprintf('best fit r^2 =  %2.4f \n\n',bestfit_rsquared);

bestfit_SSE = f_val;
fprintf('best fit SSE =  %1.2e \n\n',bestfit_SSE);

% need to get eigen proportion direction
eigen_direction_fixedpropasymp = get_eigendirection_SEIR_twodiseases_fixedpropasymp(params);

R0_fixedpropasymp = get_R0_SEIR_twodiseases_fixedpropasymp(params);

% calculate the intrinsic proportion of asymptomatic transmission
total_transmission_z = proportion_asymp*beta_a/gamma_a+(1-proportion_asymp)*beta_s/gamma_s;
asymp_transmission_z = proportion_asymp*beta_a/gamma_a;
proportion_asymp_transmission_z = asymp_transmission_z/total_transmission_z;

% calculate the realized proportion of asymptomatic transmission
total_transmission_q = beta_a*eigen_direction_fixedpropasymp(4)+beta_s*eigen_direction_fixedpropasymp(5);
asymp_transmission_q = beta_a*eigen_direction_fixedpropasymp(4);
proportion_asymp_transmission_q = asymp_transmission_q/total_transmission_q;

fprintf('Basic reproductive number \n');
fprintf('R_0 =  %2.4f \n\n',R0_fixedpropasymp);

r_fixedpropasymp = get_r_SEIR_twodiseases_fixedpropasymp(params);
fprintf('Exponential growth rate \n');
fprintf('r =  %2.4f \n\n',r_fixedpropasymp);


%% now simulate the two disease SEIR model
% asymptomatic and symptomatic infections

% mitigation parameters
params.t_m1 = 70;
params.t_min = 30;
params.t_m2 = params.t_m1+params.t_min;
params.mitigation_level = 1;

% 200 days is about 6-7 months
t_start = 0; t_end = t_end_burnin; % burn in time

dt=0.01;
params.dt=dt;
params.t_span = t_start:dt:t_end;


perturb = 1e-11;
if eigen_direction_fixedpropasymp(1)<0
    init_conds = [1;0;0;0;0;0;0] + perturb*eigen_direction_fixedpropasymp;
else
    init_conds = [1;0;0;0;0;0;0] - perturb*eigen_direction_fixedpropasymp;
end

options = odeset('RelTol',1e-10,'AbsTol',1e-12);

[t,y_traj_burnin] = ode45(@(t,y)simulate_SEIR_twodiseases_fixedpropasymp(t,y,params), params.t_span, init_conds,options);

t_start = 0; t_end = 250;

params.t_span = t_start:0.01:t_end;

init_conds = transpose(y_traj_burnin(end,:));

[t,y_traj] = ode45(@(t,y)simulate_SEIR_twodiseases_fixedpropasymp_mitigation(t,y,params), params.t_span, init_conds,options);

Rt_fixedpropasymp = get_Rt_SEIR_twodiseases_fixedpropasymp(params,y_traj);
results.Rt_fixedpropasymp=Rt_fixedpropasymp;

S_traj = y_traj(:,1);
E_a_traj = y_traj(:,2); E_s_traj = y_traj(:,3);
I_a_traj = y_traj(:,4); I_s_traj = y_traj(:,5);
R_a_traj = y_traj(:,6); R_s_traj = y_traj(:,7);

I_tot = I_a_traj + I_s_traj;
% results.I_tot=I_tot;

% calculate the total incidence
for count=1:length(params.t_span)
    this_t=params.t_span(count);
    [beta_a_traj(count,1), beta_s_traj(count,1)]= mitigation_function(this_t,params);
end

total_incidence = beta_a_traj.*(I_a_traj.*S_traj)+beta_s_traj.*(I_s_traj.*S_traj);
results.total_incidence=total_incidence;

% calculate generation interval distribution
g_a = @(t) gamma_e*gamma_a/(gamma_e-gamma_a)*(exp(-gamma_a*t)-exp(-gamma_e*t));
g_s = @(t) gamma_e*gamma_s/(gamma_e-gamma_s)*(exp(-gamma_s*t)-exp(-gamma_e*t));

GI_distribution_asymptomatic = feval(g_a,params.t_span);
GI_distribution_symptomatic = feval(g_s,params.t_span);


%%
f1 = figure(1); set(f1, 'Position', [400 250 450 850]);
subplot(4,1,1);
% q = semilogy(params.t_span, results.I_tot,'Color',cbf_colors,'LineWidth',2); hold on;
q = semilogy(params.t_span, total_incidence,'Color',cbf_colors,'LineWidth',2); hold on;
axis([0 params.t_span(end) 10^(-6) 1]);
xlabel('Time (days)'); ylabel({'Total'; 'incidence'});
% title(this_title);
% title(['p = ',num2str(proportion_asymp)])
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';

legend(q,{'All infections'},'FontWeight','normal','FontSize',11);
legend boxoff


% calculate the proportion of asymptomatic transmission
total_transmission = beta_a_traj.*(I_a_traj.*S_traj)+beta_s_traj.*(I_s_traj.*S_traj);
asymp_transmission = beta_a_traj.*(I_a_traj.*S_traj);
proportion_asymp_transmission = asymp_transmission./total_transmission;
results.proportion_asymp_transmission=proportion_asymp_transmission;

figure(1); subplot(4,1,2);
r(1) = plot(params.t_span, results.proportion_asymp_transmission,'Color',cbf_colors,'LineWidth',2); hold on;
% r(2) = plot(params.t_span, proportion_asymp*ones(size(params.t_span)),'k--','LineWidth',2); hold on;
axis([0 params.t_span(end) 0 1]);
xlabel('Time (days)'); ylabel({'Proportion'; 'asymptomatic'; 'transmission'});
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';

% calculate the proportion of asymptomatic incidence
asymp_incidence = proportion_asymp*(beta_a_traj.*(I_a_traj.*S_traj)+beta_s_traj.*(I_s_traj.*S_traj));
proportion_asymp_incidence = asymp_incidence./total_incidence;
results.proportion_asymp_incidence=proportion_asymp_incidence;

figure(1); subplot(4,1,3);
r(1) = plot(params.t_span, results.proportion_asymp_incidence,'Color',cbf_colors,'LineWidth',2); hold on;
% r(2) = plot(params.t_span, proportion_asymp*ones(size(params.t_span)),'k--','LineWidth',2); hold on;
axis([0 params.t_span(end) 0 1]);
xlabel('Time (days)'); ylabel({'Proportion'; 'asymptomatic'; 'incidence'});
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';



figure(1); subplot(4,1,4);
semilogy(params.t_span,results.Rt_fixedpropasymp,'Color',cbf_colors,'LineWidth',2);
axis([0 params.t_span(end) 10^-1 10^1]);
xlabel('Time (days)'); ylabel({'Effective'; 'reproduction'; 'number'});
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';


ind = find(total_incidence>10^-6);
fprintf('Time of first appearance: \n'); % want to be close to 25 days in
fprintf('%2.2f days \n\n',params.t_span(ind(1)));

init_proportion_asymp_incidence = proportion_asymp_incidence(1);
fprintf('Initial proportion asymptomatic incidence: \n'); % want to be close to 25 days in
fprintf('%2.2f \n\n',init_proportion_asymp_incidence);

% if you want also plot generation interval distribution
if 1
    
    f2 = figure(2); set(f2, 'Position', [1000   378   560   420]);
    r(1) = plot(params.t_span, GI_distribution_asymptomatic,'Color',cbf_colors,'LineWidth',2); hold on;
    r(2) = plot(params.t_span, GI_distribution_symptomatic,'Color',cbf_colors_db,'LineWidth',2); hold on;
    axis([0 21 0 0.2]);
    xlabel('Time (days)'); ylabel({'Probability'});
    title('Generation interval distributions');
    f2=gca;
    f2.LineWidth = 1;
    f2.FontSize = 14;
    f2.FontWeight = 'normal';
    
    legend(r,{['Asymptomatic GI distribution, T_a = ' num2str(1/params.gamma_a)],['Symptomatic GI distribution, T_s = ' num2str(1/params.gamma_s)]},'Location','NorthEast');
    legend box off
    
end


%%
% save simulated data
if save_ans==1
    folder_location = '../../Code_plt_ms_figures/sim_data/';
    save(strcat(folder_location,filename),'params','results');
    
    fprintf('Saved to file: \n'); % want to be close to 25 days in
    fprintf(strcat(filename,'\n'));
    
end
