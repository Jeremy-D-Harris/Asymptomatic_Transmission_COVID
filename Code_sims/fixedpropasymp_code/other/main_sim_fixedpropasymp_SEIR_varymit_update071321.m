
%% simulate SEIR model with fixed proportion of asymptomatic incidence, p
% same reproduction numbers, R0_a = R0_s
% vary mitigation level

clear all; close all; clc;


%% want to save?
save_ans = 0;
% 0: don't save
% 1: save

%% which set of time scales?
which_mitigation_level = 3; % 1,2,3
% 1: half mitigation, 1/5 contacts
% 2: same as Figure 1, 1/10 contacts
% 3: twice mitigation, 1/20 contacts


%% set up colors and parameters
cbf_colors_vector_blk = [0, 0, 0]; % black
cbf_colors_vector_lb = [133,192,249]/255; % light blue
cbf_colors_vector_gray = [0.5, 0.5, 0.5]; % gray

cbf_colors_vector = [cbf_colors_vector_blk;cbf_colors_vector_lb;cbf_colors_vector_gray];


% decay rates, days^-1
gamma_a=1/8; gamma_s=1/5;

% set betas s.t. R0,a=R0,s are the same and r=0.14
beta_a = 0.3287; beta_s = 0.5260;

% burnin time depends on parameters
t_end_burnin = 71.42;

if which_mitigation_level==1
    
    filename = 'SEIR_fixedpropasymp_twodiseases_sameR0s_071321_T5and8_mit1.mat';
    
    mitigation_level=1/5; % 1/5 baseline contact rates
    t_min = 10;
    fprintf('With mitigation level, 1/5 baseline contacts \n\n');
    
elseif which_mitigation_level==2
    filename = 'SEIR_fixedpropasymp_twodiseases_sameR0s_071321_T5and8_mit2.mat';
    
    mitigation_level=1/10; % 1/5 baseline contact rates
    t_min = 10;
    fprintf('With mitigation level, 1/10 baseline contacts \n\n');
    
else
    filename = 'SEIR_fixedpropasymp_twodiseases_sameR0s_071321_T5and8_mit3.mat';
    
    mitigation_level=1/20; 
    t_min = 10;
    fprintf('With mitigation level, 1/20 baseline contacts \n\n');
end


%% simulate the two disease SEIR model with two infectious compartments:
% asymptomatic and symptomatic infections

% parameters
gamma_e=1/3; % 3 day exposure period

params.beta_a = beta_a;
params.beta_s =beta_s;
params.gamma_a = gamma_a;
params.gamma_s = gamma_s;
params.gamma_e = gamma_e;

% mitigation parameters
params.t_m1 = 70;
params.t_min = t_min;
params.t_m2 = params.t_m1+params.t_min;
params.mitigation_level = mitigation_level;

% 200 days is about 6-7 months
t_start = 0; t_end = t_end_burnin; % burn in time

dt=0.01;
params.dt=dt;
params.t_span = t_start:dt:t_end;

% p is the proportion of asymptomatic incidence
proportion_asymp = 0.4;
params.p = proportion_asymp;

% need to get eigen proportion direction
eigen_direction_fixedpropasymp = get_eigendirection_SEIR_twodiseases_fixedpropasymp(params);

R0_fixedpropasymp = get_R0_SEIR_twodiseases_fixedpropasymp(params);
fprintf('Basic reproductive number \n');
fprintf('R_0 =  %2.2f \n\n',R0_fixedpropasymp);

r_fixedpropasymp = get_r_SEIR_twodiseases_fixedpropasymp(params);
fprintf('Exponential growth rate \n');
fprintf('r =  %2.2f \n\n',r_fixedpropasymp);

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
results.I_tot=I_tot;

% calculate the total incidence
for count=1:length(params.t_span)
    this_t=params.t_span(count);
    [beta_a_traj(count,1), beta_s_traj(count,1)]= mitigation_function(this_t,params);
end

total_incidence = beta_a_traj.*(I_a_traj.*S_traj)+beta_s_traj.*(I_s_traj.*S_traj);

% calculate generation interval distribution
g_a = @(t) gamma_e*gamma_a/(gamma_e-gamma_a)*(exp(-gamma_a*t)-exp(-gamma_e*t));
g_s = @(t) gamma_e*gamma_s/(gamma_e-gamma_s)*(exp(-gamma_s*t)-exp(-gamma_e*t));

GI_distribution_asymptomatic = feval(g_a,params.t_span);
GI_distribution_symptomatic = feval(g_s,params.t_span);


%%
f1 = figure(1); set(f1, 'Position', [400 250 450 850]);
subplot(4,1,1);
semilogy(params.t_span, results.I_tot,'Color',cbf_colors_vector(which_mitigation_level,:),'LineWidth',2); hold on;
axis([0 params.t_span(end) 10^(-6) 1]);
xlabel('Time (days)'); ylabel({'Fraction'; 'infections'});
% title(['p = ',num2str(proportion_asymp)])
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';

% calculate the proportion of asymptomatic incidence
asymp_incidence = proportion_asymp*(beta_a_traj.*(I_a_traj.*S_traj)+beta_s_traj.*(I_s_traj.*S_traj));
proportion_asymp_incidence = asymp_incidence./total_incidence;
results.proportion_asymp_incidence=proportion_asymp_incidence;

figure(1); subplot(4,1,2);
r(1) = plot(params.t_span, results.proportion_asymp_incidence,'Color',cbf_colors_vector(which_mitigation_level,:),'LineWidth',2); hold on;
% r(2) = plot(params.t_span, proportion_asymp*ones(size(params.t_span)),'k--','LineWidth',2); hold on;
axis([0 params.t_span(end) 0 1]);
xlabel('Time (days)'); ylabel({'Proportion'; 'asymptomatic'; 'incidence'});
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';

% calculate the proportion of asymptomatic transmission
total_transmission = beta_a_traj.*(I_a_traj.*S_traj)+beta_s_traj.*(I_s_traj.*S_traj);
asymp_transmission = beta_a_traj.*(I_a_traj.*S_traj);
proportion_asymp_transmission = asymp_transmission./total_transmission;
results.proportion_asymp_transmission=proportion_asymp_transmission;

figure(1); subplot(4,1,3);
r(1) = plot(params.t_span, results.proportion_asymp_transmission,'Color',cbf_colors_vector(which_mitigation_level,:),'LineWidth',2); hold on;
axis([0 params.t_span(end) 0 1]);
xlabel('Time (days)'); ylabel({'Proportion'; 'asymptomatic'; 'transmission'});
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';



figure(1); subplot(4,1,4);
semilogy(params.t_span,results.Rt_fixedpropasymp,'Color',cbf_colors_vector(which_mitigation_level,:),'LineWidth',2); hold on;
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


%% save simulated data
if save_ans==1
    folder_location = '../../Code_plt_ms_figures/sim_data/';
    save(strcat(folder_location,filename));
    
    fprintf('Saved to file: \n'); % want to be close to 25 days in
    fprintf(strcat(filename,'\n'));
    
end
