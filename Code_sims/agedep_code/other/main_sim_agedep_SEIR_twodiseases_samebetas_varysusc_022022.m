
%% simulate age-dependent SEIR model with asymptomatic and symptomatic transmission

clear all; close all; clc;


%% want to save?
save_ans = 0;
% 0: don't save
% 1: save


%% mitigation or not?
with_mitigation = 1;
% 0: no mitigation
% 1: with mitigation


%% which set of time scales?
which_timescales = 2; % 1,2
% 1: same time scales: Ta=Ts=5 days
% 2: longer time scales of asymptomatic transmission: Ta=8,Ts=5 days


%% variation in symptomaticity?
variation_symptomaticity_yesno =1; % 0,1
params.variation_symptomaticity_yesno=variation_symptomaticity_yesno;
% 0: p_n = p
% 1: p_n varies with age

% proportion_symptomatic_incidence = 0.3521;

%% load age-dependent contacts
load('agedep_data_Shanghai_update071521.mat');
% contact matrix: zhang et al. 2020, Table S8; NxN, N=12 --> compressed to N=8 in order to interface with age-distribution data
% susceptibility and symptomaticity estimates: from Davies et al. 2020, Fig 4 overall
% age distrubition of Shanghai: from Zhang et al. 2020, Table S2; originally derived from 2016 census


%% set up colors and parameters
cbf_colors_db = [15,32,128]/255; % dark blue - same time scales
cbf_colors_lb = [133,192,249]/255; % light blue - even longer time scale of asymptomatc

cbf_colors_vector = [cbf_colors_db;cbf_colors_lb];

if which_timescales==1
    
    cbf_colors = cbf_colors_vector(1,:);
    
    % decay rates, days^-1
    gamma_a=1/5; gamma_s=1/5;
    
    % same transmission rates
    this_beta = 0.0587;
    beta_a = this_beta; beta_s = this_beta;
    t_end_burnin = 163.29; % burn-in time
    
    
    % file names of simulations if saved
    if with_mitigation==0
        
        filename = 'SEIR_agedep_twodiseases_variationsusceptibility_T5and5.mat';
        
        fprintf('Susceptible Depletion ... \n\n');
        this_title = 'Susceptible Depletion';
        
        mitigation_level=1;
        
    else
        
        filename = 'SEIR_agedep_twodiseases_variationsusceptibility_T5and5_mit.mat';
        
        fprintf('With Mitigation ... \n\n');
        this_title = 'With Mitigation';
        
        mitigation_level=0.1465;
        
    end
    
elseif which_timescales==2
    
    cbf_colors = cbf_colors_vector(2,:);
    
    % decay rates, days^-1
    gamma_a=1/8; gamma_s=1/5;
    
    % same transmission rates
    this_beta = 0.0506;
    beta_a = this_beta; beta_s = this_beta;
    t_end_burnin = 164.83; % burn-in time
    
    
    if with_mitigation==0
        
        filename = 'SEIR_agedep_twodiseases_variationsusceptibility_T5and8.mat';
        
        fprintf('Susceptible Depletion ... \n\n');
        this_title = 'Susceptible Depletion';
        
        mitigation_level=1;
        
    else
        
        filename = 'SEIR_agedep_twodiseases_variationsusceptibility_T5and8_mit.mat';
        
        fprintf('With Mitigation ... \n\n');
        this_title = 'With Mitigation';
        
        mitigation_level=0.1108;
        
    end
    
end



%% deal with the contact matrix, age distributions, etc.
age_group_numbers = data_Shanghai_Davies.age_group_numbers;
N = length(age_group_numbers); % 8 age groups

% age distribution of Shanghai
age_distribution = data_Shanghai_Davies.age_distribution;
params.age_distribution = age_distribution;
N_total = sum(age_distribution); % total population of Shanghai

% Shanghai baseline contact matrix
Shanghai_contact_matrix = data_Shanghai_Davies.contact_matrix;

% % create a matrix of the diagonal elements
% diag_M = diag(diag(Shanghai_contact_matrix(1:N,1:N)));
%
% % create a matrix of the off-diagonal elements
% offdiag_M = Shanghai_contact_matrix(1:N,1:N)-diag_M;

% contact matrix
% M = diag_M+offdiag_M;
M = Shanghai_contact_matrix(1:N,1:N);
params.M = M;

% posterior means based off of Table 4 - Davies et al. nature letters 2020
susceptibility_byage= data_Shanghai_Davies.susc_infection;

% symptomatic probabiility by age
probability_symptomatic_byage= data_Shanghai_Davies.prob_symp;

params.sigma = susceptibility_byage;
params.prob_symp = probability_symptomatic_byage;
params.N = N;

% rescale contact matrix by N_j
% see supplemental info of Zhang et al. Science (2020)
for j=1:N
    
    N_j = age_distribution(j);
    M_tilde(:,j) = M(:,j)/N_j;
    
end
params.M_tilde = M_tilde;


%% parameters
gamma_e=1/3; % 3 day exposure period

params.beta_a = beta_a;
params.beta_s =beta_s;
params.gamma_a = gamma_a;
params.gamma_s = gamma_s;
params.gamma_e = gamma_e;

% mitigation parameters
params.t_m1 = 70;
params.t_min = 30;
params.t_m2 = params.t_m1+params.t_min;
params.mitigation_level = mitigation_level;

% 200 days is about 6-7 months
t_start = 0; t_end = t_end_burnin; % burn in time

dt=0.01;
params.dt=dt;
params.t_span = t_start:dt:t_end;

% get eigen proportion direction
eigen_direction_agedep_popN = get_eigendirection_SEIR_agedep(params);

R0_agedep = get_R0_SEIR_agedep(params);
fprintf('Basic reproductive number \n');
fprintf('R_0 =  %2.2f \n\n',R0_agedep);

r_agedep = get_r_SEIR_agedep(params);
results.r_agedep=r_agedep;
fprintf('Exponential growth rate \n'); % want to be close to 25 days in
fprintf('r =  %2.4f \n\n',r_agedep);

% pause;
perturb = 10^-11;
if eigen_direction_agedep_popN(1)>0
    init_conds = [age_distribution;zeros(6*N,1)] + perturb*(eigen_direction_agedep_popN)/norm(eigen_direction_agedep_popN)+perturb*ones(size(eigen_direction_agedep_popN));
else
    init_conds = [age_distribution;zeros(6*N,1)] - perturb*(eigen_direction_agedep_popN)/norm(eigen_direction_agedep_popN)+perturb*ones(size(eigen_direction_agedep_popN));
end

% pause;

options = odeset('RelTol',1e-12,'AbsTol',1e-12);

[t,y_burnin] = ode45(@(t,y)simulate_SEIR_agedep_twodiseases_popN(t,y,params), params.t_span, init_conds,options);

t_start = 0; t_end = 250;

params.t_span = t_start:dt:t_end;

init_conds = transpose(y_burnin(end,:));

[t,y_traj] = ode45(@(t,y)simulate_SEIR_agedep_twodiseases_popN_mitigation(t,y,params), params.t_span, init_conds,options);

S = y_traj(:,1:N);
E_a = y_traj(:,(N+1):(2*N)); E_s = y_traj(:,(2*N+1):(3*N));
I_a = y_traj(:,(3*N+1):(4*N)); I_s = y_traj(:,(4*N+1):(5*N));
R_a = y_traj(:,(5*N+1):(6*N)); R_s = y_traj(:,(6*N+1):(7*N));
results.y_traj=y_traj;

%% some calculations based on trajectories

I_all_byage = I_a + I_s;

I_a_total = sum(I_a,2); % total asymptomatic infections
I_s_total = sum(I_s,2); % total symptomatic infections
I_all_total = I_a_total+I_s_total; % total infections
I_all_total_fraction = I_all_total./N_total;
results.I_all_total_fraction=I_all_total_fraction;

% proportions of total infections
for count=1:length(params.t_span)
    I_a_proportion_infections(count,:) = I_a(count,:)/I_all_total(count);
    I_s_proportion_infections(count,:) = I_s(count,:)/I_all_total(count);
    I_all_proportion_infections(count,:) = I_all_byage(count,:)/I_all_total(count);
    
end

% get incidence: asymptomatic, symptomatic
[i_a_byage,i_s_byage] = get_incidence_SEIR_agedep(params,y_traj);

% total incidence by age
i_all_byage = i_a_byage + i_s_byage;

i_a_total = sum(i_a_byage,2); % total asymptomatic incidence across ages
i_s_total = sum(i_s_byage,2); % total symptomatic incidence across ages
i_all_total = sum(i_all_byage,2); % total incidence across age and symptom type

proportion_asymptomatic_incidence = i_a_total./i_all_total;
results.proportion_asymptomatic_incidence=proportion_asymptomatic_incidence;
fprintf('Initial proportion asymptomatic incidence: \n'); % want to be close to 25 days in
fprintf('p =  %2.4f \n\n',proportion_asymptomatic_incidence(1));

% total incidence as a fraction of the total population
i_all_total_fraction = i_all_total/N_total;
results.i_all_total_fraction=i_all_total_fraction;

% get proportion asymptomatic transmission
[transmission_a_byage,transmission_s_byage] = get_transmission_SEIR_agedep(params,y_traj);

% total incidence by age
transmission_all_byage = transmission_a_byage + transmission_s_byage;

transmission_a_total = sum(transmission_a_byage,2); % total asymptomatic transmission across ages
transmission_s_total = sum(transmission_s_byage,2); % total symptomatic transmission across ages
transmission_all_total = transmission_a_total+transmission_s_total; % total transmission across age and symptom type

proportion_asymptomatic_transmission = transmission_a_total./transmission_all_total;
results.proportion_asymptomatic_transmission=proportion_asymptomatic_transmission;
fprintf('Initial proportion asymptomatic transmission: \n'); % want to be close to 25 days in
fprintf('q =  %2.4f \n\n',proportion_asymptomatic_transmission(1));

% get Rt of the system
Rt_agedep = get_Rt_SEIR_agedep(params,y_traj);
results.Rt_agedep=Rt_agedep;

% calculate generation interval distribution
g_a = @(t) gamma_e*gamma_a/(gamma_e-gamma_a)*(exp(-gamma_a*t)-exp(-gamma_e*t));
g_s = @(t) gamma_e*gamma_s/(gamma_e-gamma_s)*(exp(-gamma_s*t)-exp(-gamma_e*t));

GI_distribution_asymptomatic = feval(g_a,params.t_span);
GI_distribution_symptomatic = feval(g_s,params.t_span);


%%
f1 = figure(1); set(f1, 'Position', [400 250 450 850]);
subplot(5,1,1);
if variation_symptomaticity_yesno==0
    semilogy(params.t_span, results.i_all_total_fraction,'Color',cbf_colors,'LineWidth',2); hold on;
else
    semilogy(params.t_span, results.i_all_total_fraction,'--','Color',cbf_colors,'LineWidth',2); hold on;
end
axis([0 params.t_span(end) 10^(-6) 1]);
xlabel('Time (days)'); ylabel({'Fraction'; 'infections'});
title(this_title);
% title(['p = ',num2str(proportion_asymp)])
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';



figure(1); subplot(5,1,2);
if variation_symptomaticity_yesno==0
    plot(params.t_span, results.proportion_asymptomatic_transmission,'Color',cbf_colors,'LineWidth',2); hold on;
else
    plot(params.t_span, results.proportion_asymptomatic_transmission,'--','Color',cbf_colors,'LineWidth',2); hold on;
end
axis([0 params.t_span(end) 0 1]);
xlabel('Time (days)'); ylabel({'Proportion'; 'asymptomatic'; 'transmission'});
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';


figure(1); subplot(5,1,3);
if variation_symptomaticity_yesno==0
    plot(params.t_span, results.proportion_asymptomatic_incidence,'Color',cbf_colors,'LineWidth',2); hold on;
else
    plot(params.t_span, results.proportion_asymptomatic_incidence,'--','Color',cbf_colors,'LineWidth',2); hold on;
end
axis([0 params.t_span(end) 0 1]);
xlabel('Time (days)'); ylabel({'Proportion'; 'asymptomatic'; 'incidence'});
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';




%% average age over time
age_vector = transpose(5:10:75);
for count=1:length(params.t_span)
    
    % find the average age of infecion
    ave_age_i_all(1,count) = (i_all_byage(count,1:N)*age_vector)/i_all_total(count);
    
end

results.ave_age_i_all=ave_age_i_all;

figure(1); subplot(5,1,4);
if variation_symptomaticity_yesno==0
    plot(params.t_span, results.ave_age_i_all,'Color',cbf_colors,'LineWidth',2); hold on;
else
    plot(params.t_span, results.ave_age_i_all,'--','Color',cbf_colors,'LineWidth',2); hold on;
end
axis([0 params.t_span(end) 10 50]);
xlabel('Time (days)'); ylabel({'Average';'age infection'});
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';



figure(1); subplot(5,1,5);
if variation_symptomaticity_yesno==0
    semilogy(params.t_span,results.Rt_agedep,'Color',cbf_colors,'LineWidth',2);
else
    semilogy(params.t_span,results.Rt_agedep,'--','Color',cbf_colors,'LineWidth',2);
end

axis([0 params.t_span(end) 10^-1 10^1]);
xlabel('Time (days)'); ylabel({'Effective'; 'reproduction'; 'number'});
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';


ind = find(i_all_total_fraction>10^-6);
fprintf('Time of first appearance: \n'); % want to be close to 25 days in
fprintf('%2.2f days \n\n',params.t_span(ind(1)));

% init_proportion_asymp_incidence = proportion_asymptomatic_incidence(1);
% fprintf('Initial proportion asymptomatic incidence: \n');
% fprintf('%2.4f \n\n',init_proportion_asymp_incidence);

% if you want also plot generation interval distribution
if 0
    
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
    
    
    
end

%% compute p(a|a), p(a|s) vs. time

Q_a = diag(ones(size(probability_symptomatic_byage))-probability_symptomatic_byage);
Q_s = diag(probability_symptomatic_byage);

A = diag(susceptibility_byage);


for count=1:length(params.t_span)
    
    
    this_I_a_vector = I_a(count,:);
    this_I_s_vector = I_s(count,:);
    this_S_vector = S(count,:);
    
    w_a(1,count) = this_S_vector*Q_a*A*M_tilde*this_I_a_vector'+this_S_vector*Q_s*A*M_tilde*this_I_a_vector';
    w_s(1,count) = this_S_vector*Q_a*A*M_tilde*this_I_s_vector'+this_S_vector*Q_s*A*M_tilde*this_I_s_vector';
    
    K_aa(1,count) = this_S_vector*Q_a*A*M_tilde*this_I_a_vector';
    K_as(1,count) = this_S_vector*Q_a*A*M_tilde*this_I_s_vector';
    
    % p(a|a) vs. time
    p_aa(1,count) = K_aa(1,count)/w_a(1,count);
    p_sa(1,count) = 1 - p_aa(1,count);
    
    % p(a|s) vs. time
    p_as(1,count) = K_as(1,count)/w_s(1,count);
    p_ss(1,count) = 1 - p_as(1,count);
    
    
end

results.p_aa=p_aa;
results.p_as=p_as;

if 0
    figure(2);
    this_p(1) = plot(params.t_span,results.p_aa,'LineWidth',2); hold on;
    this_p(2) = plot(params.t_span,results.p_as,'LineWidth',2); hold on;
    axis([0 params.t_span(end) 0 1]);
    xlabel('time (days)'); %ylabel('probability');
    f2=gca;
    f2.LineWidth = 1;
    f2.FontSize = 14;
    f2.FontWeight = 'normal';
    
    legend(this_p,{' p(a|a)',' p(a|s)'},'Fontsize',16,'Location','NorthWest')
    legend boxoff
end

Rt_agedep(end)

%%
% save simulated data
if save_ans==1
    
    folder_location = '../../Code_plt_ms_figures/sim_data/';
    save(strcat(folder_location,filename),'params','results','data_Shanghai_Davies');
    
    fprintf('Saved to file: \n');
    fprintf(strcat(filename,'\n'));
    
else
    
    fprintf('Not Saved. \n');
    
end

