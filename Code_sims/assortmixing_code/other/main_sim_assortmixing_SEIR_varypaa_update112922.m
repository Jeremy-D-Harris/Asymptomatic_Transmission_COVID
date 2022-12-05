
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

filename = 'SEIR_assortmixing_twodiseases_varyTa_112922.mat';


%% set up colors and parameters
% cbf_colors_db = [15,32,128]/255; % dark blue - same time scales
% cbf_colors_v = [169,90,161]/255; % violet - longer time scale of asymptomatic
% cbf_colors_lb = [133,192,249]/255; % light blue - even longer time scale of asymptomatc
%
% cbf_colors_vector = [cbf_colors_db;cbf_colors_v;cbf_colors_lb];

% cbf_colors = cbf_colors_vector(2,:);
cbf_colors = [0.5,0.5,0.5];


%% burnin time depends on parameters
t_end_burnin = 72.02;

% parameters
gamma_e=1/3; % 3 day exposure period

% params.beta_a = beta_a_init;
% params.beta_s =beta_s_init;
params.gamma_e = gamma_e;

% p is the proportion of asymptomatic incidence
proportion_asymp = 0.4;
params.p = proportion_asymp;
% along the line = b-m*x
% proportion_asymp = 0.7;
% p_as = 0.2;



% fixed exponential growth rate
fixed_r = 0.14;
params.fixed_r = fixed_r;

%% vary relative infectious periods, decay rates, days^-1
k_vector_relT = linspace(1,8/5,10); % ratio of Ta/Ts
params.k_vector_relT = k_vector_relT;
vary_p_aa = linspace(0.4,1,10);

% Ts_vector=5*ones(size(k_vector_relT));
% Ta_vector=k_vector_relT.*Ts_vector;

relR0 = 1;
params.relR0 = relR0;

% give an initial guess
beta_a_init = 0.4149;

% fixed infectious periods: Ts=5, Ta=8
gamma_a = 1/5; gamma_s = 1/8;
params.gamma_a = gamma_a;
params.gamma_s = gamma_s;

%% next fit beta_a to little r
for count=1:length(vary_p_aa)
    
    
    
    
    if count == 1
        
        p_aa = vary_p_aa(count); p_as = p_aa;
        p_as_init = p_aa;
%         p_aa = proportion_asymp; p_as = proportion_asymp;
        params.p_aa = p_aa; params.p_as = p_as;
        
        % inital conditions to findint betas such that r=0.14
        
        x0=beta_a_init;
        beta_s_init = relR0*(beta_a_init/gamma_a)*gamma_s;
        
        fprintf('finding minimum wrt transmission rates... \n\n');
        fprintf('count = ',count);
        
        [x_soln,f_val] = fminsearch(@(x)growthrate_objective_function_general(x,params),x0);
        
        beta_a_fixedp = x_soln(1);
        beta_s_fixedp = relR0*(beta_a_fixedp/gamma_a)*gamma_s;
        
        beta_a = beta_a_fixedp;
        beta_s = beta_s_fixedp;
        params.beta_a=beta_a;
        params.beta_s=beta_s;
        
        
    else
        
        p_aa = vary_p_aa(count);
        x0=p_as_init;
        
        params.p_aa = p_aa;
        params.p_as = p_as;
        
        
        fprintf('finding minimum wrt conditionals, p_aa, p_as... \n\n');
        fprintf('count = ',count);
        
        [x_soln,f_val] = fminsearch(@(x)growthrate_objective_function_conditionals(x,params),x0);
        
        p_as = x_soln(1);
        params.p_as = p_as;
    
        % initialize for next point
        p_as_init = p_as;   

        
    end
    
%     % initialize for next point
%     p_as_init = p_as;
    
    params.p_aa = p_aa;
    params.p_as = p_as;
    
    fprintf('p_aa =  %2.4f \n\n',p_aa);
    fprintf('p_as =  %2.4f \n\n',p_as);
    
    % bestfit_rsquared = 1-(x_soln-fixed_r)^2;
    % fprintf('best fit r^2 =  %2.4f \n\n',bestfit_rsquared);
    
    bestfit_SSE = f_val;
    fprintf('best fit SSE =  %1.2e \n\n',bestfit_SSE);
    
    results.bestfit_SSE=bestfit_SSE;
    
    % need to get eigen proportion direction
    eigen_direction_assortmixing = get_eigendirection_SEIR_twodiseases_assortmixing(params);
    
    R0_assortmixing = get_R0_SEIR_twodiseases_assortmixing(params);
    
    % calculate the intrinsic proportion of asymptomatic transmission
    total_transmission_z = proportion_asymp*beta_a/gamma_a+(1-proportion_asymp)*beta_s/gamma_s;
    asymp_transmission_z = proportion_asymp*beta_a/gamma_a;
    proportion_asymp_transmission_z = asymp_transmission_z/total_transmission_z;
    results.proportion_asymp_transmission_z=proportion_asymp_transmission_z;
    
    % calculate the realized proportion of asymptomatic transmission
    total_transmission_q = beta_a*eigen_direction_assortmixing(4)+beta_s*eigen_direction_assortmixing(5);
    asymp_transmission_q = beta_a*eigen_direction_assortmixing(4);
    proportion_asymp_transmission_q = asymp_transmission_q/total_transmission_q;
    results.proportion_asymp_transmission_q=proportion_asymp_transmission_q;
    
    fprintf('Basic reproductive number \n');
    fprintf('R_0 =  %2.4f \n\n',R0_assortmixing);
    
    r_assortmixing = get_r_SEIR_twodiseases_assortmixing(params);
    fprintf('Exponential growth rate \n');
    fprintf('r =  %2.4f \n\n',r_assortmixing);
    
    
    %% now simulate the two disease SEIR model
    % asymptomatic and symptomatic infections
    
    % mitigation parameters
    params.t_m1 = 70;
    params.t_min = 30;
    params.t_m2 = params.t_m1+params.t_min;
    params.mitigation_level = 1; % 1 = unmitigated
    
    % 200 days is about 6-7 months
    t_start = 0; t_end = t_end_burnin; % burn in time
    
    dt=0.01;
    params.dt=dt;
    params.t_span = t_start:dt:t_end;
    
    
    perturb = 1e-11;
    if eigen_direction_assortmixing(1)<0
        init_conds = [1;0;0;0;0;0;0] + perturb*eigen_direction_assortmixing;
    else
        init_conds = [1;0;0;0;0;0;0] - perturb*eigen_direction_assortmixing;
    end
    
    options = odeset('RelTol',1e-10,'AbsTol',1e-12);
    
    [t,y_traj_burnin] = ode45(@(t,y)simulate_SEIR_twodiseases_assortmixing(t,y,params), params.t_span, init_conds,options);
    
    t_start = 0; t_end = 250;
    
    params.t_span = t_start:0.01:t_end;
    
    init_conds = transpose(y_traj_burnin(end,:));
    
    [t,y_traj] = ode45(@(t,y)simulate_SEIR_twodiseases_assortmixing_mitigation(t,y,params), params.t_span, init_conds,options);
    
    this_Rt_assortmixing = get_Rt_SEIR_twodiseases_assortmixing(params,y_traj);
    results.Rt_assortmixing=this_Rt_assortmixing;
    
    S_traj = y_traj(:,1);
    E_a_traj = y_traj(:,2); E_s_traj = y_traj(:,3);
    I_a_traj = y_traj(:,4); I_s_traj = y_traj(:,5);
    R_a_traj = y_traj(:,6); R_s_traj = y_traj(:,7);
    
    I_tot = I_a_traj + I_s_traj;
    % results.I_tot=I_tot;
    
    % calculate the total incidence
    for counter=1:length(params.t_span)
        this_t=params.t_span(counter);
        [beta_a_traj(counter,1), beta_s_traj(counter,1)]= mitigation_function(this_t,params);
    end
    
    
    % calculate the proportion of asymptomatic incidence
    this_total_incidence = beta_a_traj.*(I_a_traj.*S_traj)+beta_s_traj.*(I_s_traj.*S_traj);
    this_asymp_incidence = p_aa*(beta_a_traj.*(I_a_traj.*S_traj)+p_as*beta_s_traj.*(I_s_traj.*S_traj));
    this_proportion_asymp_incidence = this_asymp_incidence./this_total_incidence;
    
    results.total_incidence=this_total_incidence;
    results.proportion_asymp_incidence=this_proportion_asymp_incidence;
    
    % calculate the proportion of asymptomatic transmission
    this_total_transmission = beta_a_traj.*(I_a_traj.*S_traj)+beta_s_traj.*(I_s_traj.*S_traj);
    this_asymp_transmission = beta_a_traj.*(I_a_traj.*S_traj);
    this_proportion_asymp_transmission = this_asymp_transmission./this_total_transmission;
    
    results.total_transmission = this_total_transmission;
    results.proportion_asymp_transmission=this_proportion_asymp_transmission;
    
    
    %% print a few important points
    ind = find(this_total_incidence>10^-6);
    fprintf('Time of first appearance: \n'); % want to be close to 25 days in
    fprintf('%2.2f days \n\n',params.t_span(ind(1)));
    
    init_proportion_asymp_incidence = this_proportion_asymp_incidence(1);
    fprintf('Initial proportion asymptomatic incidence: \n'); % want to be close to 25 days in
    fprintf('%2.2f \n\n',init_proportion_asymp_incidence);
    
    %% collect parameters and results
    params_collect(count) = params;
    results_collect(count) = results;
    
    %% plot figure
    f1 = figure(1); set(f1, 'Position', [400 250 450 850]);
    
    subplot(4,1,1);
    % q = semilogy(params.t_span, results.I_tot,'Color',cbf_colors,'LineWidth',2); hold on;
    %     this_h(count) = semilogy(params.t_span, this_total_incidence,'Color',cbf_colors,'LineWidth',2); hold on;
    %     this_h(count).Color(4) = 1-0.12*(count);
    this_h(count) = semilogy(params.t_span, this_total_incidence,'Color',cbf_colors,'LineWidth',2); hold on;
    %     this_transparency =
    %     increment = 1/(1+length(relR0));
    %     this_h(count).Color(4) = 1-increment*count;
    
    axis([0 params.t_span(end) 10^(-6) 1]);
    xlabel('Time (days)'); ylabel({'Total'; 'incidence'});
    % title(this_title);
    % title(['p = ',num2str(proportion_asymp)])
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    set(f1, 'YScale', 'log');
    
    %     if count==4
    %         legend(this_h,{'1','2','3','4'},'FontWeight','normal','FontSize',11);
    %         legend boxoff
    %     end
    
    
    figure(1); subplot(4,1,2);
    this_p = plot(params.t_span, this_proportion_asymp_transmission,'Color',cbf_colors,'LineWidth',2); hold on;
    %     this_p.Color(4) = 1-increment*count;
    
    % r(2) = plot(params.t_span, proportion_asymp*ones(size(params.t_span)),'k--','LineWidth',2); hold on;
    axis([0 params.t_span(end) 0 1]);
    xlabel('Time (days)'); ylabel({'Proportion'; 'asymptomatic'; 'transmission'});
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    
    
    figure(1); subplot(4,1,3);
    this_p = plot(params.t_span, this_proportion_asymp_incidence,'Color',cbf_colors,'LineWidth',2); hold on;
    %     this_p.Color(4) = 1-increment*count;
    % r(2) = plot(params.t_span, proportion_asymp*ones(size(params.t_span)),'k--','LineWidth',2); hold on;
    axis([0 params.t_span(end) 0 1]);
    xlabel('Time (days)'); ylabel({'Proportion'; 'asymptomatic'; 'incidence'});
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    
    
    figure(1); subplot(4,1,4);
    this_p = semilogy(params.t_span,this_Rt_assortmixing,'Color',cbf_colors,'LineWidth',2); hold on;
    %     this_p.Color(4) = 1-increment*count;
    axis([0 params.t_span(end) 10^-1 10^1]);
    xlabel('Time (days)'); ylabel({'Effective'; 'reproduction'; 'number'});
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    
    
end




%%
% save simulated data
if save_ans==1
    folder_location = '../../Code_plt_ms_figures/sim_data/';
    save(strcat(folder_location,filename),'params_collect','results_collect');
    
    fprintf('Saved to file: \n'); % want to be close to 25 days in
    fprintf(strcat(filename,'\n'));
    
end
