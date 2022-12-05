
%% simulate SEIR model with fixed proportion of asymptomatic incidence, p
% same reproduction numbers, R0_a = R0_s

clear all; close all; clc;

%% want to save?
save_ans = 0;
% 0: don't save
% 1: save


%% include correlations?
include_correlations= 0;
% 0: don't save
% 1: save

%% which set of time scales?
relR0 = 1;
params.relR0 = relR0;
% 1: Rs = Ra
% 4: Rs = 4*Ra

%% set up colors and parameters
cbf_colors_gray = [0.5 0.5 0.5];
cbf_colors_black = [0 0 0];
cbf_colors_vector = [cbf_colors_black;cbf_colors_gray];

if include_correlations == 1
    if relR0 == 4
        filename = 'SEIR_assortmixing_twodiseases_Rs4timesRa_varyTa_120222.mat';
        
    else
        filename = 'SEIR_assortmixing_twodiseases_sameR0s_varyTa_120222.mat';
    end
    
else
    if relR0 == 4
        filename = 'SEIR_fixedp_twodiseases_Rs4timesRa_varyTa_120222.mat';
    else
        filename = 'SEIR_fixedp_twodiseases_sameR0s_varyTa_120222.mat';
    end
    
end



%% set up colors and parameters
cbf_colors = [0.5,0.5,0.5];


%% vary relative infectious periods, decay rates, days^-1
n_pts = 4;
% k_vector_relT_Tavary = linspace(1,8/5,n_pts); % ratio of Ta/Ts - Ts fixed
% k_vector_relT_Tsvary = linspace(5/8,1,n_pts); % ratio of Ta/Ts - Tas fixed
k_vector_relT_Tavary = [linspace(1,6/5,n_pts/2),linspace(6/5,8/5,n_pts/2)]; % ratio of Ta/Ts - Ts fixed
k_vector_relT_Tsvary = [linspace(5/8,5/6,n_pts/2),linspace(5/6,1,n_pts/2)]; % ratio of Ta/Ts - Ts fixed

k_vector_relT = [k_vector_relT_Tsvary(1:(end-1)),k_vector_relT_Tavary];
% k_vector_relT = [k_vector_relT_Tsvary(1:(end-1)),k_vector_relT_Tavary(2:end)];
% k_vector_relT = [k_vector_relT_Tsvary,k_vector_relT_Tavary];
params.k_vector_relT = k_vector_relT;

Ts_vector_fixed=5*ones(size(k_vector_relT_Tavary));
% Ts_vector_fixed=5*ones(size(k_vector_relT_Tavary(2:end)));
% Ta_vector_vary=k_vector_relT_Tavary(2:end).*Ts_vector_fixed;
Ta_vector_vary=k_vector_relT_Tavary.*Ts_vector_fixed;

Ta_vector_fixed=5*ones(size(k_vector_relT_Tsvary(1:(end-1))));
% Ta_vector_fixed=5*ones(size(k_vector_relT_Tsvary));
Ts_vector_vary=Ta_vector_fixed./k_vector_relT_Tsvary(1:(end-1));
% Ts_vector_vary=Ta_vector_fixed./k_vector_relT_Tsvary;

Ta_vector = [Ta_vector_fixed,Ta_vector_vary];
Ts_vector = [Ts_vector_vary,Ts_vector_fixed];


%% burnin time depends on parameters
t_end_burnin = 71.77;

% parameters
gamma_e=1/3; % 3 day exposure period

% params.beta_a = beta_a_init;
% params.beta_s =beta_s_init;
params.gamma_e = gamma_e;

% p is the proportion of asymptomatic incidence
if include_correlations==1
    
    p_aa = 0.5; p_as = 0.25;
    params.p_aa = p_aa; params.p_as = p_as;
    
else
    
    proportion_asymp = 0.4;
    params.p = proportion_asymp;
    p_aa = proportion_asymp; p_as = proportion_asymp;
    params.p_aa = p_aa; params.p_as = p_as;
    
end

% proportion_asymp_incidence_intrinsic_vary = [0.681622900724801;0.672884566673193;0.672884566673193;0.666666666666667;0.660370550381307;0.660370550381299;0.649875583725111];

% fixed exponential growth rate
fixed_r = 0.14;
params.fixed_r = fixed_r;

%% vary relative infectious periods, decay rates, days^-1
n_pts = (length(k_vector_relT)+1)/2;
% k_vector_relT_Tavary = linspace(1,8/5,n_pts); % ratio of Ta/Ts - Ts fixed
% k_vector_relT_Tsvary = linspace(5/8,1,n_pts); % ratio of Ta/Ts - Tas fixed
k_vector_relT_Tavary = [linspace(1,6/5,n_pts/2),linspace(6/5,8/5,n_pts/2)]; % ratio of Ta/Ts - Ts fixed
k_vector_relT_Tsvary = [linspace(5/8,5/6,n_pts/2),linspace(5/6,1,n_pts/2)]; % ratio of Ta/Ts - Ts fixed

k_vector_relT = [k_vector_relT_Tsvary(1:(end-1)),k_vector_relT_Tavary];
% k_vector_relT = [k_vector_relT_Tsvary(1:(end-1)),k_vector_relT_Tavary(2:end)];
% k_vector_relT = [k_vector_relT_Tsvary,k_vector_relT_Tavary];
params.k_vector_relT = k_vector_relT;

Ts_vector_fixed=5*ones(size(k_vector_relT_Tavary));
% Ts_vector_fixed=5*ones(size(k_vector_relT_Tavary(2:end)));
% Ta_vector_vary=k_vector_relT_Tavary(2:end).*Ts_vector_fixed;
Ta_vector_vary=k_vector_relT_Tavary.*Ts_vector_fixed;

Ta_vector_fixed=5*ones(size(k_vector_relT_Tsvary(1:(end-1))));
% Ta_vector_fixed=5*ones(size(k_vector_relT_Tsvary));
Ts_vector_vary=Ta_vector_fixed./k_vector_relT_Tsvary(1:(end-1));
% Ts_vector_vary=Ta_vector_fixed./k_vector_relT_Tsvary;

Ta_vector = [Ta_vector_fixed,Ta_vector_vary];
Ts_vector = [Ts_vector_vary,Ts_vector_fixed];


% give an initial guess
beta_a_init = 0.4149;
proportion_asymp_init = 2/3;

%% next fit beta_a to little r
for count=1:length(k_vector_relT)
    
    
    % gammas from infectious periods: Ts=5, Ta=5,...,8
    gamma_a = 1/Ta_vector(count); gamma_s = 1/Ts_vector(count);
    params.gamma_a = gamma_a;
    params.gamma_s = gamma_s;
    
    %% inital conditions to findint betas such that r=0.14
    x0=beta_a_init;
    beta_s_init = relR0*(beta_a_init/gamma_a)*gamma_s;
    
    fprintf('finding minimum wrt transmission rates... \n\n');
    fprintf('count = %d \n\n',count);
    fprintf('Ta/Ts = %2.4f \n\n',k_vector_relT(count));
    
    options = optimset('TolFun',10^-14,'TolX',10^-14);
    [x_soln,f_val] = fminsearch(@(x)growthrate_objective_function_general(x,params),x0,options);
    
    beta_a = x_soln(1);
    beta_s = relR0*(beta_a/gamma_a)*gamma_s;
    
    % initialize for next point
    beta_a_init = beta_a;
    
    params.beta_a = beta_a;
    params.beta_s = beta_s;
    
    fprintf('beta_a =  %2.4f \n\n',beta_a);
    fprintf('beta_s =  %2.4f \n\n',beta_s);
    
    % bestfit_rsquared = 1-(x_soln-fixed_r)^2;
    % fprintf('best fit r^2 =  %2.4f \n\n',bestfit_rsquared);
    
    bestfit_SSE = f_val;
    fprintf('best fit SSE =  %1.2e \n\n',bestfit_SSE);
    
    results.bestfit_SSE_betas=bestfit_SSE;
    
    % need to get eigen proportion direction
    eigen_direction_assortmixing = get_eigendirection_SEIR_twodiseases_assortmixing(params);
    
    
    %% need to find the corresponding intrinsic proportion, p, to find z
    % calculate the intrinsic proportion of asymptomatic transmission
    if include_correlations == 1
        if k_vector_relT(count)==1
            
            proportion_asymp = 2/3;
            params.p = proportion_asymp;
            
            % initialize for next point
            proportion_asymp_init = proportion_asymp;
            
            fprintf('proportion asymptomatic =  %2.4f \n\n',proportion_asymp);
            
            
        else
            x0=proportion_asymp_init;
            
            fprintf('finding minimum wrt p... \n\n');
            
            options = optimset('TolFun',10^-14,'TolX',10^-14);
            [x_soln,f_val] = fminsearch(@(x)growthrate_objective_function_propasymp(x,params),x0,options);
            
            proportion_asymp = x_soln(1);
            params.p = proportion_asymp;
            
            % initialize for next point
            proportion_asymp_init = proportion_asymp;
            
            fprintf('proportion asymptomatic =  %2.4f \n\n',proportion_asymp);
            
            bestfit_SSE = f_val;
            fprintf('best fit SSE =  %1.2e \n\n',bestfit_SSE);
            
            results.bestfit_SSE_p=bestfit_SSE;
            
        end
        

    end
    
    %% now calculate intrinsic proportion asymp transmission
    %     total_transmission_z = p_aa*beta_a/gamma_a+(1-p_aa)*beta_s/gamma_s+p_as*beta_a/gamma_a+(1-p_as)*beta_s/gamma_s;
    %     asymp_transmission_z = p_aa*beta_a/gamma_a+p_as*beta_s/gamma_s;
    total_transmission_z = proportion_asymp*beta_a/gamma_a+(1-proportion_asymp)*beta_s/gamma_s;
    asymp_transmission_z = proportion_asymp*beta_a/gamma_a;
    
    proportion_asymp_transmission_z = asymp_transmission_z/total_transmission_z;
    results.proportion_asymp_transmission_z=proportion_asymp_transmission_z;
    
    % calculate the realized proportion of asymptomatic transmission
    total_transmission_q = beta_a*eigen_direction_assortmixing(4)+beta_s*eigen_direction_assortmixing(5);
    asymp_transmission_q = beta_a*eigen_direction_assortmixing(4);
    proportion_asymp_transmission_q = asymp_transmission_q/total_transmission_q;
    results.proportion_asymp_transmission_q=proportion_asymp_transmission_q;
    
    R0_assortmixing = get_R0_SEIR_twodiseases_assortmixing(params);
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
    this_asymp_incidence = p_aa*beta_a_traj.*(I_a_traj.*S_traj)+p_as*beta_s_traj.*(I_s_traj.*S_traj);
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
    if 1
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
    
    
end




%%
% save simulated data
if save_ans==1
    folder_location = '../../Code_plt_ms_figures/sim_data/';
    save(strcat(folder_location,filename),'params_collect','results_collect');
    
    fprintf('Saved to file: \n'); % want to be close to 25 days in
    fprintf(strcat(filename,'\n'));
    
end
