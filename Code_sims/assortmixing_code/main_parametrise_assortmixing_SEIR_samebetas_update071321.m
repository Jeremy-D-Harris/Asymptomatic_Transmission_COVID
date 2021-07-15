
%% paramtrize model wrt p(a|a) and p(a|s) to find:
% (1) initial proportion asymptomatic incidence
% (2) exponential growth rate
% beta_a = beta_s = beta

clear all; close all; clc;


%% want to save?
save_ans = 0;
% 0: don't save
% 1: save

%% want to save?
plt_figure = 1;
% 0: don't plot
% 1: plot


%% which set of time scales?
which_timescales = 1; % 1,2
% 1: same time scales: Ta=Ts=5 days
% 2: longer time scales of asymptomatic transmission: Ta=8,Ts=5 days

n_pts = 512; % more points means finer mesh

if which_timescales==1
    filename = 'parametrise_assortmixing_same_samebetas_SEIR_071321';
    
    gamma_a=1/5; gamma_s=1/5;
    beta_a = 0.4840;
    beta_s = beta_a;
    
else
    filename = 'parametrise_assortmixing_diff_samebetas_SEIR_071321';
    
    %% different time scales
    gamma_a=1/8; gamma_s=1/5;
    beta_a = 0.43467;
    beta_s = beta_a;
    
end

gamma_e = 1/3;

params.beta_a = beta_a;
params.beta_s =beta_s;
params.gamma_a = gamma_a;
params.gamma_s = gamma_s;
params.gamma_e = gamma_e;

%% range over p(a|a) and p(a|s)
upper_bound = 1;
lower_bound = 0;
params.upper_bound=upper_bound;
params.lower_bound=lower_bound;

vary_p_aa = linspace(lower_bound,upper_bound,n_pts);
vary_p_as = linspace(lower_bound,upper_bound,n_pts);
params.vary_p_aa=vary_p_aa;
params.vary_p_as=vary_p_as;

% set levels to find
level_prop_asymp_incidence = 0.4; % initial proportion asymptomatic incidence, p_e
level_r = 0.14; % exponential growth rate, r
params.level_prop_asymp_incidence=level_prop_asymp_incidence;
params.level_r=level_r;

perturb = 1e-11;

t_end_burnin =100;
t_start = 0; t_end = t_end_burnin; % burn in time

dt=0.01;
params.dt=dt;
params.t_span = t_start:dt:t_end;

%% vary over p(a|a) and p(a|s)
for count_row = 1:length(vary_p_aa)
    
    count_row
    this_p_aa = vary_p_aa(count_row);
    params.p_aa = this_p_aa;
    
    for count_col =1:length(vary_p_as)
        
        this_p_as = vary_p_as(count_col);
        params.p_as =this_p_as;
        
        eigen_direction_assortmixing = get_eigendirection_SEIR_twodiseases_assortmixing(params);
        r_assortmixing(count_row,count_col) = get_r_SEIR_twodiseases_assortmixing(params);
        
        
        if eigen_direction_assortmixing(1)>0
            init_conds = [1;0;0;0;0;0;0] - perturb*eigen_direction_assortmixing;
        else
            init_conds = [1;0;0;0;0;0;0] + perturb*eigen_direction_assortmixing;
        end
        
        S_traj_end = init_conds(1);
        I_a_traj_end = init_conds(4); I_s_traj_end = init_conds(5);
        
        init_asymp_incidence = params.p_aa*params.beta_a*S_traj_end*I_a_traj_end+params.p_as*params.beta_s*S_traj_end*I_s_traj_end;
        init_total_incidence = params.beta_a*S_traj_end*I_a_traj_end+beta_s*S_traj_end*I_s_traj_end;
        
        init_prop_asymp_incidence(count_row,count_col) = init_asymp_incidence/init_total_incidence;
        
        if count_row>count_col
            r_assortmixing_plt(count_row,count_col) = r_assortmixing(count_row,count_col);
            eigen_prop_asymp_incidence_plt(count_row,count_col) = init_prop_asymp_incidence(count_row,count_col);
        else
            r_assortmixing_plt(count_row,count_col) = 0;
            eigen_prop_asymp_incidence_plt(count_row,count_col) = 0;
            
        end
        
        
    end
    
end


%% need to pad for some reason?
vary_p_aa_plt = [vary_p_aa,vary_p_aa(end)+0.1];
vary_p_as_plt = [vary_p_as,vary_p_as(end)+0.1];
params.vary_p_aa_plt=vary_p_aa_plt;
params.vary_p_as_plt=vary_p_as_plt;

r_assortmixing_plt = [r_assortmixing_plt; zeros(1,length(r_assortmixing_plt(1,:)))];
r_assortmixing_plt = [r_assortmixing_plt, zeros(length(r_assortmixing_plt(:,1)),1)];
results.r_assortmixing_plt=r_assortmixing_plt;

eigen_prop_asymp_incidence_plt = [eigen_prop_asymp_incidence_plt; zeros(1,length(eigen_prop_asymp_incidence_plt(1,:)))];
eigen_prop_asymp_incidence_plt = [eigen_prop_asymp_incidence_plt, zeros(length(eigen_prop_asymp_incidence_plt(:,1)),1)];
results.eigen_prop_asymp_incidence_plt=eigen_prop_asymp_incidence_plt;


if plt_figure==1
    %% plot three panels
    figure(1); set(gcf,'Position',[50 400 1650 350]);
    subplot(1,3,1);
    
    if which_timescales==1
        % x variable: goes down the rows
        % y variable: goes across columns
        s1 = pcolor(params.vary_p_as_plt,params.vary_p_aa_plt,results.r_assortmixing_plt); hold on;
    else
        s1 = pcolor(params.vary_p_as_plt,params.vary_p_aa_plt,results.r_assortmixing_plt); hold on;
        ind_cont_r = contourc(r_assortmixing_plt,[level_r level_r]);
        contour_r = [vary_p_as_plt(round(ind_cont_r(1,2:end))); vary_p_aa_plt(round(ind_cont_r(2,2:end)))];
        
        ind_keep_cont_r = find(contour_r(1,:)<=level_prop_asymp_incidence & contour_r(2,:)< 1); % find P(a|a)>0.4
        contour_r_filtered = contour_r(:,ind_keep_cont_r);
        results.contour_r_filtered = contour_r_filtered;
        plot(results.contour_r_filtered(1,:),results.contour_r_filtered(2,:),'k--','LineWidth',1.5); hold on;
        
    end
    set(s1, 'EdgeColor', 'none');
    axis([params.lower_bound params.upper_bound params.lower_bound params.upper_bound]);
    caxis(gca,[0.1 0.16]);
    cb1 = colorbar;
    set(get(cb1,'title'),'string','$r$','Interpreter','Latex','FontSize',16);
    cb1.Limits = [0.1 0.16];
    cb1.Ticks = linspace(0.1, 0.16, 7);
    cb1.TickLabels = linspace(0.1, 0.16, 7);
    colormap(parula);
    xlabel('$p(a|s)$','Interpreter','Latex'); ylabel('$p(a|a)$','Interpreter','Latex');
    xticks(0:0.2:1); yticks(0:0.2:1);
    title('Exponential Growth Rate');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 16;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    figure(1);
    subplot(1,3,2);
    
    % x variable: goes down the rows
    % y variable: goes across columns
    s2 = pcolor(params.vary_p_as_plt,params.vary_p_aa_plt,results.eigen_prop_asymp_incidence_plt); hold on;
    ind_cont_pe = contourc(eigen_prop_asymp_incidence_plt,[level_prop_asymp_incidence level_prop_asymp_incidence]);
    contour_pe = [vary_p_as_plt(round(ind_cont_pe(1,2:end))); vary_p_aa_plt(round(ind_cont_pe(2,2:end)))];
    
    ind_keep_cont_pe = find(contour_pe(1,:)<level_prop_asymp_incidence & contour_pe(2,:)< 1); % find P(a|s)<0.4
    contour_pe_filtered = contour_pe(:,ind_keep_cont_pe);
    results.contour_pe_filtered=contour_pe_filtered;
    
    plot(results.contour_pe_filtered(1,:),results.contour_pe_filtered(2,:),'k','LineWidth',1.5); hold on;
    
    
    % find increased assortivity with p(a|a) = 0.7
    level_increased_assortativity = 0.7;
    if which_timescales==1
    ind_ia = find(contour_pe_filtered(2,:)>=level_increased_assortativity);
    else
    ind_ia = find(contour_pe_filtered(2,:)<=level_increased_assortativity);
    end
    contour_pe_pas_ia = contour_pe_filtered(1,ind_ia(1));
    contour_pe_paa_ia = contour_pe_filtered(2,ind_ia(1));
    results.contour_pe_pas_ia=contour_pe_pas_ia;
    results.contour_pe_paa_ia=contour_pe_paa_ia;
    
    plot(results.contour_pe_pas_ia,results.contour_pe_paa_ia,'k.','MarkerSize',20); hold on;
    set(s2, 'EdgeColor', 'none');
    axis([params.lower_bound params.upper_bound params.lower_bound params.upper_bound]);
    caxis(gca,[0 1]);
    cb2 = colorbar;
    set(get(cb2,'title'),'string','$p_e$','Interpreter','Latex','FontSize',16);
    cb2.Limits = [0 1];
    % cb2.Ticks = linspace(0, 1, 11);
    % cb2.TickLabels = linspace(0, 1, 11);
    colormap(parula);
    xlabel('$p(a|s)$','Interpreter','Latex'); ylabel('$p(a|a)$','Interpreter','Latex');
    xticks(0:0.2:1); yticks(0:0.2:1);
    title('Initial Prop. Asymp. Incidence');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 16;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    
    figure(1);
    subplot(1,3,3);
    plot(params.vary_p_aa,params.vary_p_aa,'Color',[1/2 1/2 1/2]); hold on;
    if which_timescales==2
        plot(results.contour_r_filtered(1,:),results.contour_r_filtered(2,:),'k--','LineWidth',1.5); hold on;
    end
    plot(results.contour_pe_filtered(1,:),results.contour_pe_filtered(2,:),'k','LineWidth',1.5); hold on;
    plot(results.contour_pe_pas_ia,results.contour_pe_paa_ia,'k.','MarkerSize',20); hold on;
    axis([params.lower_bound params.upper_bound params.lower_bound params.upper_bound]);
    xlabel('$p(a|s)$','Interpreter','Latex'); ylabel('$p(a|a)$','Interpreter','Latex');
    xticks(0:0.2:1); yticks(0:0.2:1);
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 16;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    %     legend(this_p,{'Initial prop asymp, $p_e = 0.4$'},'Interpreter','Latex','Location','SouthEast');
    % legend boxoff
    
    
    
    
end


%% save figure
if save_ans
    
    folder_location = './../../Code_plt_ms_figures/sim_data/';
    save(strcat(folder_location,filename),'params','results');
    
    fprintf('File saved:\n'); 
    fprintf(strcat(filename,'\n\n'));
    
    fprintf('Location:\n'); 
    fprintf(strcat(folder_location,'\n\n'));
    
else
    
    fprintf('File not saved.\n'); 
    
end
