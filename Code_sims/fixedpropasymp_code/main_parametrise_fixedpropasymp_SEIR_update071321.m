
%% now vary R0 vs. r

clear all; close all; clc;

%% save Figure ?
save_ans_Fig = 0;
% 0: don't save
% 1: save

figure_name = 'FigS5_parametrise_fixedpropasymp_071321';

% dark blue, violet, light blue
cat_cbf_colors = [15,32,128;169,90,161;133,192,249]/255;

%%%% vary R0 vs. r
lower_boundR0 = 1;
upper_boundR0 = 3;
n_pts=501;
vary_R0 = linspace(lower_boundR0,upper_boundR0,n_pts);

% set levels to find
level_R0 = 2.4;
level_r = 0.14;

% exposed period
gamma_e = 1/3;
params.gamma_e = gamma_e;

% symptomatic infectious period, Ts = 5 days
gamma_s=1/5;
params.gamma_s = gamma_s;

% vary asymptomatic infectious period, Ts = 5 days
gamma_a1=1/5; gamma_a2=1/6; gamma_a3=1/8;
vary_gamma_a = [gamma_a1 gamma_a2 gamma_a3];


%%%%% one p %%%%%
proportion_asymp = 0.4;
params.p = proportion_asymp;

%% same R0s
for count_row = 1:length(vary_gamma_a)
    
    this_gamma_a = vary_gamma_a(count_row);
    params.gamma_a = this_gamma_a;
    
    for count_col = 1:length(vary_R0)
        
        this_R0 = vary_R0(count_col);
        this_R0_a = this_R0;
        this_R0_s = this_R0;
        this_beta_a = this_R0_a*(params.gamma_a);
        params.beta_a = this_beta_a;
        this_beta_s = this_R0_s*(params.gamma_s);
        params.beta_s =this_beta_s;
        
        
        %%% need to get eigen proportion direction
        eigen_direction_fixedpropasymp = get_eigendirection_SEIR_twodiseases_fixedpropasymp(params);
        total_incidence(count_row,count_col) = this_beta_a*eigen_direction_fixedpropasymp(4)+this_beta_s*eigen_direction_fixedpropasymp(5);
        asymp_incidence(count_row,count_col) = (params.p)*total_incidence(count_col);
        asymp_prop_incidence(count_row,count_col) = asymp_incidence(count_col)/total_incidence(count_col);
        
        r_fixedpropasymp(count_row,count_col) = get_r_SEIR_twodiseases_fixedpropasymp(params);
        
        
    end
    
    
    ind_r = find(r_fixedpropasymp(count_row,:)>level_r);
    R0_a_fixedpropasymp(count_row) = vary_R0(ind_r(1));
    R0_s_fixedpropasymp(count_row) = vary_R0(ind_r(1));
    beta_a_fixedpropasymp(count_row) = R0_a_fixedpropasymp(count_row)*params.gamma_a;
    beta_s_fixedpropasymp(count_row) = R0_s_fixedpropasymp(count_row)*params.gamma_s;
    level_r_fixedpropasymp(count_row) = r_fixedpropasymp(count_row,ind_r(1));
    
    params.beta_a = beta_a_fixedpropasymp(count_row);
    params.beta_s = beta_s_fixedpropasymp(count_row);
    R0_fixedpropasymp(count_row) = get_R0_SEIR_twodiseases_fixedpropasymp(params);
    
    
    
    f1 = figure(1); set(gcf,'Position',[240   300   1120   420]);
    subplot(1,2,1);
    p(2*count_row-1) = plot(vary_R0, r_fixedpropasymp(count_row,:), 'Color',cat_cbf_colors(count_row,:),'LineWidth',2); hold on;
    p(2*count_row) = plot(R0_a_fixedpropasymp(count_row), level_r_fixedpropasymp(count_row), 'k.','Color',cat_cbf_colors(count_row,:),'MarkerSize',20); hold on;
    axis([0.5 3 0 0.2]);
    xlabel('$\mathcal{R}_0 \,(= \mathcal{R}_{0,a} = \mathcal{R}_{0,s})$','Interpreter','latex');
    ylabel('Exponential Growth Rate, $r$','Interpreter','latex');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
end

leg_char1=[' $T_a = ', num2str(1/gamma_a1),'$'];
leg_char2=['$\beta_a = ', num2str(beta_a_fixedpropasymp(1),'%2.2f'),', \beta_s = ', num2str(beta_s_fixedpropasymp(1),'%2.2f'),', \mathcal{R}_0 = ', num2str(R0_fixedpropasymp(1),'%2.2f'),'$'];
leg_char3=['$T_a = ', num2str(1/gamma_a2),'$'];
leg_char4=['$\beta_a = ', num2str(beta_a_fixedpropasymp(2),'%2.2f'),', \beta_s = ', num2str(beta_s_fixedpropasymp(2),'%2.2f'), ', \mathcal{R}_0 = ', num2str(R0_fixedpropasymp(2),'%2.2f'),'$'];
leg_char5=['$T_a = ', num2str(1/gamma_a3),'$'];
leg_char6=['$\beta = ', num2str(beta_a_fixedpropasymp(3),'%2.2f'),', \beta_s = ', num2str(beta_s_fixedpropasymp(3),'%2.2f'),', \mathcal{R}_0 = ', num2str(R0_fixedpropasymp(3),'%2.2f'),'$'];
legend(p,{leg_char1,leg_char2,leg_char3,leg_char4,leg_char5,leg_char6}, 'Location','NorthWest','Interpreter','Latex');
legend boxoff

txt = {'A'};
text(0.025,1.025,txt,'Units','normalized','FontSize',14,'FontWeight','bold');


%% same betas
% now vary beta vs. r

upper_boundbeta = 0.8;
lower_boundbeta = 0;
vary_beta_a = linspace(lower_boundbeta,upper_boundbeta,n_pts);



for count_row = 1:length(vary_gamma_a)
    
    this_gamma_a = vary_gamma_a(count_row);
    params.gamma_a = this_gamma_a;
    
    for count_col = 1:length(vary_beta_a)
        
        this_beta_a = vary_beta_a(count_col);
        params.beta_a = this_beta_a;
        this_beta_s = this_beta_a;
        params.beta_s =this_beta_s;
        
        
        %%% need to get eigen proportion direction
        eigen_direction_fixedpropasymp = get_eigendirection_SEIR_twodiseases_fixedpropasymp(params);
        total_incidence(count_row,count_col) = this_beta_a*eigen_direction_fixedpropasymp(4)+this_beta_a*eigen_direction_fixedpropasymp(5);
        asymp_incidence(count_row,count_col) = (params.p)*total_incidence(count_col);
        asymp_prop_incidence(count_row,count_col) = asymp_incidence(count_col)/total_incidence(count_col);
        
        r_fixedpropasymp(count_row,count_col) = get_r_SEIR_twodiseases_fixedpropasymp(params);
        %         R0_fixedpropasymp(n) = get_R0_SEIR_twodiseases_fixedpropasymp(params);
        
    end
    
    
    ind_r = find(r_fixedpropasymp(count_row,:)>level_r);
    beta_r_fixedpropasymp(count_row) = vary_beta_a(ind_r(1));
    level_r_fixedpropasymp(count_row) = r_fixedpropasymp(count_row,ind_r(1));
    
    params.beta_a = beta_r_fixedpropasymp(count_row);
    params.beta_s = beta_r_fixedpropasymp(count_row);
    R0_fixedpropasymp(count_row) = get_R0_SEIR_twodiseases_fixedpropasymp(params);
    
    
    f1 = figure(1); set(gcf,'Position',[240   300   1120   420]);
    subplot(1,2,2);
    q(2*count_row-1) = plot(vary_beta_a, r_fixedpropasymp(count_row,:), 'Color',cat_cbf_colors(count_row,:),'LineWidth',2); hold on;
    q(2*count_row) = plot(beta_r_fixedpropasymp(count_row), level_r_fixedpropasymp(count_row), 'k.','Color',cat_cbf_colors(count_row,:),'MarkerSize',20); hold on;
    axis([0.1 0.6 0 0.2])
    xlabel('$\beta\,(= \beta_a = \beta_s)$','Interpreter','latex');
    ylabel('Exponential Growth Rate, $r$','Interpreter','latex');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
end


leg_char1_q=['$T_a = ', num2str(1/gamma_a1),'$'];
leg_char2_q=['$\beta = ', num2str(beta_r_fixedpropasymp(1),'%2.2f'),', \mathcal{R}_0 = ', num2str(R0_fixedpropasymp(1),'%2.2f'),'$'];
leg_char3_q=['$T_a = ', num2str(1/gamma_a2),'$'];
leg_char4_q=['$\beta = ', num2str(beta_r_fixedpropasymp(2),'%2.2f'),', \mathcal{R}_0 = ', num2str(R0_fixedpropasymp(2),'%2.2f'),'$'];
leg_char5_q=['$T_a = ', num2str(1/gamma_a3),'$'];
leg_char6_q=['$\beta = ', num2str(beta_r_fixedpropasymp(3),'%2.2f'),', \mathcal{R}_0 = ', num2str(R0_fixedpropasymp(3),'%2.2f'),'$'];
legend(q,{leg_char1_q,leg_char2_q,leg_char3_q,leg_char4_q,leg_char5_q,leg_char6_q}, 'Location','NorthWest','Interpreter','Latex');
legend boxoff

txt = {'B'};
text(0.025,1.025,txt,'Units','normalized','FontSize',16,'FontWeight','bold');



%% save figure
if save_ans_Fig
    
    folder_location = './../../Figures_ms_all/supp/';
    saveas(f1,strcat(folder_location,figure_name),'epsc');
    
    fprintf('Figure saved:\n');
    fprintf(strcat(figure_name,'\n\n'));
    
    fprintf('Location:\n');
    fprintf(strcat(folder_location,'\n\n'));
    
else
    
    fprintf('Figure not saved.\n');
end