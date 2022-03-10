
%% Figure 3
% A-E: no mit, increase time-scale differences
% F-J: with mit, increase time-scale differences
% simulate the two disease SIR model with two infectious compartments:
% asymptomatic and symptomatic infectious persons

clear all; close all; clc;


%% save figure?
save_ans = 0;
% 0: don't save
% 1: save

% viridis color palette
colors_rgb_vector = [253, 231, 37;...
    160, 218, 57;...
    74, 193, 109;...
    31, 161, 135;...
    39, 127, 142;...
    39, 127, 142;...
    70, 50, 126;...
    68, 1, 84]/255;


frac_spacing = 0.5;
frac_scaling = 0.4;


fprintf('No mitigation... \n\n');

switch_over_var = [1,2]; % 1 or 2

%% load no mitigation files
for counter=1
    
    which_var = switch_over_var(counter);
    
    switch which_var
        
        case 1
            % Ta = Ts = 5 (regular assortivity)
            infile = 'SEIR_agedep_twodiseases_variationsymptomaticity_T5and5_mit.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [15,32,128]/255; % dark blue
            gamma_a1 = params.gamma_a;
            gamma_s1 = params.gamma_s;
            beta_a1 = params.beta_a;
            beta_s1 = params.beta_s;
            R0_agedep1 = results.Rt_agedep(1);
            supfig_name = 'Shanghai: intervention'
            figure_name = 'Figure_Sn_Shanghai_intervention_combine';

            
        case 2
            
            infile = 'SEIR_agedep_twodiseases_variationsymptomaticity_T5and8_mit.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [133,192,249]/255; % light blue
            gamma_a2 = params.gamma_a;
            gamma_s2 = params.gamma_s;
            beta_a2 = params.beta_a;
            beta_s2 = params.beta_s;
            R0_agedep2 = results.Rt_agedep(1);
            %             params.variation_symptomaticity_yesno=1; % to get dashed line
%             supfig_name = 'Shanghai: different G.I.'
%             figure_name = 'Figure_Sn_Shanghai_diffGI';
    end
    
    
    fprintf('Opened file: \n');
    fprintf(strcat(infile,'\n\n'));
    
    fprintf('Infectious periods: \n');
    fprintf('T_a = %1i days \n',1/params.gamma_a);
    fprintf('T_s = %1i days \n\n',1/params.gamma_s);
    
    fprintf('Basic reproductive number \n');
    fprintf('R_0 =  %2.2f \n\n',results.Rt_agedep(1));
    
    fprintf('Exponential growth rate \n');
    fprintf('r =  %2.2f \n\n',results.r_agedep);
    
    fprintf('Initial proportion asymptomatic incidence \n');
    fprintf('%2.2f \n\n',results.proportion_asymptomatic_incidence(1));
    
    ind = find(results.i_all_total_fraction>10^-6);
    fprintf('Time of first appearance: \n'); % want to be close to 25 days in
    fprintf('%2.2f days \n\n',params.t_span(ind(1)));
    
    fprintf('-------------------------------------- \n\n');
    
    %% plot panels A-E: no mitigation
    f1 = figure(1); set(f1, 'Position', [100 500 1000 650]);
    sgtitle(supfig_name);
    subplot(2,2,1);
    if params.variation_symptomaticity_yesno==0
        this_p = semilogy(params.t_span, results.i_all_total_fraction,'Color',cbf_colors,'LineWidth',2); hold on;
        this_p.Color(4) = 1-0.18*2;
    else
        semilogy(params.t_span, results.i_all_total_fraction,'--','Color',cbf_colors,'LineWidth',2); hold on;
    end
    axis([0 params.t_span(end) 10^(-6) 1]);
    %     xlabel('Time (days)');
    ylabel({'Total'; 'Incidence, $i(t)$'},'Interpreter','Latex');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    title('Equal Generation Intervals','FontSize',16);
    
    %     if counter==1
    txt = {'A'};
    text(0.025,0.925,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    
    old_pos = get(f1, 'Position');
    set(f1,'Position',[old_pos(1),frac_spacing,old_pos(3),frac_scaling]);
    old_pos = get(f1, 'Position');
    
    yticks([10^-6 10^-4 10^-2 10^0 ]);
    set(f1,'xticklabel',{[]},'yticklabel',[{' '},{'10^{-4}'},{'10^{-2}'},{'10^0'}]);
    box('off');
    
    txt = {'10^{-6}'};
    text(-0.11,0.05,txt,'Units','normalized',...
        'FontSize',14,'FontWeight','normal','FontName', 'Times');
    %     end
    
    
    
    figure(1); subplot(2,2,3);
    %         this_p = plot(params.t_span, results.ave_age_i_all,'Color',cbf_colors,'LineWidth',2); hold on;
    for count=1:8
        %         this_p = plot(params.t_span, results.ave_age_i_all,'Color',colors_rgb_vector,'LineWidth',2); hold on;
        this_h(count) = semilogy(params.t_span, results.y_traj(:,count)/data_Shanghai_Davies.age_distribution(count),'Color',colors_rgb_vector(count,:),'LineWidth',2); hold on;
        this_h(count).Color(4) = 1-0.18*2.5*(count/8);
        %     plot(params.t_span, results.ave_age_i_all,'--','Color',colors_rgb_vector,'LineWidth',2); hold on;
    end
    axis([0 params.t_span(end) 0.991 1]);
    %     xlabel('Time (days)');
    ylabel({'$S_n/A_n$'},'Interpreter','Latex');
    yticks([0.991 1]);
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    
    h_legend_1 = ['$S_1$'];
    h_legend_2 = ['$S_2$'];
    h_legend_3 = ['$S_3$'];
    h_legend_4 = ['$S_4$'];
    h_legend_5 = ['$S_5$'];
    h_legend_6 = ['$S_6$'];
    h_legend_7 = ['$S_7$'];
    h_legend_8 = ['$S_8$'];
    
    legend(this_h,{h_legend_1,h_legend_2,h_legend_3,h_legend_4,h_legend_5,h_legend_6,h_legend_7,h_legend_8}, 'Interpreter','Latex','Location','SouthWest','FontSize',10.5);
    legend boxoff
    
    %     if counter==2
    txt = {'B'};
    text(0.025,0.925,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    
    set(f1,'Position',[old_pos(1), old_pos(2)-frac_scaling, old_pos(3), old_pos(4)])
    old_pos = get(f1, 'Position');
    
    set(f1,'yticklabel',[{''},{''}]);
    box('off');
    
    
    txt = {'0.991'};
    text(-0.095,0.045,txt,'Units','normalized',...
        'FontSize',14,'FontWeight','normal','FontName', 'Times');
    
    txt = {'1'};
    text(-0.07,0.93,txt,'Units','normalized',...
        'FontSize',14,'FontWeight','normal','FontName', 'Times');
    
    
    
end



%%
fprintf('With mitigation... \n\n');

%% load with mitigation files
for counter=2
    
    which_var = switch_over_var(counter);
    
    switch which_var
        
        case 1
            % Ta = Ts = 5 (regular assortivity)
            infile = 'SEIR_agedep_twodiseases_variationsymptomaticity_T5and5_mit.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [15,32,128]/255; % dark blue
            gamma_a1 = params.gamma_a;
            gamma_s1 = params.gamma_s;
            beta_a1 = params.beta_a;
            beta_s1 = params.beta_s;
            R0_agedep1 = results.Rt_agedep(1);
            
            
        case 2
            % Ta = Ts = 5 (4x assortivity)
            infile = 'SEIR_agedep_twodiseases_variationsymptomaticity_T5and8_mit.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [133,192,249]/255; % light blue
            gamma_a2 = params.gamma_a;
            gamma_s2 = params.gamma_s;
            beta_a2 = params.beta_a;
            beta_s2 = params.beta_s;
            R0_agedep2 = results.Rt_agedep(1);
            params.variation_symptomaticity_yesno=1; % to get dashed line
            
    end
    
    
    fprintf('Opened file: \n');
    fprintf(strcat(infile,'\n\n'));
    
    fprintf('Infectious periods: \n');
    fprintf('T_a = %1i days \n',1/params.gamma_a);
    fprintf('T_s = %1i days \n\n',1/params.gamma_s);
    
    fprintf('Basic reproductive number \n');
    fprintf('R_0 =  %2.2f \n\n',results.Rt_agedep(1));
    
    fprintf('Exponential growth rate \n');
    fprintf('r =  %2.2f \n\n',results.r_agedep);
    
    fprintf('Initial proportion asymptomatic incidence \n');
    fprintf('%2.2f \n\n',results.proportion_asymptomatic_incidence(1));
    
    ind = find(results.i_all_total_fraction>10^-6);
    fprintf('Time of first appearance: \n'); % want to be close to 25 days in
    fprintf('%2.2f days \n\n',params.t_span(ind(1)));
    
    fprintf('-------------------------------------- \n\n');
    
    
    
    %% plot panels E-I: with mitigation
    figure(1);
    subplot(2,2,2);
    if params.variation_symptomaticity_yesno==0
        semilogy(params.t_span, results.i_all_total_fraction,'Color',cbf_colors,'LineWidth',2); hold on;
        
    else
        semilogy(params.t_span, results.i_all_total_fraction,'--','Color',cbf_colors,'LineWidth',2); hold on;
        
    end
    axis([0 params.t_span(end) 10^(-6) 1]);
    
    ylabel({'Total'; 'Incidence, $i(t)$'},'Interpreter','Latex');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    title('Different Generation Intervals','FontSize',16);
    
    %     if counter==2
    txt = {'C'};
    text(0.025,0.925,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    
    old_pos = get(f1, 'Position');
    set(f1,'Position',[old_pos(1),frac_spacing,old_pos(3),frac_scaling]);
    old_pos = get(f1, 'Position');
    
    yticks([10^-6 10^-4 10^-2 10^0 ]);
    set(f1,'xticklabel',{[]},'yticklabel',[{' '},{'10^{-4}'},{'10^{-2}'},{'10^0'}]);
    box('off');
    
    txt = {'10^{-6}'};
    text(-0.11,0.05,txt,'Units','normalized',...
        'FontSize',14,'FontWeight','normal','FontName', 'Times');
    
    
    
    
    
    
    
    figure(1); subplot(2,2,4);
    for count=1:8
        %         this_p = plot(params.t_span, results.ave_age_i_all,'Color',colors_rgb_vector,'LineWidth',2); hold on;
        this_p(count)=semilogy(params.t_span, results.y_traj(:,count)/data_Shanghai_Davies.age_distribution(count),'Color',colors_rgb_vector(count,:),'LineWidth',2); hold on;
        this_p(count).Color(4) = 1-0.18*2.5*(count/8);
        %         this_p.Color(4) = 1-0.18*2;
        %     plot(params.t_span, results.ave_age_i_all,'--','Color',colors_rgb_vector,'LineWidth',2); hold on;
    end
    
%     axis([0 params.t_span(end) 35 45]);
    axis([0 params.t_span(end) 0.991 1]);
    xlabel('Time (days)');
    ylabel({'$S_n/A_n$'},'Interpreter','Latex');
    yticks([0.991 1]);
    
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    %     if counter==2
    txt = {'D'};
    text(0.025,0.925,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    
    set(f1,'Position',[old_pos(1), old_pos(2)-frac_scaling, old_pos(3), old_pos(4)])
    old_pos = get(f1, 'Position');
    
    set(f1,'yticklabel',[{''},{''}]);
    box('off');
    
    
    txt = {'0.991'};
    text(-0.095,0.045,txt,'Units','normalized',...
        'FontSize',14,'FontWeight','normal','FontName', 'Times');
    
    txt = {'1'};
    text(-0.07,0.93,txt,'Units','normalized',...
        'FontSize',14,'FontWeight','normal','FontName', 'Times');
    
    
end






%% save figure
if save_ans
    
    folder_location = './../Figures_ms_all/supp/';
    saveas(f1,strcat(folder_location,figure_name),'epsc');
    
    fprintf('Figure saved:\n');
    fprintf(strcat(figure_name,'\n\n'));
    
    fprintf('Location:\n');
    fprintf(strcat(folder_location,'\n\n'));
    
end
