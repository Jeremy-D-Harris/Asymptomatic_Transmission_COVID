

%% Figure 2 - parametrised with R_0,a = R_0,s
% A-D: no mitigation, increase time-scale differences
% E-H: with mitigation, increase time-scale differences

clear all; close all; clc;

%% save figure?
save_ans_Fig = 0;
% 0: don't save
% 1: save

figure_name = 'Figure2_correlations_sameR0s_fix_120122';

%% load no mitigation files

switch_over_var = [1,2,3,4];

frac_spacing = 0.74;
frac_scaling = 0.2;

fprintf('No mitigation... \n\n');

% A-D: no mitigation
for counter=1:length(switch_over_var)
    
    this_file = switch_over_var(counter);
    
    switch this_file
        
        case 1
            % NO assortative mixing
            % load: Ta=5 (Ts=5)
            infile = 'SEIR_fixedpropasymp_twodiseases_sameR0s_110922_T5and5.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [15,32,128]/255; % dark blue
            gamma_s1 = params.gamma_s;
            gamma_a1 = params.gamma_a;
            p1 = params.p;
            
            case 2
            % NO assortative mixing
            % load: Ta=8 (Ts=5)
            infile = 'SEIR_fixedpropasymp_twodiseases_sameR0s_110922_T5and8.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [133,192,249]/255; % light blue
            gamma_s2 = params.gamma_s;
            gamma_a2 = params.gamma_a;
            p2 = params.p;
            
            
        case 3
            % WITH assortative mixing
            % load: Ta=5 (same),
            infile = 'SEIR_assortmixing_twodiseases_sameR0s_120122_fix_T5and5.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [15,32,128]/255; % dark blue
            gamma_s3 = params.gamma_s;
            gamma_a3 = params.gamma_a;
            p_aa3=params.p_aa; p_as3=params.p_as; 
            
        case 4
            % WITH assortative mixing
            % load: Ta=8 (Ta=5)
            infile = 'SEIR_assortmixing_twodiseases_sameR0s_120122_fix_T5and8.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [133,192,249]/255; % light blue
            gamma_s4 = params.gamma_s;
            gamma_a4 = params.gamma_a;
            p_aa4=params.p_aa; p_as4=params.p_as; 
            
    end
    
    fprintf('Opened file: \n'); % want to be close to 25 days in
    fprintf(strcat(infile,'\n\n'));
    
    fprintf('Infectious periods: \n');
    fprintf('T_a = %1i days \n',1/params.gamma_a);
    fprintf('T_s = %1i days \n\n',1/params.gamma_s);
    fprintf('-------------------------------------- \n\n');
    
    
    %% plot first column - no mitigation
    
    f1 = figure(1); set(f1, 'Position', [100 500 800 650]);
    subplot(4,2,1);
    if ismember(counter,[1,2])
        this_p = semilogy(params.t_span, results.total_incidence,'Color',cbf_colors,'LineWidth',2); hold on;
        this_p.Color(4) = 1-0.18*(4-counter);
    else
        this_p = semilogy(params.t_span, results.total_incidence,'--','Color',cbf_colors,'LineWidth',2); hold on;
        this_p.Color(4) = 1-0.18*(4-counter);
    end
    axis([0 params.t_span(end) 10^(-6) 1]);
    %     xlabel('Time (days)');
    ylabel({'Total'; 'Incidence, $i(t)$'},'Interpreter','Latex');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    title('Susceptible Depletion','FontSize',16);
    
    if counter==4
        txt = {'A'};
        text(0.025,0.925,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
        
        old_pos = get(f1, 'Position');
        set(f1,'Position',[old_pos(1),frac_spacing,old_pos(3),frac_scaling]);
        old_pos = get(f1, 'Position');
        
        set(f1,'xticklabel',{[]},'yticklabel',[{' '},{'10^{-4}'},{'10^{-2}'},{'10^0'}]);
        box('off');
        
        txt = {'10^{-6}'};
        text(-0.11,0.05,txt,'Units','normalized',...
            'FontSize',14,'FontWeight','normal','FontName', 'Times');
    end
    
    
    figure(1); subplot(4,2,3);
    if ismember(counter,[1,2])
        this_p=plot(params.t_span, results.proportion_asymp_transmission,'Color',cbf_colors,'LineWidth',2); hold on;
        this_p.Color(4) = 1-0.18*(4-counter);
    else
        this_p=plot(params.t_span, results.proportion_asymp_transmission,'--','Color',cbf_colors,'LineWidth',2); hold on;
        this_p.Color(4) = 1-0.18*(4-counter);
    end
    
    
    axis([0 params.t_span(end) 0 1]);
    %     xlabel('Time (days)');
    ylabel({'Proportion'; 'Asymptomatic'; 'Transmission, $q(t)$'},'Interpreter','Latex');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    if counter==4
        txt = {'B'};
        text(0.025,0.925,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
        
        set(f1,'Position',[old_pos(1), old_pos(2)-frac_scaling, old_pos(3), old_pos(4)])
        old_pos = get(f1, 'Position');
        
        set(f1,'xticklabel',{[]},'yticklabel',[{' '},{'0.5'},{' '}]);
        box('off');
        
        txt = {'1'};
        text(-0.045,0.93,txt,'Units','normalized',...
            'FontSize',14,'FontWeight','normal','FontName', 'Times');
        
        txt = {'0'};
        text(-0.045,0.045,txt,'Units','normalized',...
            'FontSize',14,'FontWeight','normal','FontName', 'Times');
        
    end
    
    
    figure(1); subplot(4,2,5);
    if ismember(counter,[1,2])
        this_p = plot(params.t_span, results.proportion_asymp_incidence,'Color',cbf_colors,'LineWidth',2); hold on;
        this_p.Color(4) = 1-0.18*(4-counter);
    else
        this_p = plot(params.t_span, results.proportion_asymp_incidence,'--','Color',cbf_colors,'LineWidth',2); hold on;
        this_p.Color(4) = 1-0.18*(4-counter);
    end
    axis([0 params.t_span(end) 0 1]);
    %     xlabel('Time (days)');
    ylabel({'Proportion'; 'Asymptomatic'; 'Incidence, $p(t)$'},'Interpreter','Latex');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    if counter==4
        txt = {'C'};
        text(0.025,0.925,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
        
        set(f1,'Position',[old_pos(1), old_pos(2)-frac_scaling, old_pos(3), old_pos(4)])
        old_pos = get(f1, 'Position');
        
        set(f1,'xticklabel',{[]},'yticklabel',[{' '},{'0.5'},{' '}]);
        box('off');
        
        txt = {'1'};
        text(-0.045,0.93,txt,'Units','normalized',...
            'FontSize',14,'FontWeight','normal','FontName', 'Times');
        
        txt = {'0'};
        text(-0.045,0.045,txt,'Units','normalized',...
            'FontSize',14,'FontWeight','normal','FontName', 'Times');
        
    end
    
    
    figure(1); subplot(4,2,7);
    
    if ismember(counter,[1,2])
        semilogy(params.t_span,ones(size(params.t_span)),'k','LineWidth',0.5); hold on;
        this_p=semilogy(params.t_span,results.Rt_fixedpropasymp,'Color',cbf_colors,'LineWidth',2); hold on;
        this_p.Color(4) = 1-0.18*(4-counter);
    else
        this_p=semilogy(params.t_span,results.Rt_assortmixing,'--','Color',cbf_colors,'LineWidth',2); hold on;
        this_p.Color(4) = 1-0.18*(4-counter);
    end
    
    axis([0 params.t_span(end) 0.2 4]);
    yticks([0.25 0.5 1 2 4]);
    xlabel('Time (days)');
    ylabel({'Effective'; 'Reproduction'; 'Number, $\mathcal R_t$'},'Interpreter','Latex');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    
    if counter==4
        txt = {'D'};
        text(0.025,0.925,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
        
        set(f1,'Position',[old_pos(1), old_pos(2)-frac_scaling, old_pos(3), old_pos(4)])
        old_pos = get(f1, 'Position');
        
        set(f1,'yticklabel',[{'0.25'},{'0.5'},{'1'},{'2'},{''},{''}]);
        set(f1,'YminorTick','off');
        box('off');
        
        txt = {'4'};
        text(-0.045,0.95,txt,'Units','normalized',...
            'FontSize',14,'FontWeight','normal','FontName', 'Times');
    end
    
end




%% load no mitigation files

fprintf('With mitigation... \n\n');

% D-F: with mitigation
for counter=1:length(switch_over_var)
    
    this_file = switch_over_var(counter);
    
    
    
    switch this_file
        
       case 1
            % NO assortative mixing
            % load: Ta=5 (Ts=5)
            infile = 'SEIR_fixedpropasymp_twodiseases_sameR0s_110922_T5and5_mit.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [15,32,128]/255; % dark blue
            gamma_s1 = params.gamma_s;
            gamma_a1 = params.gamma_a;
            p1 = params.p;
            
            case 2
            % NO assortative mixing
            % load: Ta=8 (Ts=5)
            infile = 'SEIR_fixedpropasymp_twodiseases_sameR0s_110922_T5and8_mit.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [133,192,249]/255; % light blue
            gamma_s2 = params.gamma_s;
            gamma_a2 = params.gamma_a;
            p2 = params.p;
            
            
        case 3
            % WITH assortative mixing
            % load: Ta=5 (same),
            infile = 'SEIR_assortmixing_twodiseases_sameR0s_120122_fix_T5and5_mit.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [15,32,128]/255; % dark blue
            gamma_s3 = params.gamma_s;
            gamma_a3 = params.gamma_a;
            p_aa3=params.p_aa; p_as3=params.p_as; 
            
        case 4
            % WITH assortative mixing
            % load: Ta=8 (Ta=5)
            infile = 'SEIR_assortmixing_twodiseases_sameR0s_120122_fix_T5and8_mit.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [133,192,249]/255; % light blue
            gamma_s4 = params.gamma_s;
            gamma_a4 = params.gamma_a;
            p_aa4=params.p_aa; p_as4=params.p_as; 
            
    end
    
    
    
    fprintf('Opened file: \n'); % want to be close to 25 days in
    fprintf(strcat(infile,'\n\n'));
    
    fprintf('Infectious periods: \n');
    fprintf('T_a = %1i days \n',1/params.gamma_a);
    fprintf('T_s = %1i days \n\n',1/params.gamma_s);
    fprintf('-------------------------------------- \n\n');
    
    %% plot second column - with mitigation
    
    figure(1);
    subplot(4,2,2);
    if ismember(counter,[1,2])
        h_leg(counter) = semilogy(params.t_span, results.total_incidence,'Color',cbf_colors,'LineWidth',2); hold on;
        h_leg(counter).Color(4) = 1-0.18*(4-counter); % transparency
    else
        h_leg(counter) = semilogy(params.t_span, results.total_incidence,'--','Color',cbf_colors,'LineWidth',2); hold on;
        h_leg(counter).Color(4) = 1-0.18*(4-counter); % transparency
    end
    
    
    axis([0 params.t_span(end) 10^(-6) 1]);
    %     xlabel('Time (days)');
    ylabel({'Total'; 'Incidence, $i(t)$'},'Interpreter','Latex');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    title('Intervention','FontSize',16);
    
    if counter==4
        txt = {'E'};
        text(0.025,0.925,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
        
        old_pos = get(f1, 'Position');
        set(f1,'Position',[old_pos(1),frac_spacing,old_pos(3),frac_scaling]);
        old_pos = get(f1, 'Position');
        
        set(f1,'xticklabel',{[]},'yticklabel',[{' '},{'10^{-4}'},{'10^{-2}'},{'10^0'}]);
        box('off');
        
        txt = {'10^{-6}'};
        text(-0.11,0.05,txt,'Units','normalized',...
            'FontSize',14,'FontWeight','normal','FontName', 'Times');
    end
    
    
    if counter==4
        
        legend_char1 = ['$T_s = ', num2str(1/gamma_s1),'$, $T_a = ', num2str(1/gamma_a1),'$, $p = ', num2str(p1,'%0.2f'), '$'];
        legend_char2 = ['$T_s = ', num2str(1/gamma_s2),'$, $T_a = ', num2str(1/gamma_a2),'$, $p = ', num2str(p2,'%0.2f'), '$'];
        legend_char3 = ['$T_s = T_a = ', num2str(1/gamma_a3),'$, $p_{a|a} = ', num2str(p_aa3,'%0.2f'), ', p_{a|s} = ', num2str(p_as3,'%0.2f'), '$'];
        legend_char4 = ['$T_s = ', num2str(1/gamma_s4),'$, $T_a = ', num2str(1/gamma_a4),'$, $p_{a|a} = ', num2str(p_aa4,'%0.2f'), ', p_{a|s} = ', num2str(p_as4,'%0.2f'), '$'];
%         legend_char3 = ['$T_a = ', num2str(1/gamma_a3),'$'];
%         legend(h_leg,{legend_char1,legend_char2,legend_char3}, 'Interpreter','Latex','Location','NorthEast','FontSize',10.5);
        legend(h_leg,{legend_char1,legend_char2,legend_char3,legend_char4}, 'Interpreter','Latex','Position',[0.70 0.865 0.1 0.05],'FontSize',10);
        
        
        legend boxoff
    end
    
    figure(1); subplot(4,2,4);
    if ismember(counter,[1,2])
        this_p=plot(params.t_span, results.proportion_asymp_transmission,'Color',cbf_colors,'LineWidth',2); hold on;
        this_p.Color(4) = 1-0.18*(4-counter);
    else
        this_p=plot(params.t_span, results.proportion_asymp_transmission,'--','Color',cbf_colors,'LineWidth',2); hold on;
        this_p.Color(4) = 1-0.18*(4-counter);
    end
    %     this_p.Color(4) = 1-0.18*(counter);
    
    axis([0 params.t_span(end) 0 1]);
    %     xlabel('Time (days)');
    ylabel({'Proportion'; 'Asymptomatic'; 'Transmission, $q(t)$'},'Interpreter','Latex');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    if counter==4
        txt = {'F'};
        text(0.025,0.925,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
        
        set(f1,'Position',[old_pos(1), old_pos(2)-frac_scaling, old_pos(3), old_pos(4)])
        old_pos = get(f1, 'Position');
        
        set(f1,'xticklabel',{[]},'yticklabel',[{' '},{'0.5'},{' '}]);
        box('off');
        
        txt = {'1'};
        text(-0.045,0.93,txt,'Units','normalized',...
            'FontSize',14,'FontWeight','normal','FontName', 'Times');
        
        txt = {'0'};
        text(-0.045,0.045,txt,'Units','normalized',...
            'FontSize',14,'FontWeight','normal','FontName', 'Times');
        
    end
    
    
    figure(1); subplot(4,2,6);
    if ismember(counter,[1,2])
        this_p=plot(params.t_span, results.proportion_asymp_incidence,'Color',cbf_colors,'LineWidth',2); hold on;
        this_p.Color(4) = 1-0.18*(4-counter);
    else
        this_p=plot(params.t_span, results.proportion_asymp_incidence,'--','Color',cbf_colors,'LineWidth',2); hold on;
        this_p.Color(4) = 1-0.18*(4-counter);
        
    end
    axis([0 params.t_span(end) 0 1]);
    %     xlabel('Time (days)');
    ylabel({'Proportion'; 'Asymptomatic'; 'Incidence, $p(t)$'},'Interpreter','Latex');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    if counter==4
        txt = {'G'};
        text(0.025,0.925,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
        
        set(f1,'Position',[old_pos(1), old_pos(2)-frac_scaling, old_pos(3), old_pos(4)])
        old_pos = get(f1, 'Position');
        
        set(f1,'xticklabel',{[]},'yticklabel',[{' '},{'0.5'},{' '}]);
        box('off');
        
        txt = {'1'};
        text(-0.045,0.93,txt,'Units','normalized',...
            'FontSize',14,'FontWeight','normal','FontName', 'Times');
        
        txt = {'0'};
        text(-0.045,0.045,txt,'Units','normalized',...
            'FontSize',14,'FontWeight','normal','FontName', 'Times');
        
    end
    
    figure(1); subplot(4,2,8);
    semilogy(params.t_span,ones(size(params.t_span)),'k','LineWidth',0.5); hold on;
    if ismember(counter,[1,2])
        this_p=semilogy(params.t_span,results.Rt_fixedpropasymp,'Color',cbf_colors,'LineWidth',2); hold on;
        this_p.Color(4) = 1-0.18*(4-counter);
    else
        this_p=semilogy(params.t_span,results.Rt_assortmixing,'--','Color',cbf_colors,'LineWidth',2); hold on;
        this_p.Color(4) = 1-0.18*(4-counter);
    end
    
    
    axis([0 params.t_span(end) 0.2 4]);
    
    xlabel('Time (days)');
    ylabel({'Effective'; 'Reproduction'; 'Number, $\mathcal R_t$'},'Interpreter','Latex');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    if counter==4
        txt = {'H'};
        text(0.025,0.925,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
        
        set(f1,'Position',[old_pos(1), old_pos(2)-frac_scaling, old_pos(3), old_pos(4)])
        old_pos = get(f1, 'Position');
        
        yticks([0.25 0.5 1 2 4]);
        set(f1,'yticklabel',[{'0.25'},{'0.5'},{'1'},{'2'},{''},{''}]);
        set(f1,'YminorTick','off');
        box('off');
        
        txt = {'4'};
        text(-0.045,0.95,txt,'Units','normalized',...
            'FontSize',14,'FontWeight','normal','FontName', 'Times');
        
    end
    
end


%% save figure
if save_ans_Fig
    
    folder_location = '../Figures_ms_all/main/';
    saveas(f1,strcat(folder_location,figure_name),'epsc');
    
    fprintf('Figure saved:\n');
    fprintf(strcat(figure_name,'\n\n'));
    
    fprintf('Location:\n');
    fprintf(strcat(folder_location,'\n\n'));
    
else
    
    fprintf('Figure not saved.\n');
    
end


