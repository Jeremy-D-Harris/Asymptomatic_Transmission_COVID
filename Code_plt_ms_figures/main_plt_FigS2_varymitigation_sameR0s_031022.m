
%% Figure 1
% A-D: no mit, increase time-scale differences
% E-H: with mit, increase time-scale differences

clear all; close all; clc;

%% save figure 1?
save_ans_Fig = 0;
% 0: don't save
% 1: save

frac_spacing = 0.74;
frac_scaling = 0.2;

figure_name = 'FigureS2_varymitigation_sameR0s_031022';

f1 = figure(1); set(f1, 'Position', [100 500 800 650]);

%% A-D: vary time to mitigation

for counter=1:3
    
    switch counter
        
        case 1
            
            
            % load:
            infile = 'SEIR_fixedpropasymp_twodiseases_sameR0s_012422_T5and8_mit2_20days.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [0 0 0]; % black
            R_inf(1)=results.Rt_fixedpropasymp(end);
            t_min(1)=results.t_min;
            %             gamma_a1 = params.gamma_a;
            
            
        case 2
            % load:
            infile = 'SEIR_fixedpropasymp_twodiseases_sameR0s_012422_T5and8_mit2.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [133,192,249]/255; % light blue
            R_inf(2)=results.Rt_fixedpropasymp(end);
            t_min(2)=results.t_min;
            %             gamma_a2 = params.gamma_a;
            
        case 3
            
            % load:
            infile = 'SEIR_fixedpropasymp_twodiseases_sameR0s_012422_T5and8_mit2_40days.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [0.5 0.5 0.5]; % gray
            R_inf(3)=results.Rt_fixedpropasymp(end);
            t_min(3)=results.t_min;
            %             gamma_a3 = params.gamma_a;
            
            
            
            
            
    end
    
    fprintf('Opened file: \n'); % want to be close to 25 days in
    fprintf(strcat(infile,'\n\n'));
    
    
    %% plot first column - no mitigation
    
    f1 = figure(1); set(f1, 'Position', [100 500 800 650]);
    subplot(4,2,1);
    %     this_p = semilogy(params.t_span, results.total_incidence,'Color',cbf_colors,'LineWidth',2); hold on;
    this_q(counter) = semilogy(params.t_span, results.total_incidence,'Color',cbf_colors,'LineWidth',2); hold on;
    this_q(counter).Color(4) = 1-0.18*(counter-1); % transparency
    %     this_q(counter).Color(4) = 0.64+0.18*(counter-1);
    axis([0 params.t_span(end) 10^(-6) 1]);
    
    %     xlabel('Time (days)');
    ylabel({'Total'; 'Incidence, $i(t)$'},'Interpreter','Latex');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    title('Varying the Mitigation Onset Time','FontSize',16);
    
    if counter==3
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
    
    if counter==3
        
        legend_char1 = ['$T_{m} = ', num2str(t_min(1),'%2d'),'$ days'];
        legend_char2 = ['$T_{m} = ', num2str(t_min(2),'%2d'),'$ days'];
        legend_char3 = ['$T_{m} = ', num2str(t_min(3),'%2d'),'$ days'];
        %     legend(this_h,{'1','2','3'}, 'Interpreter','Latex');
        legend(this_q,{legend_char1,legend_char2,legend_char3}, 'Interpreter','Latex','Location','NorthEast','FontSize',12);
        legend boxoff
        
    end
    
    figure(1);
    subplot(4,2,3);
    this_p=plot(params.t_span, results.proportion_asymp_transmission,'Color',cbf_colors,'LineWidth',2); hold on;
    this_p.Color(4) = 1-0.18*(counter-1); % transparency
    %     this_p.Color(4) = 0.64+0.18*(counter-1);
    
    axis([0 params.t_span(end) 0 1]);
    %     xlabel('Time (days)');
    ylabel({'Proportion'; 'Asymptomatic'; 'Transmission, $q(t)$'},'Interpreter','Latex');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    if counter==3
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
    this_p = plot(params.t_span, results.proportion_asymp_incidence,'Color',cbf_colors,'LineWidth',2); hold on;
    this_p.Color(4) = 1-0.18*(counter-1); % transparency
    %     this_p.Color(4) = 0.64+0.18*(counter-1);
    
    axis([0 params.t_span(end) 0 1]);
    %     xlabel('Time (days)');
    ylabel({'Proportion'; 'Asymptomatic'; 'Incidence, $p(t)$'},'Interpreter','Latex');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    if counter==3
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
    %     plot(params.t_span,Rt_fixedpropasymp,'Color',cbf_colors,'LineWidth',2); hold on;
    semilogy(params.t_span,ones(size(params.t_span)),'k','LineWidth',0.5); hold on;
    this_p=semilogy(params.t_span,results.Rt_fixedpropasymp,'Color',cbf_colors,'LineWidth',2); hold on;
    R_inf(counter) = results.Rt_fixedpropasymp(end);
    this_p.Color(4) = 1-0.18*(counter-1); % transparency
    %     this_p.Color(4) = 0.64+0.18*(counter-1);
    
    axis([0 params.t_span(end) 0.1 4]);
    yticks([0.125 0.25 0.5 1 2 4]);
    xlabel('Time (days)'); ylabel({'Effective'; 'Reproduction'; 'Number, $\mathcal R_t$'},'Interpreter','Latex');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    %     this_p.Color(4) = 1-0.18*(counter);
    
    if counter==3
        
        figure(1); subplot(4,2,7);
        
        txt = {'D'};
        text(0.025,0.95,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
        
        
        set(f1,'Position',[old_pos(1), old_pos(2)-frac_scaling, old_pos(3), old_pos(4)])
        old_pos = get(f1, 'Position');
        
        set(f1,'yticklabel',[{'0.125'},{'0.25'},{'0.5'},{'1'},{'2'},{''},{''}]);
        set(f1,'YminorTick','off');
        box('off');
        
        txt = {'4'};
        text(-0.045,0.95,txt,'Units','normalized',...
            'FontSize',14,'FontWeight','normal','FontName', 'Times');
        
    end
    
end




%% E-H: vary decay rates
for counter=1:3
    
    switch counter
        
        case 1
            
            % load:
            infile = 'SEIR_fixedpropasymp_twodiseases_sameR0s_012422_T5and8_mit3.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [0 0 0]; % black
            R_inf(3)=results.Rt_fixedpropasymp(end);
            %             gamma_a3 = params.gamma_a;
            
            
            
        case 2
            % load:
            infile = 'SEIR_fixedpropasymp_twodiseases_sameR0s_012422_T5and8_mit2.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [133,192,249]/255; % light blue
            R_inf(2)=results.Rt_fixedpropasymp(end);
            %             gamma_a2 = params.gamma_a;
            
        case 3
            
            % load:
            infile = 'SEIR_fixedpropasymp_twodiseases_sameR0s_012422_T5and8_mit1.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [0.5 0.5 0.5]; % gray
            R_inf(1)=results.Rt_fixedpropasymp(end);
            %             gamma_a1 = params.gamma_a;
            
            
    end
    
    fprintf('Opened file: \n'); % want to be close to 25 days in
    fprintf(strcat(infile,'\n\n'));
    
    
    %% plot second column - with mitigation
    
    figure(1);
    subplot(4,2,2);
    this_h(counter) = semilogy(params.t_span, results.total_incidence,'Color',cbf_colors,'LineWidth',2); hold on;
    this_h(counter).Color(4) = 1-0.18*(counter-1); % transparency
    %     this_h(counter).Color(4) = 0.64+0.18*(counter-1); % transparency
    
    axis([0 params.t_span(end) 10^(-6) 1]);
    %     xlabel('Time (days)');
    ylabel({'Total'; 'Incidence, $i(t)$'},'Interpreter','Latex');
    
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    title('Varying the Mitigation Intensity','FontSize',16);
    
    if counter==3
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
    
    if counter==3
        
        %     figure(1); subplot(4,2,2);
        legend_char1 = ['$R_{\infty} = ', num2str(R_inf(3),'%2.2f'),'$'];
        legend_char2 = ['$R_{\infty} = ', num2str(R_inf(2),'%2.2f'),'$'];
        legend_char3 = ['$R_{\infty} = ', num2str(R_inf(1),'%2.2f'),'$'];
        legend(this_h,{legend_char1,legend_char2,legend_char3}, 'Interpreter','Latex','Location','NorthEast','FontSize',12);
        legend boxoff
    end
    
    figure(1); subplot(4,2,4);
    this_p=plot(params.t_span, results.proportion_asymp_transmission,'Color',cbf_colors,'LineWidth',2); hold on;
    this_p.Color(4) = 1-0.18*(counter-1); % transparency
%     this_p.Color(4) = 0.64+0.18*(counter-1);
    
    axis([0 params.t_span(end) 0 1]);
    %     xlabel('Time (days)');
    ylabel({'Proportion'; 'Asymptomatic'; 'Transmission, $q(t)$'},'Interpreter','Latex');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    
    if counter==3
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
    this_p=plot(params.t_span, results.proportion_asymp_incidence,'Color',cbf_colors,'LineWidth',2); hold on;
    this_p.Color(4) = 1-0.18*(counter-1); % transparency
%     this_p.Color(4) = 0.64+0.18*(counter-1);
    axis([0 params.t_span(end) 0 1]);
    %     xlabel('Time (days)');
    ylabel({'Proportion'; 'Asymptomatic'; 'Incidence, $p(t)$'},'Interpreter','Latex');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    if counter==3
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
    this_p=semilogy(params.t_span,results.Rt_fixedpropasymp,'Color',cbf_colors,'LineWidth',2); hold on;
    R_inf(counter) = results.Rt_fixedpropasymp(end);
    this_p.Color(4) = 1-0.18*(counter-1); % transparency
%     this_p.Color(4) = 0.64+0.18*(counter-1);
    
    axis([0 params.t_span(end) 0.1 4]);
    yticks([0.125 0.25 0.5 1 2 4]);
    xlabel('Time (days)'); ylabel({'Effective'; 'Reproduction'; 'Number, $\mathcal R_t$'},'Interpreter','Latex');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    if counter==3
        txt = {'H'};
        text(0.025,0.95,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
        
        set(f1,'Position',[old_pos(1), old_pos(2)-frac_scaling, old_pos(3), old_pos(4)])
        old_pos = get(f1, 'Position');
        
        set(f1,'yticklabel',[{'0.125'},{'0.25'},{'0.5'},{'1'},{'2'},{''},{''}]);
        set(f1,'YminorTick','off');
        box('off');
        
        txt = {'4'};
        text(-0.045,0.95,txt,'Units','normalized',...
            'FontSize',14,'FontWeight','normal','FontName', 'Times');
    end
    
end




%% save figure
if save_ans_Fig
    
    folder_location = './../Figures_ms_all/supp/';
    saveas(f1,strcat(folder_location,figure_name),'epsc');
    
    fprintf('Figure saved:\n'); % want to be close to 25 days in
    fprintf(strcat(figure_name,'\n\n'));
    
    fprintf('Location:\n'); % want to be close to 25 days in
    fprintf(strcat(folder_location,'\n\n'));
    
end


