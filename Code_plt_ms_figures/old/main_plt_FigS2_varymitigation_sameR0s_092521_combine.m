
%% Figure 1
% A-D: no mit, increase time-scale differences
% E-H: with mit, increase time-scale differences

clear all; close all; clc;

%% save figure 1?
save_ans_Fig = 0;
% 0: don't save
% 1: save

figure_name = 'FigS2_varymitigation_sameR0s_092521_combine';

f1 = figure(1); set(f1, 'Position', [100 500 800 650]);

%% A-D: vary strength of mitigation
for counter=1:3
    
    switch counter
        
        case 1
            % load:
            infile = 'SEIR_fixedpropasymp_twodiseases_sameR0s_072221_T5and8_mit1.mat';
            load(strcat('./sim_data/',infile));
            %             cbf_colors_vector = [15,32,128]/255; % dark blue
            %             gamma_a1 = params.gamma_a;
            
        case 2
            % load:
            infile = 'SEIR_fixedpropasymp_twodiseases_sameR0s_072221_T5and8_mit2.mat';
            load(strcat('./sim_data/',infile));
            %             gamma_a2 = params.gamma_a;
            
        case 3
            % load:
            infile = 'SEIR_fixedpropasymp_twodiseases_sameR0s_072221_T5and8_mit3.mat';
            load(strcat('./sim_data/',infile));
            %             gamma_a3 = params.gamma_a;
            
    end
    
    fprintf('Opened file: \n'); % want to be close to 25 days in
    fprintf(strcat(infile,'\n\n'));
    
    
    %% plot
    
    subplot(4,2,1);
    this_q(counter) = semilogy(params.t_span, results.I_tot,'Color',cbf_colors_vector(counter,:),'LineWidth',2); hold on;
    this_q(counter).Color(4) = 1-0.18*(counter);
    axis([0 params.t_span(end) 10^(-6) 1]);
    xlabel('Time (days)');
    ylabel({'Fraction'; 'Infections'});
    title('Fixing the Decay Rate of Infections');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 12;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    if counter==1
        txt = {'A'};
        text(0.025,1.075,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
    
    figure(1); subplot(4,2,3);
    this_p = plot(params.t_span, results.proportion_asymp_incidence,'Color',cbf_colors_vector(counter,:),'LineWidth',2); hold on;
    this_p.Color(4) = 1-0.18*(counter);
    
    axis([0 params.t_span(end) 0 1]);
    xlabel('Time (days)');
    ylabel({'Proportion'; 'Asymptomatic'; 'Incidence'});
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    if counter==1
        txt = {'B'};
        text(0.025,1.075,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
    
    figure(1); subplot(4,2,5);
    this_p=plot(params.t_span, results.proportion_asymp_transmission,'Color',cbf_colors_vector(counter,:),'LineWidth',2); hold on;
    this_p.Color(4) = 1-0.18*(counter);
    
    axis([0 params.t_span(end) 0 1]);
    xlabel('Time (days)');
    ylabel({'Proportion'; 'Asymptomatic'; 'Transmission'});
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    if counter==1
        txt = {'C'};
        text(0.025,1.075,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
    figure(1); subplot(4,2,7);
    %     plot(params.t_span,Rt_fixedpropasymp,'Color',cbf_colors_vector,'LineWidth',2); hold on;
    semilogy(params.t_span,ones(size(params.t_span)),'k','LineWidth',0.5); hold on;
    this_p=semilogy(params.t_span,results.Rt_fixedpropasymp,'Color',cbf_colors_vector(counter,:),'LineWidth',2); hold on;
    R_inf(counter)=Rt_fixedpropasymp(end);
    this_p.Color(4) = 1-0.18*(counter);
    
    axis([0 params.t_span(end) 0.1 10]);
    xlabel('Time (days)'); ylabel({'Effective'; 'Reproduction'; 'Number'});
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    %     this_p.Color(4) = 1-0.18*(counter);
    
    
    if counter==1
        txt = {'D'};
        text(0.025,1.075,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
end

if counter==3
    
    figure(1); subplot(4,2,1);
    legend_char1 = ['$R_{\infty} = ', num2str(R_inf(1),'%2.2f'),'$'];
    legend_char2 = ['$R_{\infty} = ', num2str(R_inf(2),'%2.2f'),'$'];
    legend_char3 = ['$R_{\infty} = ', num2str(R_inf(3),'%2.2f'),'$'];
    legend(this_q,{legend_char1,legend_char2,legend_char3}, 'Interpreter','Latex','Location','NorthEast','FontSize',12);
    legend boxoff
end


%% E-H: vary speed of mitigation
for counter=1:3
    
    switch counter
        
        case 1
            % load:
            infile = 'SEIR_fixedpropasymp_twodiseases_sameR0s_071321_T5and8_mit1.mat';
            load(strcat('./sim_data/',infile));
            %             cbf_colors_vector = [15,32,128]/255; % dark blue
            %             gamma_a1 = params.gamma_a;
            
        case 2
            % load:
            infile = 'SEIR_fixedpropasymp_twodiseases_sameR0s_071321_T5and8_mit2.mat';
            load(strcat('./sim_data/',infile));
            %             gamma_a2 = params.gamma_a;
            
        case 3
            % load:
            infile = 'SEIR_fixedpropasymp_twodiseases_sameR0s_071321_T5and8_mit3.mat';
            load(strcat('./sim_data/',infile));
            %             gamma_a3 = params.gamma_a;
            
    end
    
    fprintf('Opened file: \n'); % want to be close to 25 days in
    fprintf(strcat(infile,'\n\n'));
    
    
    %% plot
    
    subplot(4,2,2);
    this_q(counter) = semilogy(params.t_span, results.I_tot,'Color',cbf_colors_vector(counter,:),'LineWidth',2); hold on;
    this_q(counter).Color(4) = 1-0.18*(counter);
    axis([0 params.t_span(end) 10^(-6) 1]);
    xlabel('Time (days)');
    ylabel({'Fraction'; 'Infections'});
    title('Varying the Decay Rate of Infections');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 12;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    if counter==1
        txt = {'E'};
        text(0.025,1.075,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
    
    subplot(4,2,4);
    this_p = plot(params.t_span, results.proportion_asymp_incidence,'Color',cbf_colors_vector(counter,:),'LineWidth',2); hold on;
    this_p.Color(4) = 1-0.18*(counter);
    
    axis([0 params.t_span(end) 0 1]);
    xlabel('Time (days)');
    ylabel({'Proportion'; 'Asymptomatic'; 'Incidence'});
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    if counter==1
        txt = {'F'};
        text(0.025,1.075,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
    
    figure(1); subplot(4,2,6);
    this_p=plot(params.t_span, results.proportion_asymp_transmission,'Color',cbf_colors_vector(counter,:),'LineWidth',2); hold on;
    this_p.Color(4) = 1-0.18*(counter);
    
    axis([0 params.t_span(end) 0 1]);
    xlabel('Time (days)');
    ylabel({'Proportion'; 'Asymptomatic'; 'Transmission'});
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    if counter==1
        txt = {'G'};
        text(0.025,1.075,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
    figure(1); subplot(4,2,8);
    %     plot(params.t_span,Rt_fixedpropasymp,'Color',cbf_colors_vector,'LineWidth',2); hold on;
    semilogy(params.t_span,ones(size(params.t_span)),'k','LineWidth',0.5); hold on;
    this_p=semilogy(params.t_span,results.Rt_fixedpropasymp,'Color',cbf_colors_vector(counter,:),'LineWidth',2); hold on;
    R_inf(counter)=Rt_fixedpropasymp(end);
    this_p.Color(4) = 1-0.18*(counter);
    
    axis([0 params.t_span(end) 0.1 10]);
    xlabel('Time (days)'); ylabel({'Effective'; 'Reproduction'; 'Number'});
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    %     this_p.Color(4) = 1-0.18*(counter);
    
    
    if counter==1
        txt = {'H'};
        text(0.025,1.075,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
end

if counter==3
    
    figure(1); subplot(4,2,2);
    legend_char1 = ['$R_{\infty} = ', num2str(R_inf(1),'%2.2f'),'$'];
    legend_char2 = ['$R_{\infty} = ', num2str(R_inf(2),'%2.2f'),'$'];
    legend_char3 = ['$R_{\infty} = ', num2str(R_inf(3),'%2.2f'),'$'];
    legend(this_q,{legend_char1,legend_char2,legend_char3}, 'Interpreter','Latex','Location','NorthEast','FontSize',12);
    legend boxoff
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


