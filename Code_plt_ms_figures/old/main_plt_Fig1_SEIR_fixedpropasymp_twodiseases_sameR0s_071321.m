
%% Figure 1 - parametrised with R_0,a = R_0,s
% A-D: no mit, increase time-scale differences
% E-H: with mit, increase time-scale differences

clear all; close all; clc;

%% save figure?
save_ans_Fig = 0;
% 0: don't save
% 1: save

figure_name = 'Figure1_fixedpropasymp_sameR0s_071321';

%% load no mitigation files

switch_over_Ta = [5,6,8];

fprintf('No mitigation... \n\n');

% A-D: no mitigation
for counter=1:length(switch_over_Ta)
    
    which_Ta = switch_over_Ta(counter);
    
    switch which_Ta
        
        case 5
            % load: Ta=5 (same)
            infile = 'SEIR_fixedpropasymp_twodiseases_sameR0s_071321_T5and5.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [15,32,128]/255; % dark blue
            gamma_a1 = params.gamma_a;
            
        case 6
            % load: Ta=6
            infile = 'SEIR_fixedpropasymp_twodiseases_sameR0s_071321_T5and6.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [169,90,161]/255; % violet
            gamma_a2 = params.gamma_a;
            
        case 8
            % load: Ta=8
            infile = 'SEIR_fixedpropasymp_twodiseases_sameR0s_071321_T5and8.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [133,192,249]/255; % light blue
            gamma_a3 = params.gamma_a;
            
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
    this_p = semilogy(params.t_span, results.I_tot,'Color',cbf_colors,'LineWidth',2); hold on;
    this_p.Color(4) = 1-0.18*(counter);
    axis([0 params.t_span(end) 10^(-6) 1]);
    xlabel('Time (days)');
    ylabel({'Fraction'; 'Infections'});
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    title('Susceptible Depletion','FontSize',16);
    
    if counter==1
        txt = {'A'};
        text(0.025,1.075,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
    
    figure(1); subplot(4,2,3);
    this_p=plot(params.t_span, results.proportion_asymp_transmission,'Color',cbf_colors,'LineWidth',2); hold on;
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
        txt = {'B'};
        text(0.025,1.075,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
    
    figure(1); subplot(4,2,5);
    this_p = plot(params.t_span, results.proportion_asymp_incidence,'Color',cbf_colors,'LineWidth',2); hold on;
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
        txt = {'C'};
        text(0.025,1.075,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
    figure(1); subplot(4,2,7);
    %     plot(params.t_span,Rt_fixedpropasymp,'Color',cbf_colors,'LineWidth',2); hold on;
    semilogy(params.t_span,ones(size(params.t_span)),'k','LineWidth',0.5); hold on;
    this_p=semilogy(params.t_span,results.Rt_fixedpropasymp,'Color',cbf_colors,'LineWidth',2); hold on;
    
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



%% load no mitigation files

fprintf('With mitigation... \n\n');

% D-F: with mitigation
for counter=1:length(switch_over_Ta)
    
    which_Ta = switch_over_Ta(counter);
    
    switch which_Ta
        
        case 5
            % load: Ta=5 (same)
            infile = 'SEIR_fixedpropasymp_twodiseases_sameR0s_071321_T5and5_mit.mat';
            load(strcat('./sim_data/',infile));
            %             cbf_colors = [0.5,0.5,0.5]; % gray
            cbf_colors = [15,32,128]/255; % dark blue
            gamma_a1 = params.gamma_a;
            
        case 6
            % load: Ta=6
            infile = 'SEIR_fixedpropasymp_twodiseases_sameR0s_071321_T5and6_mit.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [169,90,161]/255; % violet
            gamma_a2 = params.gamma_a;
            
            
        case 8
            % load: Ta=8
            infile = 'SEIR_fixedpropasymp_twodiseases_sameR0s_071321_T5and8_mit.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [133,192,249]/255; % light blue
            gamma_a3 = params.gamma_a;
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
    h(counter) = semilogy(params.t_span, results.I_tot,'Color',cbf_colors,'LineWidth',2); hold on;
    h(counter).Color(4) = 1-0.18*(counter); % transparency
    
    axis([0 params.t_span(end) 10^(-6) 1]);
    xlabel('Time (days)');
    ylabel({'Fraction'; 'Infections'});
    
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    title('Intervention','FontSize',16);
    
    if counter==1
        txt = {'E'};
        %     text(5,9.4,txt,
        text(0.025,1.075,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
    if counter==3
        
        legend_char1 = ['$T_a = ', num2str(1/gamma_a1),'$'];
        legend_char2 = ['$T_a = ', num2str(1/gamma_a2),'$'];
        legend_char3 = ['$T_a = ', num2str(1/gamma_a3),'$'];
        legend(h,{legend_char1,legend_char2,legend_char3}, 'Interpreter','Latex','Location','NorthEast','FontSize',12);
        
        legend boxoff
    end
    
   
    figure(1); subplot(4,2,4);
    this_p=plot(params.t_span, results.proportion_asymp_transmission,'Color',cbf_colors,'LineWidth',2); hold on;
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
        txt = {'F'};
        text(0.025,1.075,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
    
    figure(1); subplot(4,2,6);
    this_p=plot(params.t_span, results.proportion_asymp_incidence,'Color',cbf_colors,'LineWidth',2); hold on;
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
        txt = {'G'};
        text(0.025,1.075,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
    figure(1); subplot(4,2,8);
    semilogy(params.t_span,ones(size(params.t_span)),'k','LineWidth',0.5); hold on;
    this_p=semilogy(params.t_span,results.Rt_fixedpropasymp,'Color',cbf_colors,'LineWidth',2); hold on;
    this_p.Color(4) = 1-0.18*(counter);
    
    axis([0 params.t_span(end) 0.1 10]);
    xlabel('Time (days)');
    ylabel({'Effective'; 'Reproduction'; 'Number'});
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    if counter==1
        txt = {'H'};
        text(0.025,1.075,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
end


%% save figure
if save_ans_Fig
    
    folder_location = './../Figures_ms_all/main/';
    saveas(f1,strcat(folder_location,figure_name),'epsc');
    
    fprintf('Figure saved:\n'); 
    fprintf(strcat(figure_name,'\n\n'));
    
    fprintf('Location:\n'); 
    fprintf(strcat(folder_location,'\n\n'));
    
    else
    
    fprintf('Figure not saved.\n');
    
end


