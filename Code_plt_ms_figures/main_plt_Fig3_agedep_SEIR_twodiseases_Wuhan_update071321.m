
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

figure_name = 'Figure3_agedep_071321';



fprintf('No mitigation... \n\n');

switch_over_var = [1,2,3,4];

%% load no mitigation files
for counter=1:length(switch_over_var)
    
    which_var = switch_over_var(counter);
    
    switch which_var
        
        case 1
            % Ta = Ts = 5 (regular assortivity)
            infile = 'SEIR_agedep_twodiseases_071321_T5and5.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [15,32,128]/255; % dark blue
            gamma_a1 = params.gamma_a;
            gamma_s1 = params.gamma_s;
            beta_a1 = params.beta_a;
            beta_s1 = params.beta_s;
            R0_agedep1 = results.Rt_agedep(1);
            
            
        case 2
            % Ta = Ts = 5 (4x assortivity)
            infile = 'SEIR_agedep_twodiseases_071321_T5and5_ia.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [15,32,128]/255; % dark blue
            gamma_a2 = params.gamma_a;
            gamma_s2 = params.gamma_s;
            beta_a2 = params.beta_a;
            beta_s2 = params.beta_s;
            R0_agedep2 = results.Rt_agedep(1);
            
        case 3
            % Ta = 8, Ts = 5 (regular assortivity)
            infile = 'SEIR_agedep_twodiseases_071321_T5and8.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [133,192,249]/255; % light blue
            gamma_a3 = params.gamma_a;
            gamma_s3 = params.gamma_s;
            beta_a3 = params.beta_a;
            beta_s3 = params.beta_s;
            R0_agedep3 = results.Rt_agedep(1);
            
            
        case 4
            
            % Ta = 8, Ts = 5 (4x assortivity)
            infile = 'SEIR_agedep_twodiseases_071321_T5and8_ia.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [133,192,249]/255; % light blue
            gamma_a4 = params.gamma_a;
            gamma_s4 = params.gamma_s;
            beta_a4 = params.beta_a;
            beta_s4 = params.beta_s;
            R0_agedep4 = results.Rt_agedep(1);
            
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
    
    ind = find(results.I_all_total_fraction>10^-6);
    fprintf('Time of first appearance: \n'); % want to be close to 25 days in
    fprintf('%2.2f days \n\n',params.t_span(ind(1)));
    
    fprintf('-------------------------------------- \n\n');
    
    %% plot panels A-E: no mitigation
    f1 = figure(1); set(f1, 'Position', [100 500 800 750]);
    subplot(5,2,1);
    if params.increased_assortativity_yesno==0
        semilogy(params.t_span, results.I_all_total_fraction,'Color',cbf_colors,'LineWidth',1); hold on;
    else
        this_p = semilogy(params.t_span, results.I_all_total_fraction,'--','Color',cbf_colors,'LineWidth',2); hold on;
        this_p.Color(4) = 1-0.18*2;
    end
    axis([0 params.t_span(end) 10^(-6) 1]);
    xlabel('Time (days)'); ylabel({'Fraction'; 'Infections'});
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    if counter==1
        title('Susceptible Depletion','FontSize',16);
        txt = {'A'};
        text(0.025,1.1,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
    
    
    figure(1); subplot(5,2,3);
    if params.increased_assortativity_yesno==0
        plot(params.t_span, results.proportion_asymptomatic_incidence,'Color',cbf_colors,'LineWidth',1); hold on;
    else
        this_p=plot(params.t_span, results.proportion_asymptomatic_incidence,'--','Color',cbf_colors,'LineWidth',2); hold on;
        this_p.Color(4) = 1-0.18*2;
    end
    axis([0 params.t_span(end) 0 1]);
    xlabel('Time (days)'); ylabel({'Proportion'; 'Asymptomatic'; 'Incidence'});
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    if ismember(counter,[1,3])
        txt = {'B'};
        text(0.025,1.1,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
    
    
    figure(1); subplot(5,2,5);
    if params.increased_assortativity_yesno==0
        plot(params.t_span, results.proportion_asymptomatic_transmission,'Color',cbf_colors,'LineWidth',1); hold on;
    else
        this_p=plot(params.t_span, results.proportion_asymptomatic_transmission,'--','Color',cbf_colors,'LineWidth',2); hold on;
        this_p.Color(4) = 1-0.18*2;
    end
    axis([0 params.t_span(end) 0 1]);
    xlabel('Time (days)'); ylabel({'Proportion'; 'Asymptomatic'; 'Transmission'});
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    if counter==1
        txt = {'C'};
        text(0.025,1.1,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
    
    figure(1); subplot(5,2,7);
    if params.increased_assortativity_yesno==0
        plot(params.t_span, results.ave_age_i_all,'Color',cbf_colors,'LineWidth',1); hold on;
    else
        this_p = plot(params.t_span, results.ave_age_i_all,'--','Color',cbf_colors,'LineWidth',2); hold on;
        this_p.Color(4) = 1-0.18*2;
    end
    axis([0 params.t_span(end) 20 50]);
    xlabel('Time (days)'); ylabel({'Average Age'; 'of Infection'});
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    if counter==1
        txt = {'D'};
        text(0.025,1.1,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
    
    
    figure(1); subplot(5,2,9);
    if counter==1
        semilogy(params.t_span,ones(size(params.t_span)),'k','LineWidth',0.5); hold on;
    end
    
    if params.increased_assortativity_yesno==0
        semilogy(params.t_span,results.Rt_agedep,'Color',cbf_colors,'LineWidth',1); hold on;
    else
        this_p=semilogy(params.t_span,results.Rt_agedep,'--','Color',cbf_colors,'LineWidth',2); hold on;
        this_p.Color(4) = 1-0.18*2;
        
    end
    
    if counter==4
        semilogy(params.t_span(1:20:end),ones(size(params.t_span(1:20:end))),'k','LineWidth',0.5); hold on;
    end
    axis([0 params.t_span(end) 10^-1 10^1]);
    xlabel('Time (days)'); ylabel({'Effective'; 'Reproduction'; 'Number'});
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    if counter==1
        txt = {'E'};
        text(0.025,1.1,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
end



%%
fprintf('With mitigation... \n\n');

%% load with mitigation files
for counter=1:length(switch_over_var)
    
    which_var = switch_over_var(counter);
    
    switch which_var
        
        case 1
            % Ta = Ts = 5 (regular assortivity)
            infile = 'SEIR_agedep_twodiseases_071321_T5and5_mit.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [15,32,128]/255; % dark blue
            gamma_a1 = params.gamma_a;
            gamma_s1 = params.gamma_s;
            beta_a1 = params.beta_a;
            beta_s1 = params.beta_s;
            R0_agedep1 = results.Rt_agedep(1);
            
            
        case 2
            % Ta = Ts = 5 (4x assortivity)
            infile = 'SEIR_agedep_twodiseases_071321_T5and5_mit_ia.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [15,32,128]/255; % dark blue
            gamma_a2 = params.gamma_a;
            gamma_s2 = params.gamma_s;
            beta_a2 = params.beta_a;
            beta_s2 = params.beta_s;
            R0_agedep2 = results.Rt_agedep(1);
            
        case 3
            % Ta = 8, Ts = 5 (regular assortivity)
            infile = 'SEIR_agedep_twodiseases_071321_T5and8_mit.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [133,192,249]/255; % light blue
            gamma_a3 = params.gamma_a;
            gamma_s3 = params.gamma_s;
            beta_a3 = params.beta_a;
            beta_s3 = params.beta_s;
            R0_agedep3 = results.Rt_agedep(1);
            
            
        case 4
            
            % Ta = 8, Ts = 5 (4x assortivity)
            infile = 'SEIR_agedep_twodiseases_071321_T5and8_mit_ia.mat';
            load(strcat('./sim_data/',infile));
            cbf_colors = [133,192,249]/255; % light blue
            gamma_a4 = params.gamma_a;
            gamma_s4 = params.gamma_s;
            beta_a4 = params.beta_a;
            beta_s4 = params.beta_s;
            R0_agedep4 = results.Rt_agedep(1);
            
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
    
    ind = find(results.I_all_total_fraction>10^-6);
    fprintf('Time of first appearance: \n'); % want to be close to 25 days in
    fprintf('%2.2f days \n\n',params.t_span(ind(1)));
    
    fprintf('-------------------------------------- \n\n');
    
    %% plot panels E-I: with mitigation
    figure(1);
    subplot(5,2,2);
    if params.increased_assortativity_yesno==0
        this_h(counter)=semilogy(params.t_span, results.I_all_total_fraction,'Color',cbf_colors,'LineWidth',1); hold on;
    else
        this_h(counter)=semilogy(params.t_span, results.I_all_total_fraction,'--','Color',cbf_colors,'LineWidth',2); hold on;
        this_h(counter).Color(4) = 1-0.18*2;
    end
    axis([0 params.t_span(end) 10^(-6) 1]);
    xlabel('Time (days)'); ylabel({'Fraction'; 'Infections'});
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    
    
    if counter==1
        
        title('Intervention','FontSize',16);
        txt = {'F'};
        text(0.025,1.1,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
    if counter==4
        
        legend_char1 = ['$T_a = ', num2str(1/gamma_a1),'$'];
        legend_char2 = ['$T_a = ', num2str(1/gamma_a2),'$ (4x)'];
        legend_char3 = ['$T_a = ', num2str(1/gamma_a3),'$'];
        legend_char4 = ['$T_a = ', num2str(1/gamma_a4),'$ (4x)'];
        legend(this_h,{legend_char1,legend_char2,legend_char3,legend_char4}, 'Interpreter','Latex','Location','NorthEast','FontSize',11);
        legend boxoff
    end

    
    
    figure(1); subplot(5,2,4);
    if params.increased_assortativity_yesno==0
        plot(params.t_span, results.proportion_asymptomatic_incidence,'Color',cbf_colors,'LineWidth',1); hold on;
    else
        this_p = plot(params.t_span, results.proportion_asymptomatic_incidence,'--','Color',cbf_colors,'LineWidth',2); hold on;
        this_p.Color(4) = 1-0.18*2;
    end
    axis([0 params.t_span(end) 0 1]);
    xlabel('Time (days)'); ylabel({'Proportion'; 'Asymptomatic'; 'Incidence'});
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    if counter==1
        txt = {'G'};
        text(0.025,1.1,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
    
    
    
    figure(1); subplot(5,2,6);
    if params.increased_assortativity_yesno==0
        plot(params.t_span, results.proportion_asymptomatic_transmission,'Color',cbf_colors,'LineWidth',1); hold on;
    else
        this_p=plot(params.t_span, results.proportion_asymptomatic_transmission,'--','Color',cbf_colors,'LineWidth',2); hold on;
        this_p.Color(4) = 1-0.18*2;
    end
    axis([0 params.t_span(end) 0 1]);
    xlabel('Time (days)'); ylabel({'Proportion'; 'Asymptomatic'; 'Transmission'});
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    if counter==1
        txt = {'H'};
        text(0.025,1.1,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
    figure(1); subplot(5,2,8);
    if params.increased_assortativity_yesno==0
        plot(params.t_span, results.ave_age_i_all,'Color',cbf_colors,'LineWidth',1); hold on;
    else
        this_p = plot(params.t_span, results.ave_age_i_all,'--','Color',cbf_colors,'LineWidth',2); hold on;
        this_p.Color(4) = 1-0.18*2;
    end
    axis([0 params.t_span(end) 40 45]);
    xlabel('Time (days)'); ylabel({'Average Age'; 'of Infection'});
    yticks([40:45]); yticklabels({'40','','','','','45'});
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    if counter==1
        txt = {'I'};
        text(0.025,1.1,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
    figure(1); subplot(5,2,10);
    
    if counter==1
        semilogy(params.t_span(1:20:end),ones(size(params.t_span(1:20:end))),'k','LineWidth',0.5); hold on;
    end
    
    if params.increased_assortativity_yesno==0
        semilogy(params.t_span,results.Rt_agedep,'Color',cbf_colors,'LineWidth',1); hold on;
    else
        this_p = semilogy(params.t_span,results.Rt_agedep,'--','Color',cbf_colors,'LineWidth',2); hold on;
        this_p.Color(4) = 1-0.18*2;
    end
    
    axis([0 params.t_span(end) 10^-1 10^1]);
    xlabel('Time (days)'); ylabel({'Effective'; 'Reproduction'; 'Number'});
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    
    
    if counter==1
        txt = {'J'};
        text(0.025,1.1,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
end



%% save figure
if save_ans
    
    folder_location = './../Figures_ms_all/main/';
    saveas(f1,strcat(folder_location,figure_name),'epsc');
    
    fprintf('Figure saved:\n');
    fprintf(strcat(figure_name,'\n\n'));
    
    fprintf('Location:\n');
    fprintf(strcat(folder_location,'\n\n'));
    
end
