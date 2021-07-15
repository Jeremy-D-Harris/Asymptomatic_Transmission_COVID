
%% Figure S13 - plot p(a|a) and p(a|s) versus time
% see methods

clear all; close all; clc;

%% save Figure ?
save_ans_Fig = 0;
% 0: don't save
% 1: save

figure_name = 'FigS13_agedep_paa_pas_071321';

switch_over_var = [1,2,3,4];

f1 = figure(1); set(f1,'Position',[240   300   1120   420]);

if 1
    % A: no mitigation
    for counter=1:length(switch_over_var)
        
        which_var = switch_over_var(counter);
        
        switch which_var
            
            case 1
                % same time scales
                infile = 'SEIR_agedep_twodiseases_071321_T5and5.mat';
                load(strcat('sim_data/',infile));
                
                cbf_colors = [15,32,128]/255; % dark blue
                gamma_a1 = params.gamma_a;
                gamma_s1 = params.gamma_s;
                beta_a1 = params.beta_a;
                beta_s1 = params.beta_s;
                R0_agedep1 = results.Rt_agedep(1);
                
                
            case 2
                % same time scales, 4x assortativity
                infile = 'SEIR_agedep_twodiseases_071321_T5and5_ia.mat';
                load(strcat('sim_data/',infile));
                cbf_colors = [15,32,128]/255; % dark blue
                
                gamma_a2 = params.gamma_a;
                gamma_s2 = params.gamma_s;
                beta_a2 = params.beta_a;
                beta_s2 = params.beta_s;
                R0_agedep2 = results.Rt_agedep(1);
                
                
            case 3
                % Ta=8, Ts=5
                infile = 'SEIR_agedep_twodiseases_071321_T5and8.mat';
                load(strcat('sim_data/',infile));
                cbf_colors = [133,192,249]/255; % light blue
                
                gamma_a3 = params.gamma_a;
                gamma_s3 = params.gamma_s;
                beta_a3 = params.beta_a;
                beta_s3 = params.beta_s;
                R0_agedep3 = results.Rt_agedep(1);
                
                
            case 4
                % Ta=8, Ts=5 with 4x assortativity
                infile = 'SEIR_agedep_twodiseases_071321_T5and8_ia.mat';
                load(strcat('sim_data/',infile));
                cbf_colors = [133,192,249]/255; % light blue
                
                gamma_a4 = params.gamma_a;
                gamma_s4 = params.gamma_s;
                beta_a4 = params.beta_a;
                beta_s4 = params.beta_s;
                R0_agedep4 = results.Rt_agedep(1);
                
        end
         
        
        subplot(1,2,1);
        if params.increased_assortativity_yesno==0
            h(counter) = plot(params.t_span, results.p_aa,'Color',cbf_colors,'LineWidth',1); hold on;
            h(counter) = plot(params.t_span, results.p_as,'Color',cbf_colors,'LineWidth',1); hold on;
        else
            h(counter) = plot(params.t_span, results.p_aa,'--','Color',cbf_colors,'LineWidth',2); hold on;
            h(counter) = plot(params.t_span, results.p_as,'--','Color',cbf_colors,'LineWidth',2); hold on;
        end
        
        axis([0 params.t_span(end) 0.5 0.75]);
        xlabel('Time (days)');
        ylabel({'Probability'});
        yticks(0.5:0.05:0.75);
        f1=gca;
        f1.LineWidth = 1;
        f1.FontSize = 16;
        f1.FontWeight = 'normal';
        f1.FontName = 'Times';
        
        title('Susceptible Depletion','FontSize',18);
        
        if counter == 4
            
            text(102,0.66,'$p_{a|a}$','Interpreter','Latex','FontSize',18);
            text(127,0.595,'$p_{a|s}$','Interpreter','Latex','FontSize',18);
            
            txt = {'A'};
            text(0.025,1.025,txt,'Units','normalized','FontSize',16,'FontWeight','bold');
        end
        
    end
    
    
    
    
end


%% same figures -- pa|a and pa|s vs. time -- with intervention
if 1
    % B: with mitigation
    for counter=1:length(switch_over_var)
        
        which_var = switch_over_var(counter);
        
        switch which_var
            
            case 1
                % same time scales - mitigation
                infile = 'SEIR_agedep_twodiseases_071321_T5and5_mit.mat';
                load(strcat('sim_data/',infile));
                cbf_colors = [15,32,128]/255; % dark blue
                
                
            case 2
                
                % same time scales, 4x assortativity - mitigation
                infile = 'SEIR_agedep_twodiseases_071321_T5and5_mit_ia.mat';
                load(strcat('sim_data/',infile));
                cbf_colors = [15,32,128]/255; % dark blue
                
                
            case 3
                
                % Ta=8, Ts=5 - mitigation
                infile = 'SEIR_agedep_twodiseases_071321_T5and8_mit.mat';
                load(strcat('sim_data/',infile));
                cbf_colors = [133,192,249]/255; % light blue
                
                
            case 4
                
                % Ta=8, Ts=5 with 4x assortativity - mitigation
                infile = 'SEIR_agedep_twodiseases_071321_T5and8_mit_ia.mat';
                load(strcat('sim_data/',infile));
                cbf_colors = [133,192,249]/255; % light blue

                
        end
        
        
        
        subplot(1,2,2);
        if params.increased_assortativity_yesno==0
            h(counter) = plot(params.t_span, results.p_aa,'Color',cbf_colors,'LineWidth',1); hold on;
            h(counter) = plot(params.t_span, results.p_as,'Color',cbf_colors,'LineWidth',1); hold on;
        else
            h(counter) = plot(params.t_span, results.p_aa,'--','Color',cbf_colors,'LineWidth',2); hold on;
            h(counter) = plot(params.t_span, results.p_as,'--','Color',cbf_colors,'LineWidth',2); hold on;
        end
        
        axis([0 params.t_span(end) 0.5 0.75]);
        xlabel('Time (days)');
        ylabel({'Probability'});
        yticks(0.5:0.05:0.75);
        % title(['p = ',num2str(proportion_asymp)])
        f1=gca;
        f1.LineWidth = 1;
        f1.FontSize = 16;
        f1.FontWeight = 'normal';
        f1.FontName = 'Times';
        
        title('Intervention','FontSize',18);
        
        if counter==1
            txt = {'B'};
            %     text(5,9.4,txt,
            text(0.025,1.025,txt,'Units','normalized','FontSize',16,'FontWeight','bold');
        end
        

        
        if counter==4
            
            ylim=get(gca,'ylim');
            xlim=get(gca,'xlim');
            text(100,0.6,'$p_{a|a}$','Interpreter','Latex','FontSize',18);
            text(100,0.572,'$p_{a|s}$','Interpreter','Latex','FontSize',18);
            
            legend_char1 = ['$T_a = ', num2str(1/gamma_a1),'$'];
            legend_char2 = ['$T_a = ', num2str(1/gamma_a2),'\, (4\times)$'];
            legend_char3 = ['$T_a = ', num2str(1/gamma_a3),' $'];
            legend_char4 = ['$T_a = ', num2str(1/gamma_a4),'\, (4\times)$'];

            legend(h,{legend_char1,legend_char2,legend_char3,legend_char4}, 'Interpreter','Latex','Location','NorthEast','FontSize',12);
            legend boxoff
        end
        
        
        
        
    end
    
    
    
end


%% save figure
if save_ans_Fig
    
    folder_location = './../Figures_ms_all/supp/';
    saveas(f1,strcat(folder_location,figure_name),'epsc');
    
    fprintf('Figure saved:\n');
    fprintf(strcat(figure_name,'\n\n'));
    
    fprintf('Location:\n');
    fprintf(strcat(folder_location,'\n\n'));
    
else
    
    fprintf('Figure not saved.\n');
    
end