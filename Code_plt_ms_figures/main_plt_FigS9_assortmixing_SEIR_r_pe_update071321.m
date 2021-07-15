
%% Figure S9 - paramtrize model wrt p(a|a) and p(a|s) to find:
% (1) initial proportion asymptomatic incidence
% (2) exponential growth rate
% beta_a = beta_s = beta
clear all; close all; clc;


%% want to save?
save_Fig = 0;
% 0: don't save
% 1: save

figure_name='FigS9_parametrise_assortmixing_071321';

%% which set of time scales?
% which_timescales = 2; % 1,2
% 1: same time scales: Ta=Ts=5 days
% 2: longer time scales of asymptomatic transmission: Ta=8,Ts=5 days
f1=figure(1); set(gcf,'Position',[50 400 1100 500]);
for which_timescales=1:2
    
    if which_timescales==1
        filename = 'parametrise_assortmixing_same_samebetas_SEIR_071321';
    else
        filename = 'parametrise_assortmixing_diff_samebetas_SEIR_071321';
    end
    
    load(strcat('./sim_data/',filename));
    
    %% plot six panels
    
    subplot(2,3,1+3*(which_timescales-1));
    
    if which_timescales==1
        % x variable: goes down the rows
        % y variable: goes across columns
        s1 = pcolor(params.vary_p_as_plt,params.vary_p_aa_plt,results.r_assortmixing_plt); hold on;
    else
        s1 = pcolor(params.vary_p_as_plt,params.vary_p_aa_plt,results.r_assortmixing_plt); hold on;

        plot(results.contour_r_filtered(1,:),results.contour_r_filtered(2,:),'k--','LineWidth',2.5); hold on;
        
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
    title('Exponential Growth Rate','FontSize',12);
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    if which_timescales==1
        txt = {'A'};
    else
        txt = {'D'};
    end
    text(-0.15,1.125,txt,'Units','normalized','FontSize',16,'FontWeight','bold');
    
    figure(1);
    subplot(2,3,2+3*(which_timescales-1));
    
    % x variable: goes down the rows
    % y variable: goes across columns
    s2 = pcolor(params.vary_p_as_plt,params.vary_p_aa_plt,results.eigen_prop_asymp_incidence_plt); hold on;
    
    plot(results.contour_pe_filtered(1,:),results.contour_pe_filtered(2,:),'k','LineWidth',2); hold on;
    
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
    title('Prop. Asymp. Incidence','FontSize',12);
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    if which_timescales==1
        txt = {'B'};
    else
        txt = {'E'};
    end
    text(-0.15,1.125,txt,'Units','normalized','FontSize',16,'FontWeight','bold');
    
    sparse_par=2;
    figure(1);
    subplot(2,3,3+3*(which_timescales-1));
    plot(params.vary_p_aa,params.vary_p_aa,'Color',[1/2 1/2 1/2]); hold on;
    if which_timescales==1
        
        this_p(1)=plot(results.contour_pe_filtered(1,1:sparse_par:end),results.contour_pe_filtered(2,1:sparse_par:end),'k','LineWidth',2); hold on;
        this_p(2)=plot(results.contour_pe_pas_ia,results.contour_pe_paa_ia,'k.','MarkerSize',20); hold on;
    else
        
        this_q(1)=plot(results.contour_r_filtered(1,1:sparse_par:end),results.contour_r_filtered(2,1:sparse_par:end),'k--','LineWidth',2.5); hold on;
        this_q(2)=plot(results.contour_pe_filtered(1,1:sparse_par:end),results.contour_pe_filtered(2,1:sparse_par:end),'k','LineWidth',2); hold on;
        this_q(3)=plot(results.contour_pe_pas_ia,results.contour_pe_paa_ia,'k.','MarkerSize',20); hold on;
    end
    
    axis([params.lower_bound params.upper_bound params.lower_bound params.upper_bound]);
    xlabel('$p(a|s)$','Interpreter','Latex'); ylabel('$p(a|a)$','Interpreter','Latex');
    xticks(0:0.2:1); yticks(0:0.2:1);
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    if which_timescales==1
        legend(this_p,{'$p_e = 0.40$',['$p(a|s)= $',num2str(results.contour_pe_pas_ia,'%2.2f')]},'Interpreter','Latex','Location','SouthEast');
    else
        legend(this_q,{'$r = 0.14$','$p_e = 0.40$',['$p(a|s)= $',num2str(results.contour_pe_pas_ia,'%2.2f')]},'Interpreter','Latex','Location','SouthEast');
        
    end
    legend boxoff
    
    if which_timescales==1
        txt = {'C'};
    else
        txt = {'F'};
    end
    text(-0.125,1.125,txt,'Units','normalized','FontSize',16,'FontWeight','bold');
    
end

%% save figure
if save_Fig
    
    folder_location = './../Figures_ms_all/supp/';
    saveas(f1,strcat(folder_location,figure_name),'epsc');
    
    fprintf('File saved:\n');
    fprintf(strcat(figure_name,'\n\n'));
    
    fprintf('Location:\n');
    fprintf(strcat(folder_location,'\n\n'));
    
else
    
    fprintf('Figure not saved.\n');
    
end
