
%% Figure S4 - plot generation interval distribution 
% for Te=3days and Ta=5,6,8days

clear all; close all; clc;


%% save Figure ?
save_ans_Fig = 0;
% 0: don't save
% 1: save

figure_name = 'FigS4_GIdistributions_071321';

% dark blue, violet, light blue
cat_cbf_colors = [15,32,128;169,90,161;133,192,249]/255;


% time units in days
dt=0.1;
t_span = 25;
time = 0:0.1:t_span; 
gamma_e = 1/3;
gamma_s=1/5;

gamma_a1=1/5; gamma_a2=1/6; gamma_a3=1/8;
vary_gamma_a = [gamma_a1 gamma_a2 gamma_a3];
% % vary_gamma_a = [gamma_a1 gamma_a2];
% vary_gamma_a = [gamma_a1];

params.gamma_s = gamma_s;
params.gamma_e = gamma_e;


f1 = figure(1); set(gcf,'Position',[240   300   560   420]);
for count = 1:length(vary_gamma_a)
    
    this_gamma_a = vary_gamma_a(count);
    params.gamma_a = this_gamma_a;
    
    gen_interval_dist(count,:) = this_gamma_a*gamma_e/(gamma_e-this_gamma_a)*(exp(-this_gamma_a*time)-exp(-gamma_e*time));
%     symp_gen_interval_dist(1,:) = this_gamma_s*gamma_e/(gamma_e-this_gamma_s)*(exp(-gamma_s*time)-exp(-gamma_e*time));
    
    if 1
        
        p(count) = plot(time, gen_interval_dist(count,:), 'Color',cat_cbf_colors(count,:),'LineWidth',2); hold on;
        h = area(time,gen_interval_dist(count,:),'FaceColor',cat_cbf_colors(count,:));
        % Choose a number between 0 (invisible) and 1 (opaque) for facealpha.  
        set(h,'facealpha',0.1)
        axis([0 t_span 0 0.1]);
        yticks([0:0.02:0.1]);
        xlabel('Time (days)');
        ylabel('Probability');
        f1=gca;
        f1.LineWidth = 1;
        f1.FontSize = 14;
        f1.FontWeight = 'normal';
        f1.FontName = 'Times';
        
        
        
    end
end

title('Generation Interval Distributions','FontSize',16);
leg_char1=['$T_a = ', num2str(1/gamma_a1),'$'];
leg_char2=['$T_a = ', num2str(1/gamma_a2),'$'];
leg_char3=['$T_a = ', num2str(1/gamma_a3),'$'];
legend(p,{leg_char1,leg_char2,leg_char3}, 'Location','NorthEast','Interpreter','Latex','FontSize',14);
% legend(p,{leg_char1,leg_char2}, 'Location','NorthEast','Interpreter','Latex','FontSize',12);
legend boxoff


%% save figure
if save_ans_Fig
    
    folder_location = './../Figures_ms_all/supp/';
    saveas(f1,strcat(folder_location,figure_name),'epsc');
    
    fprintf('Figure saved:\n');
    fprintf(strcat(figure_name,'\n\n'));
    
    fprintf('Location:\n');
    fprintf(strcat(folder_location,'\n\n'));
    
end
