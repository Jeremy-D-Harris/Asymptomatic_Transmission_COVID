%% supp fig -
% panel A: initial proportion asymptomatic transmission vs. R_0,s/R_0,a
% panel B: change in realized proportion asymptomatic transmission vs. R_0,s/R_0,a

% clear all; close all; clc;
%
%% save figure?
save_ans_Fig = 0;
% 0: don't save
% 1: save

figure_name = 'Figsupp_diffqz_vs_relT_113022';

%% set up colors
% cbf_colors_db = [15,32,128]/255; % dark blue - same time scales
cbf_colors_v = [169,90,161]/255; % violet - longer time scale of asymptomatic
cbf_colors_lb = [133,192,249]/255; % light blue - even longer time scale of asymptomatc

cbf_colors_gray = [0.5 0.5 0.5];
cbf_colors_black = [0 0 0];
cbf_colors_vector = [cbf_colors_black;cbf_colors_gray];
%%
% now load files and calculate
% change in proportion asymptomatic transmission

%% which file
% infile = 'SEIR_fixedp_twodiseases_sameR0s_varyTa_113022_refine.mat';
infile = 'SEIR_assortmixing_twodiseases_sameR0s_varyTa_113022_refine.mat';
% infile = 'SEIR_fixedp_twodiseases_Rs4timesRa_varyTa_113022_refine.mat';
% infile = 'SEIR_assortmixing_twodiseases_Rs4timesRa_varyTa_113022_refine.mat';




load(strcat('./sim_data/',infile));

%%
k_vector_relT = params_collect.k_vector_relT;
for count=1:length(k_vector_relT)
    
    % transmission
    prop_asymp_trans_sameR0s(count,:)=results_collect(count).proportion_asymp_transmission;
    
    %     diff_q_z_Ta6(count) = prop_asymp_trans_Ta6(count,:)-results_collect(count).proportion_asymp_transmission_z;
    diff_q_z_sameR0s(count) = results_collect(count).proportion_asymp_transmission_q-results_collect(count).proportion_asymp_transmission_z;
    
    change_prop_asymp_trans_sameR0s(count)=prop_asymp_trans_sameR0s(count,end)-prop_asymp_trans_sameR0s(count,1);
    
    % incidence
    prop_asymp_incidence_sameR0s(count,:)=results_collect(count).proportion_asymp_incidence;
    
%     diff_pt_p_sameR0s(count) = params_collect(count).p;
    diff_pt_p_sameR0s(count) = prop_asymp_incidence_sameR0s(count,1) - params_collect(count).p;
    
    change_prop_asymp_incidence_sameR0s(count)=prop_asymp_incidence_sameR0s(count,end)-prop_asymp_incidence_sameR0s(count,1);
end


%% Ta=8
% infile = 'SEIR_fixedpropasymp_twodiseases_Ta8_varyrelR0_smooth_110122.mat';
% load(strcat('./sim_data/',infile));

% %%
% % k_vector_relR0 = params_collect.k_vector_relR0;
% for count=1:length(k_vector_relT)
%     
%     prop_asymp_trans_Rs4timesRa(count,:)=results_collect(count).proportion_asymp_transmission;
%     
%     %     diff_q_z_Ta6(count) = prop_asymp_trans_Ta6(count,:)-results_collect(count).proportion_asymp_transmission_z;
%     diff_q_z_Rs4timesRa(count) = results_collect(count).proportion_asymp_transmission_q-results_collect(count).proportion_asymp_transmission_z;
%     
%     change_prop_asymp_trans_Rs4timesRa(count)=prop_asymp_trans_sameR0s(count,end)-prop_asymp_trans_sameR0s(count,1);
%     
%     
%     %     prop_asymp_trans_Rs4timesRa(count,:)=results_collect(count).proportion_asymp_transmission;
%     %
%     %     diff_q_z_Ta8(count) = results_collect(count).proportion_asymp_transmission_q-results_collect(count).proportion_asymp_transmission_z;
%     %
%     %     change_prop_asymp_trans_Ta8(count)=prop_asymp_trans_Rs4timesRa(count,end)-prop_asymp_trans_Rs4timesRa(count,1);
%     
% end

%% plot relationships with proportion asymptomatic transmission
% plot q-z vs. R0s/R0a
f1 = figure(1); set(f1, 'Position', [100 200 850 700]);
subplot(2,2,1);
this_h(1)=plot(k_vector_relT,diff_q_z_sameR0s,'Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;

get_relT_pts = [5/8 5/6 1 6/5 8/5];
ind_k_vector_relT(1) = 1;
ind = find(k_vector_relT<get_relT_pts(2));
ind_k_vector_relT(2) = ind(end)+1;
ind = find(k_vector_relT<get_relT_pts(3));
ind_k_vector_relT(3) = ind(end)+1;
ind = find(k_vector_relT<get_relT_pts(4));
ind_k_vector_relT(4) = ind(end)+1;
ind_k_vector_relT(5) = length(k_vector_relT);

plot(k_vector_relT(ind_k_vector_relT),diff_q_z_sameR0s(ind_k_vector_relT),'.','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
% plot(1,0,'.','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
% this_h(2)=plot(k_vector_relT,diff_q_z_Ta8,'Color',cbf_colors_vector(2,:),'LineWidth',2,'MarkerSize',18); hold on;
% plot(k_vector_relT(1:10:end),diff_q_z_Ta8(1:10:end),'.','Color',cbf_colors_vector(2,:),'LineWidth',2,'MarkerSize',18); hold on;
% axis([0 5 -0.06 0]);
xlim([5/8 8/5]);
xticks([5/8 5/6 1 6/5 8/5]);
% yticks([-0.06:0.02:0]);
xlabel('$T_a/T_s$','Interpreter','Latex');
% xlabel('$\beta_{s}/\beta_{a}$','Interpreter','Latex');
ylabel({'Difference between '; 'Intrinsic and Realized Proportions of'; 'Asymptomatic Transmission, $q(0)-z$'},'Interpreter','Latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

set(f1,'xticklabel',[{'5/8'},{'5/6'},{'1'},{'6/5'},{'8/5'}]);
% set(f1,'yticklabel',[{'-0.2'},{'-0.15'},{'-0.1'},{'-0.05'},{'0'}]);

txt = {'A'};
text(0.025,0.95,txt,'Units','normalized','FontSize',14,'FontWeight','bold');

% plot change over time: q(infty)-q(0)
subplot(2,2,2);
this_h(1)=plot(k_vector_relT,change_prop_asymp_trans_sameR0s,'Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
plot(k_vector_relT(ind_k_vector_relT),change_prop_asymp_trans_sameR0s(ind_k_vector_relT),'.','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
% plot(1,0,'.','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
% this_h(2)=plot(k_vector_relT,change_prop_asymp_trans_Ta8,'Color',cbf_colors_vector(2,:),'LineWidth',2,'MarkerSize',18); hold on;
% plot(k_vector_relT(1:10:end),change_prop_asymp_trans_Ta8(1:10:end),'.','Color',cbf_colors_vector(2,:),'LineWidth',2,'MarkerSize',18); hold on;
% axis([0 5 0 0.3]);
xlim([5/8 8/5]);
xticks([5/8 5/6 1 6/5 8/5]);
% yticks([0:0.1:0.3]);
xlabel('$T_a/T_s$','Interpreter','Latex');
ylabel({'Change in Realized Proportion of'; 'Asymptomatic Transmission, $q(\infty) - q(0)$'},'Interpreter','Latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

set(f1,'xticklabel',[{'5/8'},{'5/6'},{'1'},{'6/5'},{'8/5'}]);
% set(f1,'yticklabel',[{'0'},{'0.1'},{'0.2'},{'0.3'}]);


% legend_char1 = ['$T_s = 5$, $T_a = 6$'];
% legend_char2 = ['$T_s = 5$, $T_a = 8$'];
% legend(this_h,{legend_char1,legend_char2}, 'Interpreter','Latex','Location','NorthEast','FontSize',11);
% legend box off;

txt = {'B'};
text(0.025,0.95,txt,'Units','normalized','FontSize',14,'FontWeight','bold');

%% plot C-D:relationships with proportion asymptomatic incidence
% plot q-z vs. R0s/R0a
f1 = figure(1); set(f1, 'Position', [100 200 850 700]);
subplot(2,2,3);
this_h(1)=plot(k_vector_relT,diff_pt_p_sameR0s,'Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
plot(k_vector_relT(ind_k_vector_relT),diff_pt_p_sameR0s(ind_k_vector_relT),'.','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
% plot(1,0,'.','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
% this_h(2)=plot(k_vector_relT,diff_q_z_Ta8,'Color',cbf_colors_vector(2,:),'LineWidth',2,'MarkerSize',18); hold on;
% plot(k_vector_relT(1:10:end),diff_q_z_Ta8(1:10:end),'.','Color',cbf_colors_vector(2,:),'LineWidth',2,'MarkerSize',18); hold on;
% axis([0 5 -0.06 0]);
xlim([5/8 8/5]);
xticks([5/8 5/6 1 6/5 8/5]);
% yticks([-0.06:0.02:0]);
xlabel('$T_a/T_s$','Interpreter','Latex');
% xlabel('$\beta_{s}/\beta_{a}$','Interpreter','Latex');
ylabel({'Difference between '; 'Intrinsic and Realized Proportions of'; 'Asymptomatic Incidence, $p(0)-p$'},'Interpreter','Latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

set(f1,'xticklabel',[{'5/8'},{'5/6'},{'1'},{'6/5'},{'8/5'}]);
% set(f1,'yticklabel',[{'-0.2'},{'-0.15'},{'-0.1'},{'-0.05'},{'0'}]);

txt = {'C'};
text(0.025,0.95,txt,'Units','normalized','FontSize',14,'FontWeight','bold');

% plot change over time: q(infty)-q(0)
subplot(2,2,4);
this_h(1)=plot(k_vector_relT,change_prop_asymp_incidence_sameR0s,'Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
plot(k_vector_relT(ind_k_vector_relT),change_prop_asymp_incidence_sameR0s(ind_k_vector_relT),'.','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
% plot(1,0,'.','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
% this_h(2)=plot(k_vector_relT,change_prop_asymp_trans_Ta8,'Color',cbf_colors_vector(2,:),'LineWidth',2,'MarkerSize',18); hold on;
% plot(k_vector_relT(1:10:end),change_prop_asymp_trans_Ta8(1:10:end),'.','Color',cbf_colors_vector(2,:),'LineWidth',2,'MarkerSize',18); hold on;
% axis([0 5 0 0.3]);
xlim([5/8 8/5]);
xticks([5/8 5/6 1 6/5 8/5]);
% yticks([0:0.1:0.3]);
xlabel('$T_a/T_s$','Interpreter','Latex');
ylabel({'Change in Realized Proportion of'; 'Asymptomatic Incidence, $p(\infty) - p(0)$'},'Interpreter','Latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

set(f1,'xticklabel',[{'5/8'},{'5/6'},{'1'},{'6/5'},{'8/5'}]);
% set(f1,'yticklabel',[{'0'},{'0.1'},{'0.2'},{'0.3'}]);


txt = {'D'};
text(0.025,0.95,txt,'Units','normalized','FontSize',14,'FontWeight','bold');


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

