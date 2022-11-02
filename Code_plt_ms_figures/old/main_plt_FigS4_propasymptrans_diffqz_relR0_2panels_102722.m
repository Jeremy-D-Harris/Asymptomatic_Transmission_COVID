%% supp fig - 
% panel A: initial proportion asymptomatic transmission vs. R_0,s/R_0,a
% panel B: change in realized proportion asymptomatic transmission vs. R_0,s/R_0,a

clear all; close all; clc;

%% save figure?
save_ans_Fig = 0;
% 0: don't save
% 1: save

figure_name = 'SuppFig_fixedpropasymp_diffqz_vs_ratioR0_2panels_102722';

%% set up colors
% cbf_colors_db = [15,32,128]/255; % dark blue - same time scales
cbf_colors_v = [169,90,161]/255; % violet - longer time scale of asymptomatic
cbf_colors_lb = [133,192,249]/255; % light blue - even longer time scale of asymptomatc

cbf_colors_vector = [cbf_colors_v;cbf_colors_lb];
%%
p = 0.4; % intrinsic proportion of asymptomatic infections
t_vary = [1,5/6,5/8]; % ratio of T_s/T_a = t

k_vector = linspace(1,1,4); % ratio of beta_s/beta_a = k
% k_discrete = k_vector(1:2:end);

for count=1:length(t_vary)
    
    this_t = t_vary(count);
    
    z_discrete(count,:) = p./(p+(1-p)*k_vector*this_t); % initial proprotion of asymptomatic transmission
%     z_discrete(count,:) = z(count,:);
    
end


% now load files and calculate change in proportion asymptomatic transmission

%% Ta=6
infile = 'SEIR_fixedpropasymp_twodiseases_sameR0s_011722_T5and6.mat';
load(strcat('./sim_data/',infile));
init_prop_asymp_trans_Ta6(1)=results.proportion_asymp_transmission(1);
final_prop_asymp_trans_Ta6(1)=results.proportion_asymp_transmission(end);

%%
infile = 'SEIR_fixedpropasymp_twodiseases_Rs2timesRa_102622_T5and6.mat';
load(strcat('./sim_data/',infile));
init_prop_asymp_trans_Ta6(2)=results.proportion_asymp_transmission(1);
final_prop_asymp_trans_Ta6(2)=results.proportion_asymp_transmission(end);


infile = 'SEIR_fixedpropasymp_twodiseases_Rs3timesRa_102622_T5and6.mat';
load(strcat('./sim_data/',infile));
init_prop_asymp_trans_Ta6(3)=results.proportion_asymp_transmission(1);
final_prop_asymp_trans_Ta6(3)=results.proportion_asymp_transmission(end);


infile = 'SEIR_fixedpropasymp_twodiseases_Rs4timesRa_102622_T5and6.mat';
load(strcat('./sim_data/',infile));
init_prop_asymp_trans_Ta6(4)=results.proportion_asymp_transmission(1);
final_prop_asymp_trans_Ta6(4)=results.proportion_asymp_transmission(end);

change_prop_asymp_trans_Ta6=final_prop_asymp_trans_Ta6-init_prop_asymp_trans_Ta6;


%% Ta=8
infile = 'SEIR_fixedpropasymp_twodiseases_sameR0s_011722_T5and8.mat';
load(strcat('./sim_data/',infile));
init_prop_asymp_trans_Ta8(1)=results.proportion_asymp_transmission(1);
final_prop_asymp_trans_Ta8(1)=results.proportion_asymp_transmission(end);


infile = 'SEIR_fixedpropasymp_twodiseases_Rs2timesRa_102622_T5and8.mat';
load(strcat('./sim_data/',infile));
init_prop_asymp_trans_Ta8(2)=results.proportion_asymp_transmission(1);
final_prop_asymp_trans_Ta8(2)=results.proportion_asymp_transmission(end);


infile = 'SEIR_fixedpropasymp_twodiseases_Rs3timesRa_102622_T5and8.mat';
load(strcat('./sim_data/',infile));
init_prop_asymp_trans_Ta8(3)=results.proportion_asymp_transmission(1);
final_prop_asymp_trans_Ta8(3)=results.proportion_asymp_transmission(end);


infile = 'SEIR_fixedpropasymp_twodiseases_Rs4timesRa_102622_T5and8.mat';
load(strcat('./sim_data/',infile));
init_prop_asymp_trans_Ta8(4)=results.proportion_asymp_transmission(1);
final_prop_asymp_trans_Ta8(4)=results.proportion_asymp_transmission(end);

change_prop_asymp_trans_Ta8=final_prop_asymp_trans_Ta8-init_prop_asymp_trans_Ta8;

diff_q_z_Ta6 = init_prop_asymp_trans_Ta6-z_discrete(2,:);
diff_q_z_Ta8 = init_prop_asymp_trans_Ta8-z_discrete(3,:);

%% plot relationships with proportion asymptomatic transmission
% plot q-z vs. R0s/R0a
f1 = figure(1); set(f1, 'Position', [100 500 850 350]);
subplot(1,2,1);
this_h(1)=plot(1:4,diff_q_z_Ta6,'.-','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
this_h(2)=plot(1:4,diff_q_z_Ta8,'.-','Color',cbf_colors_vector(2,:),'LineWidth',2,'MarkerSize',18); hold on;
axis([0 5 -0.2 0]);
xticks([0:1:5]);
yticks([-0.2:0.05:0]);
xlabel('$\mathcal R_{0,s}/\mathcal R_{0,a}$','Interpreter','Latex');
% xlabel('$\beta_{s}/\beta_{a}$','Interpreter','Latex');
ylabel({'Difference between '; 'Intrinsic and Realized Proportions of'; 'Asymptomatic Transmission, $q(0)-z$'},'Interpreter','Latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

set(f1,'xticklabel',[{'0'},{'1'},{'2'},{'3'},{'4'},{'5'}]);
set(f1,'yticklabel',[{'-0.2'},{'-0.15'},{'-0.1'},{'-0.05'},{'0'}]);

txt = {'A'};
text(0.025,0.95,txt,'Units','normalized','FontSize',14,'FontWeight','bold');

% plot change over time: q(infty)-q(0)
subplot(1,2,2);
this_h(1)=plot(1:4,change_prop_asymp_trans_Ta6,'.-','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
this_h(2)=plot(1:4,change_prop_asymp_trans_Ta8,'.-','Color',cbf_colors_vector(2,:),'LineWidth',2,'MarkerSize',18); hold on;
axis([0 5 0 0.3]);
xticks([0:1:5]);
yticks([0:0.1:0.3]);
xlabel('$\mathcal R_{0,s}/\mathcal R_{0,a}$','Interpreter','Latex');
% xlabel('$\beta_{s}/\beta_{a}$','Interpreter','Latex');
ylabel({'Change in Realized Proportion of'; 'Asymptomatic Transmission, $q(\infty) - q(0)$'},'Interpreter','Latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

set(f1,'xticklabel',[{'0'},{'1'},{'2'},{'3'},{'4'},{'5'}]);
set(f1,'yticklabel',[{'0'},{'0.1'},{'0.2'},{'0.3'}]);


legend_char1 = ['$T_s = 5$, $T_a = 6$'];
legend_char2 = ['$T_s = 5$, $T_a = 8$'];
legend(this_h,{legend_char1,legend_char2}, 'Interpreter','Latex','Location','NorthEast','FontSize',11);
legend box off;

txt = {'B'};
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

