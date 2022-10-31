%% supp fig - initial proportion asymptomatic transmission vs. R_0,s/R_0,a

clear all; close all; clc;

%% save figure?
save_ans_Fig = 0;
% 0: don't save
% 1: save

figure_name = 'SuppFig_fixedpropasymp_zvsratiobeta_102622';

%% set up colors
cbf_colors_db = [15,32,128]/255; % dark blue - same time scales
cbf_colors_v = [169,90,161]/255; % violet - longer time scale of asymptomatic
cbf_colors_lb = [133,192,249]/255; % light blue - even longer time scale of asymptomatc

cbf_colors_vector = [cbf_colors_db;cbf_colors_v;cbf_colors_lb];
%%
p = 0.4; % intrinsic proportion of asymptomatic infections
t_vary = [1,5/6,5/8]; % ratio of T_s/T_a = t
beta_ratio_sims = [1,1.196149217809868,1.600243383024034];

k_vector = linspace(0,5,51); % ratio of beta_s/beta_a = k
k_discrete = k_vector(11:10:41);

for count=1:length(t_vary)
    
    this_t = t_vary(count);
    %     this_beta_ratio = beta_ratio_sims(count)
    
    %     this_ind = find(k_vector>this_beta_ratio);
    %     this_k = k_vector(this_ind(1));
    %     k_discrete(count) = this_k;
    
    z(count,:) = p./(p+(1-p)*k_vector*this_t); % initial proprotion of asymptomatic transmission
    z_discrete(count,:) = z(count,11:10:41);
    %     z_discrete(count) = p./(p+(1-p)*this_k*this_t);
    
end

%% plot figure
figure(1);
plot(k_vector,p*ones(size(k_vector)),'k--','LineWidth',1); hold on;
for counter=1:length(t_vary)
    
    this_t = t_vary(counter);
    this_p(counter)=plot(k_vector,z(counter,:),'Color',cbf_colors_vector(counter,:),'LineWidth',2); hold on;
    plot(k_discrete,z_discrete(counter,:),'.','Color',cbf_colors_vector(counter,:),'MarkerSize',20); hold on;
    this_p(counter).Color(4) = 1-0.18*(counter);
end
%     this_p=semilogy(params.t_span,results.Rt_fixedpropasymp,'Color',cbf_colors,'LineWidth',2); hold on;
%     this_p.Color(4) = 1-0.18*(counter);

axis([0 5 0 1]);
xticks([0:0.5:5]); yticks([0:0.1:1]);
xlabel('$\beta_s/\beta_a$','Interpreter','Latex');
% xlabel('$\beta_{s}/\beta_{a}$','Interpreter','Latex');
ylabel({'Intrinsic Proportion of'; 'Asymptomatic Transmission, $z$'},'Interpreter','Latex');
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

set(f1,'xticklabel',[{'0'},{''},{'1'},{''},{'2'},{''},{'3'},{''},{'4'},{''},{'5'}]);
set(f1,'yticklabel',[{'0'},{''},{'0.2'},{''},{'0.4'},{''},{'0.6'},{''},{'0.8'},{''},{'1'}]);

legend_char1 = ['$T_s = T_a = 5$ days'];
legend_char2 = ['$T_s = 5$ days, $T_a = 6 $ days'];
legend_char3 = ['$T_s = 5$ days, $T_a = 8$ days'];
legend(this_p,{legend_char1,legend_char2,legend_char3}, 'Interpreter','Latex','Location','NorthEast','FontSize',11);
legend box off

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

