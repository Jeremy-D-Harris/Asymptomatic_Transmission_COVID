ou % function void = main_plt_agedep_data_update070121(void)

% code to plot age-dependent data:
% (A) contact matrix,
% (B) susceptibility/symptomaticity estimates,
% (C) Age distribution of Wuhan population

%% clear and load data
clear all; close all; clc;


%% want to save figure?
save_fig_ans = 0;
% 0: don't save
% 1: save

if save_fig_ans==1
    figure_name = 'FigS12_agedep_data_071321';
end

%% load age-dependent data
load('agedep_data_update071321.mat');
% age distribution of Wuhan (2016)
% age groups are
%  1     2      3      4      5      6      7      8
% 0-9, 10-19, 20-29, 30-39, 40-49, 50-59, 60-69, 70-79
age_group_numbers = data_Wuhan_Davies.age_group_numbers;
age_group_numbers_plus9 = [age_group_numbers;9];

susceptibility_infection= data_Wuhan_Davies.susc_infection;
probability_sympomatic= data_Wuhan_Davies.prob_symp;
contact_matrix = data_Wuhan_Davies.contact_matrix;

% age distribution as proportions
age_distribution = data_Wuhan_Davies.age_distribution/(sum(data_Wuhan_Davies.age_distribution));

f1=figure(1);
set(gcf,'Position',[1250 400 1250 350]);

%% contact matrix, Wuhan baseline estimates
figure(1); subplot(1,3,1);

% x variable: goes down the rows
% y variable: goes across columns
% in contact matrices (in which rows: participants, columns: their contact)
s1 = pcolor(age_group_numbers_plus9,age_group_numbers_plus9,contact_matrix); hold on;
%     axis([1 9 1 9])
set(s1, 'EdgeColor', 'none');
caxis(gca,[0,5]);
% shading interp
cb1 = colorbar;
set(get(cb1,'title'),'string','# contacts per day');
cb1.Limits = [0 4];

cv2 = [2 2 135]/255;
cv1 = [230 230 255]/255;
t = transpose(linspace(0, 1, 101));
for n=1:length(t)
    cmp_lint(n,:) = (1-(t(n))^(1/1.75))*cv1+(t(n)^(1/1.75))*cv2;
end

colormap(cmp_lint);
title('Baseline');
xlabel('Age group: particpant'); ylabel('Age group: contact');


%     xtickangle(90);
shift_ticks_Y = 1.5:(length(age_group_numbers)+0.5);
shift_ticks_X = 1.5:(length(age_group_numbers)+0.5);
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';


f1.YTickLabel = {'0-9', '10-19', '20-29', '30-39', '40-49', '50-59', '60-69', '70+'};
f1.YTick = shift_ticks_Y;
f1.XTick = shift_ticks_X;

xticks(shift_ticks_X);
xtickangle(60);
f1.XTickLabel = {'0-9', '10-19', '20-29', '30-39', '40-49', '50-59', '60-69', '70+'};

txt = {'A'};
text(0.025,1.05,txt,'Units','normalized','FontSize',16,'FontWeight','bold');


%% susceptibility estimates by age (from Davies et al.)

prob_symp_age = probability_sympomatic;
results.prob_symp_age = prob_symp_age;

% prob_asymp_age = 1- prob_symp_age;
% reults.prob_asymp_age = prob_asymp_age;

prob_infection = susceptibility_infection;
results.prob_symp_age = prob_infection;


figure(1); subplot(1,3,2);
combine_prob_symp_infeciton = [prob_infection,prob_symp_age];

figure(1);
% bar(age_group_numbers, prob_symp_age); hold on;
% p = bar(age_group_numbers, combine_prob_symp_infeciton,'grouped'); hold on;
p(1) = plot(age_group_numbers, prob_symp_age,'.-','Color',[0.65 0.65 0.65],'LineWidth',1.5, 'MarkerSize',20); hold on;
p(2) = plot(age_group_numbers, prob_infection,'.--','Color',[0.65 0.65 0.65],'LineWidth',1.5, 'MarkerSize',20); hold on;
xlabel('Age group'); ylabel('Probability');
axis([0.5 8.5 0 1]);
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

% ticks_Y = age_group_numbers;
ticks_X = age_group_numbers;
f1.XTick = ticks_X ;
xticks(ticks_X );
xtickangle(60);
f1.XTickLabel = {'0-9', '10-19', '20-29', '30-39', '40-49', '50-59', '60-69', '70+'};


legend(p,{'Susceptibility to infection','Probability symptomatic infection'}, 'Location','SouthEast');
legend boxoff


txt = {'B'};
text(0.025,1.05,txt,'Units','normalized','FontSize',16,'FontWeight','bold');



%% population age distribution of Wuhan (2016)
figure(1); subplot(1,3,3);
% bar(age_group_numbers, prob_symp_age); hold on;
bar(age_group_numbers, age_distribution); hold on;
xlabel('Age group'); ylabel('Fraction of population');
axis([0 9 0 0.25]);
f1=gca;
f1.LineWidth = 1;
f1.FontSize = 14;
f1.FontWeight = 'normal';
f1.FontName = 'Times';

% {'0', '10', '20', '30', '40', '50', '60', '>70'};
xticks(age_group_numbers);
xtickangle(60);
% f2.XTickLabel = {'0','5', '10', '15', '20', '25', '30', '35', '40', '45', '50', '55','60', '65', '70','>75'};
f1.XTickLabel = {'0-9', '10-19', '20-29', '30-39', '40-49', '50-59', '60-69', '70+'};


txt = {'C'};
text(0.025,1.05,txt,'Units','normalized','FontSize',16,'FontWeight','bold');


%% save figure
if save_fig_ans==1
    
    folder_location = '../../Figures_ms_all/supp/';
    saveas(f1,strcat(folder_location,figure_name),'epsc');
    
    fprintf('Figure saved:\n'); % want to be close to 25 days in
    fprintf(strcat(figure_name,'\n\n'));
    
    fprintf('Location:\n'); % want to be close to 25 days in
    fprintf(strcat(folder_location,'\n\n'));
    
end

if 0
    %% also plot original contact matrix - Wuhan baseline
    
    % pad with zeros
    contact_matrix_original = zeros(15,15);
    contact_matrix_original(1:14,1:14)=data_Wuhan_Davies.contact_matrix_original;
    
    figure(2);
    % x variable: goes down the rows
    % y variable: goes across columns
    % in contact matrices (in which rows: participants, columns: their contact)
    s2 = pcolor(transpose(1:15),transpose(1:15),contact_matrix_original); hold on;
    %     axis([1 9 1 9])
    set(s2, 'EdgeColor', 'none');
    caxis(gca,[0,5]);
    % shading interp
    cb1 = colorbar;
    set(get(cb1,'title'),'string','# contacts per day');
    cb1.Limits = [0 4];
    
    cv2 = [2 2 135]/255;
    cv1 = [230 230 255]/255;
    t = transpose(linspace(0, 1, 101));
    for n=1:length(t)
        cmp_lint(n,:) = (1-(t(n))^(1/1.75))*cv1+(t(n)^(1/1.75))*cv2;
    end
    
    colormap(cmp_lint);
    % title('Baseline');
    xlabel('Age group: particpant'); ylabel('Age group: contact');
    
    
    %     xtickangle(90);
    shift_ticks_Y = 1.5:14.5;
    shift_ticks_X = 1.5:14.5;
    f2=gca;
    f2.LineWidth = 1;
    f2.FontSize = 14;
    f2.FontWeight = 'normal';
    f2.FontName = 'Times';
    
    
    f2.YTickLabel = {'0-4', '5-9', '10-14', '15-19', '20-24',...
        '25-29', '30-34','35-39', '40-44',...
        '45-49', '50-54', '55-59', '60-64', '65+'};
    f2.YTick = shift_ticks_Y;
    f2.XTick = shift_ticks_X;
    
    xticks(shift_ticks_X);
    xtickangle(60);
    f2.YTickLabel = {'0-4', '5-9', '10-14', '15-19', '20-24',...
        '25-29', '30-34','35-39', '40-44',...
        '45-49', '50-54', '55-59', '60-64', '65+'};
    
end