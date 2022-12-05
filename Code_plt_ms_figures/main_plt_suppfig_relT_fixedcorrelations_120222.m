%% supp fig -
% panel A: q-z vs. Ta/Ts
% panel B: change in realized proportion asymptomatic transmission vs. Ta/Ts
% panel C: p(0)-p vs. Ta/Ts
% panel D: change in realized proportion asymptomatic incidence vs. Ta/Ts

clear all; close all; clc;
%
%% save figure?
save_ans_Fig = 0;
% 0: don't save
% 1: save

figure_name = 'Figsupp_fixedcorrelations_relevance_vs_relT_120222';

%% set up colors
% cbf_colors_db = [15,32,128]/255; % dark blue - same time scales
% cbf_colors_v = [169,90,161]/255; % violet - longer time scale of asymptomatic
% cbf_colors_lb = [133,192,249]/255; % light blue - even longer time scale of asymptomatc

cbf_colors_gray = [0.5 0.5 0.5];
cbf_colors_black = [0 0 0];
cbf_colors_vector = [cbf_colors_black;cbf_colors_gray];
%%

% should have defined fixed_p
fixed_p = 0.4;

% now load files and make calculations for each

%% which file

for which_file = 1:4
    
    switch which_file
        
        case 1
            infile = 'SEIR_fixedp_twodiseases_sameR0s_varyTa_120222.mat';
            %             cbf_colors = cbf_colors_vector(1,:);
            %             linestyle = '-';
        case 2
            infile = 'SEIR_assortmixing_twodiseases_sameR0s_varyTa_120222.mat';
            %             cbf_colors = cbf_colors_vector(1,:);
            %             linesstyle = '--';
        case 3
            infile = 'SEIR_fixedp_twodiseases_Rs4timesRa_varyTa_120222.mat';
            %             cbf_colors = cbf_colors_vector(2,:);
            %             linesstyle = '-';
        case 4
            infile = 'SEIR_assortmixing_twodiseases_Rs4timesRa_varyTa_120222.mat';
            %             cbf_colors = cbf_colors_vector(2,:);
            %             linesstyle = '--';
            
    end
    
    
    
    load(strcat('./sim_data/',infile));
    
    %%
    k_vector_relT = params_collect.k_vector_relT;
    for count=1:length(k_vector_relT)
        
        % transmission
        prop_asymp_transmission(count,:)=results_collect(count).proportion_asymp_transmission;
        
        %     diff_q_z_Ta6(count) = prop_asymp_trans_Ta6(count,:)-results_collect(count).proportion_asymp_transmission_z;
        diff_q_z(count) = results_collect(count).proportion_asymp_transmission_q-results_collect(count).proportion_asymp_transmission_z;
        
        change_prop_asymp_transmission(count)=prop_asymp_transmission(count,end)-prop_asymp_transmission(count,1);
        
        % incidence
        proportion_asymp_incidence(count,:)=results_collect(count).proportion_asymp_incidence;
        init_prop_asymp_incidence(count) = proportion_asymp_incidence(count,1);
        
        %     diff_pt_p_sameR0s(count) = params_collect(count).p;
        diff_pt_p(count) = proportion_asymp_incidence(count,1) - fixed_p;
        p_collect(count) = params_collect(count).p;
        
        change_prop_asymp_incidence(count)=proportion_asymp_incidence(count,end)-proportion_asymp_incidence(count,1);
        
        % transmission rates
        beta_a_collect(count) = params_collect(count).beta_a;
        beta_s_collect(count) = params_collect(count).beta_s;
        
        reproduction_number(count) = results_collect(count).Rt_assortmixing(1);
        
    end
    
    plt_reproduction_number(which_file,:) = reproduction_number;
    
    plt_beta_a_collect(which_file,:) = beta_a_collect;
    plt_beta_s_collect(which_file,:) = beta_s_collect;
    
    plt_p_collect(which_file,:) = p_collect;
    
    plt_diff_q_z_sameR0s(which_file,:) = diff_q_z;
    

    plt_change_prop_asymp_trans_sameR0s(which_file,:) = change_prop_asymp_transmission;
    
    plt_diff_pt_p(which_file,:) = diff_pt_p;
    
    plt_init_prop_asymp_incidence(which_file,:) = init_prop_asymp_incidence;
    plt_change_prop_asymp_incidence_sameR0s(which_file,:) = change_prop_asymp_incidence;
    
    
end

%% get specific points
get_relT_pts = [5/8 5/6 1 6/5 8/5];
ind_k_vector_relT(1) = 1;
ind = find(k_vector_relT<get_relT_pts(2));
ind_k_vector_relT(2) = ind(end)+1;
ind = find(k_vector_relT<get_relT_pts(3));
ind_k_vector_relT(3) = ind(end)+1;
ind = find(k_vector_relT<get_relT_pts(4));
ind_k_vector_relT(4) = ind(end)+1;
ind_k_vector_relT(5) = length(k_vector_relT);

% reset due to singularity

for which_file = 1:4
    
    plt_diff_q_z_sameR0s(which_file,ind_k_vector_relT(3)) = 0;
    plt_diff_pt_p(which_file,ind_k_vector_relT(3)) = (plt_diff_pt_p(which_file,ind_k_vector_relT(3)-1)+plt_diff_pt_p(which_file,ind_k_vector_relT(3)+1))/2;
    
    plt_p_collect(which_file,ind_k_vector_relT(3)) = (plt_p_collect(which_file,ind_k_vector_relT(3)-1)+plt_p_collect(which_file,ind_k_vector_relT(3)+1))/2;
    
    plt_beta_a_collect(which_file,ind_k_vector_relT(3)) = (plt_beta_a_collect(which_file,ind_k_vector_relT(3)-1)+plt_beta_a_collect(which_file,ind_k_vector_relT(3)+1))/2;
    plt_beta_s_collect(which_file,ind_k_vector_relT(3)) = (plt_beta_a_collect(which_file,ind_k_vector_relT(3)-1)+plt_beta_s_collect(which_file,ind_k_vector_relT(3)+1))/2;
    
end


%% now plot things !!
for which_file = 1:4
    
    
    %% plot relationships with proportion asymptomatic transmission
    % plot q-z vs. R0s/R0a
    f1 = figure(1); set(f1, 'Position', [100 200 850 700]);
    subplot(2,2,1);
    if which_file == 1
        plot(k_vector_relT,plt_diff_q_z_sameR0s(which_file,:),'Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
        plot(k_vector_relT(ind_k_vector_relT),plt_diff_q_z_sameR0s(which_file,ind_k_vector_relT),'.','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
    elseif which_file == 2
        plot(k_vector_relT,plt_diff_q_z_sameR0s(which_file,:),'--','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
        plot(k_vector_relT(ind_k_vector_relT),plt_diff_q_z_sameR0s(which_file,ind_k_vector_relT),'.','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
    elseif which_file == 3
        plot(k_vector_relT,plt_diff_q_z_sameR0s(which_file,:),'Color',cbf_colors_vector(2,:),'LineWidth',2,'MarkerSize',18); hold on;
        plot(k_vector_relT(ind_k_vector_relT),plt_diff_q_z_sameR0s(which_file,ind_k_vector_relT),'.','Color',cbf_colors_vector(2,:),'LineWidth',2,'MarkerSize',18); hold on;
    else
        plot(k_vector_relT,plt_diff_q_z_sameR0s(which_file,:),'--','Color',cbf_colors_vector(2,:),'LineWidth',2,'MarkerSize',18); hold on;
        plot(k_vector_relT(ind_k_vector_relT),plt_diff_q_z_sameR0s(which_file,ind_k_vector_relT),'.','Color',cbf_colors_vector(2,:),'LineWidth',2,'MarkerSize',18); hold on;
    end
    
    
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
    
    if which_file ==4
        txt = {'A'};
        text(0.025,0.95,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
    % plot change over time: q(infty)-q(0)
    subplot(2,2,2);
    
    if which_file == 1
        this_h(which_file)=plot(k_vector_relT,plt_change_prop_asymp_trans_sameR0s(which_file,:),'Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
        plot(k_vector_relT(ind_k_vector_relT),plt_change_prop_asymp_trans_sameR0s(which_file,ind_k_vector_relT),'.','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
    elseif which_file == 2
        this_h(which_file)=plot(k_vector_relT,plt_change_prop_asymp_trans_sameR0s(which_file,:),'--','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
        plot(k_vector_relT(ind_k_vector_relT),plt_change_prop_asymp_trans_sameR0s(which_file,ind_k_vector_relT),'.','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
    elseif which_file == 3
        this_h(which_file)=plot(k_vector_relT,plt_change_prop_asymp_trans_sameR0s(which_file,:),'Color',cbf_colors_vector(2,:),'LineWidth',2,'MarkerSize',18); hold on;
        plot(k_vector_relT(ind_k_vector_relT),plt_change_prop_asymp_trans_sameR0s(which_file,ind_k_vector_relT),'.','Color',cbf_colors_vector(2,:),'LineWidth',2,'MarkerSize',18); hold on;
    else
        this_h(which_file)=plot(k_vector_relT,plt_change_prop_asymp_trans_sameR0s(which_file,:),'--','Color',cbf_colors_vector(2,:),'LineWidth',2,'MarkerSize',18); hold on;
        plot(k_vector_relT(ind_k_vector_relT),plt_change_prop_asymp_trans_sameR0s(which_file,ind_k_vector_relT),'.','Color',cbf_colors_vector(2,:),'LineWidth',2,'MarkerSize',18); hold on;
        
    end
    axis([5/8 8/5 -0.4 0.4]);
%     xlim([5/8 8/5]);
    xticks([5/8 5/6 1 6/5 8/5]);
    yticks([-0.4:0.2:0.4]);
    xlabel('$T_a/T_s$','Interpreter','Latex');
    ylabel({'Change in Realized Proportion of'; 'Asymptomatic Transmission, $q(\infty) - q(0)$'},'Interpreter','Latex');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    set(f1,'xticklabel',[{'5/8'},{'5/6'},{'1'},{'6/5'},{'8/5'}]);
    set(f1,'yticklabel',[{'-0.4'},{'-0.2'},{'0'},{'0.2'},{'0.4'}]);
    
    
    
    if which_file ==4
        txt = {'B'};
        text(0.025,0.95,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
        
        legend_char1 = ['$\mathcal R_{0,a} = \mathcal R_{0,s}$, $p = 0.40$'];
        legend_char2 = ['$\mathcal R_{0,a} = \mathcal R_{0,s}$, $p_{a|a} = 0.50$, $p_{a|s} = 0.25$'];
        legend_char3 = ['$\mathcal R_{0,a} = 4\,\mathcal R_{0,s}$, $p = 0.40$'];
        legend_char4 = ['$\mathcal R_{0,a} = 4\,\mathcal R_{0,s}$, $p_{a|a} = 0.50$, $p_{a|s} = 0.25$'];
        
        legend(this_h,{legend_char1,legend_char2,legend_char3,legend_char4}, 'Interpreter','Latex','Location','SouthEast','FontSize',10.5);
        legend box off;
        
        
    end
    
    %% plot C-D:relationships with proportion asymptomatic incidence
    % plot q-z vs. R0s/R0a
    f1 = figure(1); set(f1, 'Position', [100 200 850 700]);
    subplot(2,2,3);
    if which_file == 1
%         plot(k_vector_relT,plt_diff_pt_p(which_file,:),'-','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
%         plot(k_vector_relT(ind_k_vector_relT),plt_diff_pt_p(which_file,ind_k_vector_relT),'.','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
        plot(k_vector_relT,plt_init_prop_asymp_incidence(which_file,:),'-','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
        plot(k_vector_relT(ind_k_vector_relT),plt_init_prop_asymp_incidence(which_file,ind_k_vector_relT),'.','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
        
    elseif which_file == 2
%         plot(k_vector_relT,plt_diff_pt_p(which_file,:),'--','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
%         plot(k_vector_relT(ind_k_vector_relT),plt_diff_pt_p(which_file,ind_k_vector_relT),'.','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
        plot(k_vector_relT,plt_init_prop_asymp_incidence(which_file,:),'--','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
        plot(k_vector_relT(ind_k_vector_relT),plt_init_prop_asymp_incidence(which_file,ind_k_vector_relT),'.','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
        
    elseif which_file == 3
%         this_p = plot(k_vector_relT,plt_diff_pt_p(which_file,:),'Color',cbf_colors_vector(2,:),'LineWidth',3,'MarkerSize',18); hold on;
        this_p = plot(k_vector_relT,plt_init_prop_asymp_incidence(which_file,:),'Color',cbf_colors_vector(2,:),'LineWidth',3,'MarkerSize',18); hold on;
        this_p.Color(4) = 0.5;
%         this_p = plot(k_vector_relT(ind_k_vector_relT),plt_diff_pt_p_sameR0s(which_file,ind_k_vector_relT),'.','Color',[cbf_colors_vector(2,:),0.2],'MarkerSize',18); hold on;
%         scatter(k_vector_relT(ind_k_vector_relT),plt_diff_pt_p(which_file,ind_k_vector_relT),'filled','SizeData',75,'MarkerEdgeColor',cbf_colors_vector(2,:),'MarkerFaceColor',cbf_colors_vector(2,:),'MarkerFaceAlpha', .5);
        scatter(k_vector_relT(ind_k_vector_relT),plt_init_prop_asymp_incidence(which_file,ind_k_vector_relT),'filled','SizeData',75,'MarkerEdgeColor',cbf_colors_vector(2,:),'MarkerFaceColor',cbf_colors_vector(2,:),'MarkerFaceAlpha', .5);
        
    elseif which_file == 4
%         plot(k_vector_relT,plt_diff_pt_p(which_file,:),'--','Color',cbf_colors_vector(2,:),'LineWidth',2,'MarkerSize',18); hold on;
%         plot(k_vector_relT(ind_k_vector_relT),plt_diff_pt_p(which_file,ind_k_vector_relT),'.','Color',cbf_colors_vector(2,:),'LineWidth',2,'MarkerSize',18); hold on;
        plot(k_vector_relT,plt_init_prop_asymp_incidence(which_file,:),'--','Color',cbf_colors_vector(2,:),'LineWidth',2,'MarkerSize',18); hold on;
        plot(k_vector_relT(ind_k_vector_relT),plt_init_prop_asymp_incidence(which_file,ind_k_vector_relT),'.','Color',cbf_colors_vector(2,:),'LineWidth',2,'MarkerSize',18); hold on;
    end
    % plot(1,0,'.','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
    % this_h(2)=plot(k_vector_relT,diff_q_z_Ta8,'Color',cbf_colors_vector(2,:),'LineWidth',2,'MarkerSize',18); hold on;
    % plot(k_vector_relT(1:10:end),diff_q_z_Ta8(1:10:end),'.','Color',cbf_colors_vector(2,:),'LineWidth',2,'MarkerSize',18); hold on;
    axis([5/8 8/5 0.25 0.45]);
%     xlim([5/8 8/5]);
    xticks([5/8 5/6 1 6/5 8/5]);
%     yticks([-0.15:0.05:0.05]);
    xlabel('$T_a/T_s$','Interpreter','Latex');
    % xlabel('$\beta_{s}/\beta_{a}$','Interpreter','Latex');
    ylabel({'Initial Proportion of'; 'Asymptomatic Incidence, $p(0)$'},'Interpreter','Latex');
    f1=gca;
    f1.LineWidth = 1;
    f1.FontSize = 14;
    f1.FontWeight = 'normal';
    f1.FontName = 'Times';
    
    set(f1,'xticklabel',[{'5/8'},{'5/6'},{'1'},{'6/5'},{'8/5'}]);
%     set(f1,'yticklabel',[{'-0.12'},{'-0.08'},{'-0.04'},{'0'}]);
    
    if which_file ==4
        txt = {'C'};
        text(0.025,0.95,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
    end
    
    % plot change over time: q(infty)-q(0)
    subplot(2,2,4);
    if which_file == 1
    plot(k_vector_relT,plt_change_prop_asymp_incidence_sameR0s(which_file,:),'Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
    plot(k_vector_relT(ind_k_vector_relT),plt_change_prop_asymp_incidence_sameR0s(which_file,ind_k_vector_relT),'.','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
    elseif which_file == 2
            plot(k_vector_relT,plt_change_prop_asymp_incidence_sameR0s(which_file,:),'--','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
    plot(k_vector_relT(ind_k_vector_relT),plt_change_prop_asymp_incidence_sameR0s(which_file,ind_k_vector_relT),'.','Color',cbf_colors_vector(1,:),'LineWidth',2,'MarkerSize',18); hold on;
        
    elseif which_file == 3
                this_p = plot(k_vector_relT,plt_change_prop_asymp_incidence_sameR0s(which_file,:),'Color',cbf_colors_vector(2,:),'LineWidth',3,'MarkerSize',18); hold on;
                this_p.Color(4) = 0.5;
    scatter(k_vector_relT(ind_k_vector_relT),plt_change_prop_asymp_incidence_sameR0s(which_file,ind_k_vector_relT),'filled','SizeData',75,'MarkerEdgeColor',cbf_colors_vector(2,:),'MarkerFaceColor',cbf_colors_vector(2,:),'MarkerFaceAlpha', .5); hold on;
    else
            plot(k_vector_relT,plt_change_prop_asymp_incidence_sameR0s(which_file,:),'--','Color',cbf_colors_vector(2,:),'LineWidth',2,'MarkerSize',18); hold on;
    plot(k_vector_relT(ind_k_vector_relT),plt_change_prop_asymp_incidence_sameR0s(which_file,ind_k_vector_relT),'.','Color',cbf_colors_vector(2,:),'LineWidth',2,'MarkerSize',18); hold on;
                
    end
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
    
    if which_file ==4
        txt = {'D'};
        text(0.025,0.95,txt,'Units','normalized','FontSize',14,'FontWeight','bold');
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

