% function void = main_create_datastructure(void)

clear all; close all; clc;


%% want to save?
save_ans = 1;
% 0: don't save
% 1: save

filename = 'agedep_data_update071321.mat';

% define the age group numbers 1-8:
%  1     2      3      4      5      6      7      8
% 0-9, 10-19, 20-29, 30-39, 40-49, 50-59, 60-69, 70-79
num_agegroups=8;
age_group_numbers = transpose(1:num_agegroups);

%% load csv files from data directory
data_place = './../../Data/';


filename_contact = 'contactmatrix_Wuhan_baseline_ZhangetalTableS8.csv';
contact_matrix_Wuhan=load(strcat(data_place,filename_contact));
fprintf('Loaded: \n');
fprintf(strcat(filename_contact, '\n\n'));

filename_susc = 'susc_symp_byage_DaviesetalFig4.csv';
susc_symp_byage=load(strcat(data_place,filename_susc));
fprintf('Loaded: \n');
fprintf(strcat(filename_susc, '\n\n'));

filename_age = 'agedistribution_Wuhan.csv';
agedistribution_Wuhan=load(strcat(data_place,filename_age));
fprintf('Loaded: \n');
fprintf(strcat(filename_age, '\n\n'));

% separate susceptibility and probability of symptomatic infection estimates
susc_infection = susc_symp_byage(:,1); % susceptibility
prob_symp= susc_symp_byage(:,2); % probability of symptomatic infection


%% local variables to structure variables

data_Wuhan_Davies.age_group_numbers = age_group_numbers;


data_Wuhan_Davies.susc_infection = susc_infection;

data_Wuhan_Davies.prob_symp = prob_symp;

data_Wuhan_Davies.age_distribution = agedistribution_Wuhan;

% original contact matrix from Zhang et al. 2020
data_Wuhan_Davies.contact_matrix_original = contact_matrix_Wuhan;


%% reduce the contact matrix from 14x14 to 8x8

% 9x9, pad with zeros for plotting purposes
M = zeros(num_agegroups+1,num_agegroups+1);

% fill in rows 1-6
for row_count=1:6
    for col_count=1:8
        
        if col_count < 7
            % fill in 6x6 block by reducing 1-12 to 1-6 by averaging over 4 contact rates
            M(row_count,col_count) = sum(sum(contact_matrix_Wuhan((2*row_count-1):(2*row_count),(2*col_count-1):(2*col_count))))/4;
            
        elseif col_count==7
            % fill in 1-6x 7 by averaging over 2 contact rates
            M(row_count,col_count) = sum(sum(contact_matrix_Wuhan((2*row_count-1):(2*row_count),2*col_count-1)))/2; % row 1-12, col 13
            
        else
            % fill in 1-6x 8 by averaging over 2 contact rates
            M(row_count,col_count) = sum(sum(contact_matrix_Wuhan((2*row_count-1):(2*row_count),2*col_count-2)))/2; % row 1-12, col 14
        end
    end
end

% fill in rows 7-8
row_count=7;
for col_count=1:8
    
    if col_count < 7
        % reduce row 7, average over 2 contact rates
        M(row_count,col_count) = sum(contact_matrix_Wuhan((2*row_count-1),(2*col_count-1):(2*col_count)))/2;
        M(row_count+1,col_count) = sum(contact_matrix_Wuhan((2*row_count),(2*col_count-1):(2*col_count)))/2;
    elseif col_count==7
        % row 7, col 7
        M(row_count,col_count) = sum(contact_matrix_Wuhan((2*row_count-1),(2*col_count-1))); % row 13, col 13
        % row 8, col 7
        M(row_count+1,col_count) = sum(contact_matrix_Wuhan((2*row_count),(2*col_count-1))); % row 14, col 13
    else
        
        % row 7, col 8
        M(row_count,col_count) = sum(contact_matrix_Wuhan((2*row_count-1),(2*col_count-2))); % row 13, col 14
        % row 8, col 8
        M(row_count+1,col_count) = sum(contact_matrix_Wuhan((2*row_count),(2*col_count-2))); % row 14, col 14
        
    end
    
end


data_Wuhan_Davies.contact_matrix = M;


%%
% save data structure
if save_ans==1
    
    save(filename,'data_Wuhan_Davies');
    
    fprintf('Data saved to file: \n'); 
    fprintf(strcat(filename,'\n\n'));
    
else
    
    fprintf('Data not saved. \n'); 
    
end


fprintf('You can plot the data using the function: \n'); 
fprintf('main_plt_agedep_data_ ... \n');
