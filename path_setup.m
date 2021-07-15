
%% Setup the path to access all the codes, data and packages.

% change the main path variable to the directory where 
% the project is located. So, write: '/XX/Asymptomatic_Spread_COVID19'
% where XX is location of the project

main_path = '/XX/Asymptomatic_Spread_COVID19/';

folder_package = [main_path,'Code_plt_ms_figures'];
addpath(genpath(folder_package));

folder_shared_codes = [main_path,'Code_sims'];
addpath(genpath(folder_shared_codes));

folder_ECLIP_codes = [main_path,'Data'];
addpath(genpath(folder_ECLIP_codes));

folder_ECLIP_figs = [main_path,'Figures_ms_all'];
addpath(genpath(folder_ECLIP_figs));

folder_data = [main_path,'Manuscript_codereview'];
addpath(genpath(folder_data));
