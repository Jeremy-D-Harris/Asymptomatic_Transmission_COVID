# How time-scale differences in asymptomatic and symptomatic transmission shape SARS-CoV-2 outbreak dynamics
Jeremy Harris, Sang Woo Park, Jonathan Dushoff, Joshua S. Weitz, April 2022

 -- updated 04/21/22 by jdh --

**Code for:** JD Harris, SW Park, J Dushoff, and JS Weitz. "How time-scale differences in asymptomatic and symptomatic transmission shape SARS-CoV-2 outbreak dynamics." This repository contains all the MATLAB codes for running the three SEIR models of asymptomatic transmission: <br>
  (1) Fixed proportion of asymptomatic incidence (base model); <br>
  (2) Correlations between transmission and disease statuses; <br>
  (3) Assortative mixing by age using age-dependent contact rates and clinical severity estimates by age.

A preprint of the manuscript can be found on MedRxiv: [DOI](https://www.medrxiv.org/content/10.1101/2022.04.21.22274139v1)

**Instructions:** 
MATLAB was used to parametrise models, run simulations, and plot figures -- with the exception of an R-script to plot CDC surveillance data on mortality and median age of infection (May-Aug 2020). Once the Github repository is downloaded, navigate to the subdirectory 'Code_plt_ms_figures' to plot the figures in the manuscript by running the appropriate function. Simulation data can be reproduced by navigating to 'Code_sims' and running appropriate functions. See below for subfolder descriptions.

**Folder descriptions:** <br>

- **Code_sims:** All code to produce figures; subfolders organized by model with models (1)-(3) described above: <br>
  (1) 'fixedpropasymp_code' <br>
  (2) 'assortmixing_code' <br>
  (3) 'agedep_code'

Within each of these folders, the main files parametrise and simulate the model. There are several choices for the user at the top of main file scripts. For instance, in 'main_sim_fixedpropasymp_SEIR_twodiseases_sameR0s_update011722.m' the the first choice for the user is to save the simulation data using the variable 'save_ans': 0 means don't save and 1 means save. The output file will be saved to the directory 'Code_plt_ms_figures/sim_data/' so that the corresponding figure can be produced.

- **Code_plt_ms_figures:**
Read in the data from 'Code_plt_ms_figures/sim_data/' and plot the figures: Figures 1-3 and Supp Figures 2-5. If saved, the figures will save to the folder 'Figures_ms_all.' From here, they were uploaded to the Overleaf document.

- **Manuscript_forCodeReview:** Folder where you can find the draft of the manuscript that corresponds to Figures in 'Figures_ms_all.'

- **Figures_ms_all:** Figure files in the manuscript for main Figures 1-3 and Supplemental Figures 1-5.


- **Data** Folder containing data shown and used: <br>
  - 'agedependent_data.xlsx' - clinical severity estimates from Davies et al. (2020); age distribution of Shanghai and contact matrix from Shanghai (Zhang et al. 2020). Used to parametrise age-dependent SEIR model. <br>
  - 'MMWR_Boehmeretal_medianages_4USregions_MayAug2020.xlsx' - median age data from MMWR Boehmer et al. (2020) <br>
  - 'Provisional_COVID19_DeathCounts_by_Sex_Age_and_Week.csv' - CDC provisional end of week death counts from COVID-19 US 2020 <br>


**References:** <br>
[1] Boehmer, Tegan K., et al. "Changing age distribution of the COVID-19 pandemic—United States, May–August 2020." Morbidity and Mortality Weekly Report 69.39 (2020): 1404. <br>
[2] Davies, Nicholas G., et al. "Age-dependent effects in the transmission and control of COVID-19 epidemics." Nature Medicine 26.8 (2020): 1205-1211. <br>
[3] Diekmann, Odo, J. A. P. Heesterbeek, and Michael G. Roberts. "The construction of next-generation matrices for compartmental epidemic models." Journal of the Royal Society Interface 7.47 (2010): 873-885. <br>
[4] Zhang, Juanjuan, et al. "Changes in contact patterns shape the dynamics of the COVID-19 outbreak in China." Science 368.6498 (2020): 1481-1486. <br>
