# Individual-level differences in symptomatic and asymptomatictransmission shape population-level dynamics of the SARS-CoV-2outbreak
Jeremy Harris and Sang Woo Park and collaborators, January 2022

 -- updated 01/24/22 by jdh --

**Code for:** JD Harris and SW Park et al. "Individual-level differences in symptomatic and asymptomatictransmission shape population-level dynamics of the SARS-CoV-2outbreak" This repository contains all the Matlab codes for running the three SEIR-like models of asymptomatic transmission: <br>
  (1) fixing proportion of asymptomatic incidence (no assortative mixing); <br>
  (2) assortative mixing based on infection type, either asymptomatic or symptomatic; <br>
  (3) assortative mixing by age using estimates of contact rates, susceptibility, and clinical severity by age.

A preprint of the manuscript can be found on BioRxiv: [link](XXX) -- not yet!!

This code is archived on Zenodo: [DOI](XXX) -- not yet!!

**Instructions:** With the exception of an R-script to plot CDC surveillance data on mortality and median age of infection (May-Aug 2020), MATLAB was used to parametrise models, run simulations, and plot figures. To do these analyses, download the project at https://github.com/Jeremy-D-Harris/Asymptomatic_Spread_COVID.git

Once the project is downloaded, navigate to the subdirectory 'Code_sims' to produce the simulation data. Once you produce the simulation, you can plot the figures in the manuscript by running the appropriate function in 'Code_plt_ms_figures.' See below for subfolder descriptions. -- actually you should be able to plot figures, since simulation data is already produced.

**Folder descriptions:**
Folders are
- **Code_sims:** All code to produce figures; subfolders organized by model with models (1)-(3) described above: <br>
  (1) 'fixedpropasymp_code' <br>
  (2) 'assortmixing_code' <br>
  (3) 'agedep_code'

Within each of these folders, you'll find main files that simulate and parametrise the model. There are several choices for the user at the top of main file scripts. For instance, in 'main_sim_assortmixing_SEIR_twodiseases_sameR0s_update011722.m' the the first choice for the user is to save the simulation data using the variable save_ans: 0 means don't save and 1 means save. The output file will be saved to the directory 'Code_plt_ms_figures/sim_data/' so that the corresponding figure can be produced.

- **Code_plt_ms_figures:**
 Read in the data from 'Code_plt_ms_figures/sim_data/' and plot the figures: Figures 1-3 and Supp Figures 2-5. If saved, the figures will save to the folder 'Figures_ms_all.' From here, they are uploaded to the Overleaf document.

- **Manuscript_forCodeReview:** Folder where you can find the draft of the manuscript that corresponds to Figures in 'Figures_ms_all.'

- **Figures_ms_all** Figure files in the manuscript, including main Figures 1-3 and Supp Figures 1-13.


- **Data** Folder containing data shown and used: <br>
  - 'agedependent_data.xlsx' includes estimates from Davies et al., age distribution of Wuhan, and contact matrix from Wuhan (taken from Zhang et al.). (Used to parametrise age-dependent SEIR model.) <br>
  - 'MMWR_Boehmeretal_medianages_4USregions_MayAug2020.xlsx' includes median age data from MMWR Boehmer et al. (2020) <br>
  - 'Provisional_COVID19_DeathCounts_by_Sex_Age_and_Week.csv' includes CDC provisional end of week death counts from
COVID-19 US 2020


**References:**
- [1] Boehmer, Tegan K., et al. "Changing age distribution of the COVID-19 pandemic—United States, May–August 2020." Morbidity and Mortality Weekly Report 69.39 (2020): 1404.
- [2] Davies, Nicholas G., et al. "Age-dependent effects in the transmission and control of COVID-19 epidemics." Nature medicine 26.8 (2020): 1205-1211.
- [3] Diekmann, Odo, J. A. P. Heesterbeek, and Michael G. Roberts. "The construction of next-generation matrices for compartmental epidemic models." Journal of the Royal Society Interface 7.47 (2010): 873-885.
- [4] Zhang, Juanjuan, et al. "Changes in contact patterns shape the dynamics of the COVID-19 outbreak in China." Science 368.6498 (2020): 1481-1486.
