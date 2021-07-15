# Impact of time scales on proportion of asymptomatic transmission and incidence
Jeremy Harris and Sang Woo Park and collaborators, June 17, 2021

 %% updated 06/25/21 by jdh %%

**Code for:** JD Harris and SW Park et al. "Differences in time scales coupled with assortative mixing can lead tochanges in the proportion of asymptomatic transmission and incidenceover time." This repository contains all the Matlab codes necessary for the running the SEIR models of asymptomatic transmission: (1) fixing proportion of asymptomatic incidence (no assortative mixing); (2) assortative mixing based on infection type, either asymptomatic or symptomatic; (3) assortative mixing by age using estimates of contact rates, susceptibility, and clinical severity by age.

A preprint of the manuscript can be found on BioRxiv: [link](XXX)

This code is archived on Zenodo: [![DOI](XXX)

**Instructions:** all codes are written in MATLAB. To plot the figures or do the analyses, you need download the project at https://github.com/Jeremy-D-Harris/XXX

**Folder descriptions:**

- **Manuscript_forCodeReview:** Folder where you can find an updated draft version of the manuscript.

- **Figures_ms_all** Figure files in the manuscript, including main Figures 1-3 and Supp Figures 1-13. should be raw outputs????

- **Code_sims:** All code to produce the figures: three subfolder with one for each of the models

- **Code_plot_ms_figures:** Scripts used to plot the figures in the manuscript, Figures 1-3 and Supp Figures 2, 4-6, 8-10, 12-13. Takes the data from sim_data subfolder within each of the subdirectories that are labeled by the model.

- **Data** Folder containing all the data (e.g., age distribution of Wuhan), empirical estimates (e.g. contact rates by age), and scripts to plot those data????



**References:**
- [1] Diekmann, Odo, J. A. P. Heesterbeek, and Michael G. Roberts. "The construction of next-generation matrices for compartmental epidemic models." Journal of the Royal Society Interface 7.47 (2010): 873-885.
