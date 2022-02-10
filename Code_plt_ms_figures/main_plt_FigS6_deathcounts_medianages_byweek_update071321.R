# Fig S1 - open and plot cdc surveillance data for
# A) COVID-19 deaths
# B) Median age of infection

# clear console screen
cat("\014");

# uncomment to install packages; only have to run these once
# install.packages("tidyverse");  # install package
# install.packages("zoo");
# install.packages("ggsci");
# install.packages("egg");

# load a bunch of packages first!
# library(zoo);
# library(tidyverse);
# library(lubridate);
# library(ggplot2);
# library(dplyr);
# library(ggsci);
# library(egg);

# set working directory
# need to change to "/XX/Code_plt_ms_figures" where XX is the location of the project
setwd("/Users/jharris387/Dropbox (Personal)/asymptomatic_spread_clean1/Code_plt_ms_figures");

# original data downloaded 04/05/21
# https://data.cdc.gov/Case-Surveillance/COVID-19-Case-Surveillance-Public-Use-Data/vbim-akqf
my_file <- "./../Data/Provisional_COVID19_DeathCounts_by_Sex_Age_and_Week.csv";

# read in the csv data file
us_data <- read_csv(my_file);

us_data <- rename(us_data, data_as_of=`Data as of`,state=`State`,MMWR_week=`MMWR Week`,end_week=`End Week`,sex = `Sex`,age_group =`Age Group`, total_deaths = `Total Deaths`, covid19_deaths = `COVID-19 Deaths`);

# now reduce to relevant variables
# group by the month data
df_groupby_week <- us_data %>% select(MMWR_week,end_week,age_group,covid19_deaths) %>%
  filter(age_group =="All Ages") %>%
  mutate(date=mdy(end_week)) %>%
  arrange(date) %>%
  group_by(date);

df_all_ages <- df_groupby_week %>% summarise(total_deaths_per_week = sum(covid19_deaths,na.rm = TRUE)/2);

# Deaths due to COVID-19 by week from 24 May 2020 to 19 Sept 2020: data are shifted by 21 days to estimate the shape of the incident cases from the beginning of May to the end of August.

# let's start at week 22 = 2020-05-30 up to week 37 = 2020-09-12
g1 <- ggplot(df_all_ages, aes(date, total_deaths_per_week, group = 1)) +
  geom_bar(fill="#FFDDAF", width=6, stat="identity",alpha=0.75)+
  geom_line(position = "identity",colour = "#1261A0", size = 0.75,alpha=0.9)+
  labs(x = "Month", y = "COVID-19 deaths (per week)")+
  scale_x_date(limits = as.Date(c("2020-05-24","2020-09-19")),breaks = as.Date(c("2020-06-13","2020-07-12","2020-08-13","2020-09-13")), labels =c("May","June","July","Aug"),expand = c(0, 0))+
  scale_y_continuous(limits = c(0,10000),expand = c(0, 0))+
  theme_bw()+
  theme(plot.tag = element_text(),
        plot.margin = unit(c(1,0.5,0.5,0.5), "lines"),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12))


# now load the median age data
# Median age of COVID-19 infection (positive RT-PCR) beginning of week from 1 May 2020 to 28 Aug 2020
# extracted from Beohmer et al.: https://www.cdc.gov/mmwr/volumes/69/wr/mm6939e1.htm
my_file <- "./../Data/MMWR_Boehmeretal_medianages_4USregions_MayAug2020.csv";

# read in the csv data file
us_medianage_data <- read_csv(my_file);  

g2 <- ggplot(us_medianage_data, aes(x=end_week, y=median_positive, color=Region)) +
  geom_line(size = 0.75,alpha=0.9)+
  labs(x = "Month", y = "Median Age")+
  scale_x_date(limits = as.Date(c("2020-05-01","2020-08-28")),breaks = as.Date(c("2020-05-22","2020-06-22","2020-07-22","2020-08-22")), labels =c("May","June","July","Aug"),expand = c(0, 0))+
  scale_y_continuous(limits = c(30,55),expand = c(0, 0))+
  scale_color_manual(values=c('#000000','#56B3E9','#E69F00','#009E73','#D55E00'),limits = c("US","Midwest","Northeast","South","West"))+
  theme_bw()+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "lines"),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12));

g3 <- ggarrange(g1, g2, nrow=2)
# ggsave("./../Figures_ms_all/supp/FigS1_medianage_changes_US.pdf", g3,height=11, width=8)
# export as .png in plot window so that the formmatted text size is saved
# note: added panel labels A and B later in adobe illustrator and saved as .eps 
