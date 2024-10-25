###Objective of this script is to TEMPORARILY test out .rds files 
# left behind in the Results Folder of this R Project (CoastalBirdsUPDATE). 
# if the file is determined to no longer be relevant to the R Project, check for 
# an updated replacement file, then delete. 


if(!require(here)){
  install.packages("here")
  require(here)
}
library(here)


if(!require(dplyr)){
  install.packages("dplyr")
  require(dplyr)
}
library(dplyr)


#Let's start with files I already moved to the "Potential Delete" Folder 

# 1 - BodyMass_MUTI.rds 
test <- readRDS(here("Results", "Potential Delete", "BodyMass_MUTI.rds"))
print(test)
#old version of body mass and MUTI plot. 
#this is not significant anymore. 


# 2 - develop_model_UN.rds
test <- readRDS(here("Results", "Potential Delete", "develop_model_UN.rds"))
print(test)
# old verion - new version (in FIGURES_results.R script)


# 3 - Developmental_MUTI.rds
test <- readRDS(here("Results", "Potential Delete", "Developmental_MUTI.rds"))
print(test)
# old verion - new version (in FIGURES_results.R script)


#KEEP # 4 - UAI_territorial.rds
test <- readRDS(here("Results", "Potential Delete", "UAI_territorial.rds"))
print(test)
#keep for now... still in FIGURES_results.R script, and is a good little guide to what this model looks like --> colors would, of course, need to change. 


# 5 - UN_developmental.rds
test <- readRDS(here("Results", "Potential Delete", "UN_developmental.rds"))
print(test)
# just plain wrong, delete 


# 6 - UN_high_nest.rds
test <- readRDS(here("Results", "Potential Delete", "UN_high_nest.rds"))
print(test)
# just plain wrong, delete 



#now let's look at files still in the results folder... 



# KEEP # 7 - BodyMass_UAI.rds
test <- readRDS(here("Results", "BodyMass_UAI.rds"))
print(test)
#looks accurate - in FIGURES_results.R script 


# KEEP # 8 - BroodValue_UAI_plot.rds
test <- readRDS(here("Results", "BroodValue_UAI_plot.rds"))
print(test)
#in FIGURES_results.R script, WITH outlier 


# KEEP # 9 - BroodValue_UAI_plot_filtered.rds
#in FIGURES_results.R script, WITH outlier 


# KEEP # 10 - clutch_UAI_plot.rds
#in FIGURES_results.R script


# KEEP # 11 - develop_UN_plot.rds
#in FIGURES_results.R script


# KEEP # 12 - life_history_figures.rds
#in FIGURES_results.R script


# KEEP # 13 - MUTI_developmental.rds
#in FIGURES_results.R script


# KEEP # 14 - MUTI_high_nest.rds
#in FIGURES_results.R scripts


# KEEP # 15 - MUTI_low_nest.rds
#in FIGURES_results.R scripts


# KEEP # 16 - Nest_Safety_UN_plot.rds
#in FIGURES_results.R scripts


# KEEP # 17 - NestHigh_UN_Plot.rds
#in FIGURES_results.R scripts


# 18 - "other density for UN"
#old script, DELETE 


# KEEP # 19 - UAI_Density_Comparison.rds
# just organized 


# KEEP # 20 - MUTI_Density_Comparison.rds
# just organized 


# KEEP # 21 - UN_Density_Comparison.rds
# just organized 


# KEEP # 22 - UAI_high_nest.rds
#in FIGURES_results.R scripts


# KEEP # 23 - UAI_low_nest.rds
#in FIGURES_results.R scripts


# KEEP # 24 - UN_develop_stacked_barchart.rds
#in FIGURES_results.R scripts


# 25 - UN_nest_high_stacked_barchart.rds
test <- readRDS(here("Results", "UN_nest_high_stacked_barchart.rds"))
print(test)
#old version, delete 