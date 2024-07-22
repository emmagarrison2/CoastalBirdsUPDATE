#######Objective of this script: refine list of all bird species (with UAI, MUTI, and/or UN scores) to COASTAL BIRDS ONLY 

library(ape)
library(geiger)
library(nlme)
library(effects)
library(ggplot2)
library(ggeffects)
library(sjPlot)
library(dplyr)
library(here)

######## 1st - load in "UAI_MUTI_UN_final.rds" file 
# it is all species with at least one of the three urban tolerance index scores 
#created in Step2 script 

AllBirds <- readRDS(here("Outputs", "UAI_MUTI_UN_final.rds"))

View(AllBirds)
nrow(AllBirds) #4435, which is correct! 


######### 2nd - Using Birds of the World (Cornell Lab or Ornithology), we classified each family represented by AllBirds dataframe as either 
#"Coastal" or "Not-Coastal", based on whether the family habitat page mentioned coastal areas. 
#The column "Coastal" categorizes these families. 
#"Yes" = Coastal family 
#"N/A" = Non-coastal family 

coastal_families <- read.csv(here("Notes", "FamilyNames.csv"))
View(coastal_families)
