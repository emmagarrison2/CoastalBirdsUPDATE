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

########Starting point - load in "UAI_MUTI_UN_final.rds" file 
# it is all species with at least one of the three urban tolerance index scores 
#created in Step2 script 

#insert EDITED csv here 
#AllBirds <- readRDS(here("Outputs", "UAI_MUTI_UN_final.rds"))

View(AllBirds)
nrow(AllBirds) #4435, which is correct! 

#had to alter the original UAI_MUTI_UN_final.csv from Step2 Script, because there were no 
#family names written down for Species_UN (that didn't have UAI or MUTI scores)


#combine rows Family_Jetz(from UAI, along with 80 species UN-only) and family.MUTI (from MUTI)
#let's call this combined family row "Family_All"

##############Round 1###################

### Using Birds of the World (Cornell Lab or Ornithology), we classified each family represented by AllBirds dataframe as either 
#"Coastal" or "Not-Coastal", based on whether the family habitat page mentioned coastal areas. 
#The column "Coastal" categorizes these families. 
#"Yes" = Coastal family 
#"N/A" = Non-coastal family 

coastal_families <- read.csv(here("Notes", "FamilyNames.csv"))
View(coastal_families)

### For all species from a "Coastal" (Yes) family, mark as a coastal species. 

#FamilyJetz from AllBirds with #BLFamilyLatin from coastal_families

