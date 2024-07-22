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
AllBirds <- read.csv(here("Outputs", "UAI_MUTI_UN_final_edited.csv"))

View(AllBirds)
nrow(AllBirds) #4435, which is correct! 

#had to alter the original UAI_MUTI_UN_final.csv from Step2 Script, because there were no 
#family names written down for Species_UN (that didn't have UAI or MUTI scores)


#combine rows Family_Jetz(from UAI) family.MUTI (from MUTI), and Family_UN_only (from UN)
#let's call this combined family row "Family_all"

AllBirds.r <- AllBirds %>% mutate (Family_all = coalesce(Family_Jetz, family.MUTI, Family_UN_only))
View(AllBirds.r)
colnames(AllBirds.r)

unique_family_names <- unique(AllBirds.r$Family_all)
df_unique_family_names <- data.frame(Family = unique_family_names)
#View(df_unique_family_names)

#great - now lets download this df as a csv, and assess all Families using Birds of the World

##############Round 1###################

### Using Birds of the World (Cornell Lab or Ornithology), we classified each family represented by AllBirds dataframe as either 
#"Coastal" or "Not-Coastal", based on whether the family habitat page mentioned coastal areas. 
#The column "Coastal" categorizes these families. 
#"Yes" = Coastal family 
#"N/A" = Non-coastal family 

coastal_families <- read.csv(here("Notes", "FamilyNames.csv"))
colnames(coastal_families)
#View(coastal_families)
coastal_families$Family_all <- coastal_families$Family_all

### For all species from a "Coastal" (Yes) family, mark as a coastal species. 



#Family_all from AllBirds.r with #BLFamilyLatin from coastal_families

