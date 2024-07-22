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

#great - now lets download this df as a csv, and assess all Families using Birds of the World

write.csv(df_unique_family_names, here("Notes", "families_all.csv"))

##############Round 1###################

### Using Birds of the World (Cornell Lab or Ornithology), we classified each family represented by AllBirds dataframe as either 
#"Coastal" or "Not-Coastal", based on whether the family habitat page mentioned coastal areas. 
#The column "Coastal" categorizes these families. 
#"Yes" = Coastal family 
#"N/A" = Non-coastal family 


#go off of our previous FamilyNames.csv (FamilyNames_old) (we want to double check that no additional families were added via this UPDATE process... and we don't want to redo work)
#this section can be deleted eventually 
Familynames_new <- read.csv (here("Notes", "families_all_edited.csv"))
colnames(Familynames_new)

FamilyNames_old <- read.csv(here("Notes", "FamilyNames.csv"))
colnames(FamilyNames_old)
#View(coastal_families)
FamilyNames_old$Family <- FamilyNames_old$BLFamilyLatin
colnames(FamilyNames_old)


#full join 
families_joined <- full_join(Familynames_new, FamilyNames_old, by="Family")
View(families_joined)
#there are a few additional families, so let's look through Birds of the World and add their information in via editing csv! 
write.csv(families_joined, here("Notes", "families_joined.csv"))

#upload the edited version - 
Family_List <- read.csv(here("Notes", "Family_List_Coastal.csv"))
View(Family_List)
colnames(Family_List)
Family_List$Family_all <- Family_List$Family


####### For all species from a "Coastal" (Yes) family, mark them as a coastal species. 
#For AllBirds.r --> column = Family_all
#For Family_List --> column = Family 

#join them together 
Coastal_Families <- full_join(AllBirds.r, Family_List, by = "Family_all")
View(Coastal_Families)

#that worked 

#save rds of Coastal_Families, for quick recall 

saveRDS(Coastal_Families, here("Outputs", "Coastal_Families.rds"))




##############Round 1.5###################


#in this step, we will filter through the species from special-case families (such as Accipitridae) to see if we should put any of these species down as coastal. 

#I'll need to go off of our old list... cut down a csv of all species represented by old list that are from the following families... 

##Accipitridae 
##Cathartidae 
##Falconiidae 

Coastal_Fam <- readRDS(here("Outputs", "Coastal_Families.rds"))
colnames(Coastal_Fam)
refine_1 <- Coastal_Fam %>% select(Family_all == "Accipitridae", "Cathartidae", "Falconidae" )
