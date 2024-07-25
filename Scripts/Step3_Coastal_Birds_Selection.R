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
unique(Coastal_Fam$Family_all)
refine_1 <- Coastal_Fam %>% filter(Family_all %in% c("Accipitridae", "Cathartidae", "Falconidae"))
View(refine_1)
colnames(refine_1)
nrow(refine_1) #163

#great, now let's join with old coastal-species list, as I've probably already investigated most of these Accipitridae, Cathartidae, and Falconidae species 

#the Species_Jetz column contains a scientific name for all species in df "refine_1"
#duplicate column and name it "Species" 
refine_1$Species <- refine_1$Species_Jetz
#make sure Species column is formatted correctly for join 
refine_1$Species <- gsub(" ", "_", refine_1$Species)
View(refine_1)

#read in old coastal species list 
Old_List <- read.csv(here("Notes", "CoastalSpecies_11April24.csv"))
View(Old_List)
#instead of an "NA", write in "No" for CoastalSp column for all non-coastal species 
Old_List$CoastalSp <- ifelse(Old_List$CoastalSp =="Yes", "Yes", "No")
View(Old_List)
#View(refine_1)
#join with "Species" column from "CoastalSpecies_11April24.csv" (the old coastal species list)

Raptor_Search <- left_join(refine_1, Old_List, by = "Species")
View(Raptor_Search)


#great, now let's investigate all the species with NA for column "CoastalSp" -- these are species that were added to our list in the recent refined-join 

#to investigate these species on BOTW, let's export a csv 

write.csv(Raptor_Search, here("Notes", "raptors_to_BOTW.csv"))


#after going through list of NA species, and investigating on BOTW to assign as Coastal (Yes) or Non-coastal (No), read this edited csv in 
Coastal_raptors <- read.csv (here("Notes", "Coastal_raptors.csv"))
View(Coastal_raptors)
#correct column is "Coastal" 

colnames(Coastal_Families)
Coastal_Families$Species <- Coastal_Families$Species_Jetz
Coastal_Families$Species <- gsub(" ", "_", Coastal_Families$Species)

#join with Coastal_Families by column "Species_Jetz", as this was the column we used to join refine_1 with Old_List 
Coastal_Round_1 <- Coastal_Families %>%
  left_join(Coastal_raptors, by = "Species", suffix = c(".fam", ".rap"))
View(Coastal_Round_1)  

#cool now let's create a new column (Coastal) to coalesce Coastal.fam (families) and Coastal.rap (select raptors)
Coastal_Round_1.5 <- Coastal_Round_1 %>% 
  mutate(Coastal = coalesce(Coastal.rap, Coastal.fam))
View(Coastal_Round_1.5)

#simplify number of columns, it's getting overwhelming 
colnames(Coastal_Round_1.5)

Coastal_Round_1_f <- Coastal_Round_1.5 %>% 
  rename(Family_all = Family_all.fam) %>%
  select(MUTIscore, SciName_MUTI, CommonName_MUTI, Species_eBird, Species_Jetz, Family_Jetz, CommonName_UAI, aveUAI, Species_BirdLife, Species_UN, Urban, family.MUTI, Family_all, Species, English, Coastal)
View(Coastal_Round_1_f)



####################################Round 2###################################
####
##
#in this round, we will sort through the common names of species that were marked as "Yes" for Urban Tolerance 

