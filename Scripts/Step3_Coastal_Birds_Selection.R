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
library(tidyverse)

########Starting point - load in "UAI_MUTI_UN_final.rds" file 
# it is all species with at least one of the three urban tolerance index scores 

AllBirds <- readRDS(here("Outputs", "UAI_MUTI_UN_final.rds"))

View(AllBirds)
nrow(AllBirds) #4434, which is correct! 

# extract a data frame that contain a list of family names
unique_family_names <- AllBirds %>%
  distinct(Family_eBird) %>%
  tidyr::separate(., col = Family_eBird, # this splits the column Family_eBird into two columns
                  into = c("Family_Sci", "Family_English"), # the first column is called Family_Sci and the second Family_English
                  sep =  "^\\S*\\K\\s+") # split at the first space encountered

nrow(unique_family_names) # 205 families

View(unique_family_names)

# great - now lets download this df as a csv, and assess all Families using Birds of the World

write.csv(unique_family_names, here("Notes", "families_all.csv"))

##############Round 1###################

### Using Birds of the World (Cornell Lab or Ornithology), we classified each family represented by AllBirds dataframe as either 
#"Coastal" or "Not-Coastal", based on whether the family habitat page mentioned coastal areas. 
#The column "Coastal" categorizes these families. 
#"Yes" = Coastal family 
#"N/A" = Non-coastal family 


##################################################################
#################################################################
# EVERYTHING IN THIS SECTION CAN EVENTUALLY BE REMOVED
# SARAH'S SUGGESTION - ARCHIVE THIS SCRIPT (under a different file name) AND ALL THE FAMILY LIST CSV FILES THAT WE NO LONGER NEED
# I RECOMMEND YOU PUT ALL OF THEM IN A FOLDER OUTSIDE OF YOUR R PROJECT
# THAT WAY WE CAN ACCESS THEM IF WE NEED TO REVISIT THIS FOR SOME REASON 
# OTHERWISE THE STEPS SHOULD BE TO EXPORT families_all.csv AS SHOWN IN LINE 35
# ADD DESCRIPTION OF SEARCH PROCESS USING BIRDS OF THE WORLD
# AND THEN IMPORT THE CSV FILE CALLED families_all_coastal.csv
# FINALLY, JOIN THE COASTAL CLASSIFICATIONS TO AllBirds


# get a list of new families that have been added by cleaning up the joining process
#this section can be deleted eventually 


FamilyNames_old <- read.csv(here("Notes", "FamilyNames.csv"))
colnames(FamilyNames_old)
nrow(FamilyNames_old)#194
colnames(unique_family_names)
nrow(unique_family_names)#205

FamilyNames_old$Family_Sci <- FamilyNames_old$BLFamilyLatin
colnames(FamilyNames_old)

# using a left join
# joining the already investigated family names to the new list of family names (unique_family_names)
families_joined <- left_join(unique_family_names, FamilyNames_old, by="Family_Sci")
View(families_joined)
nrow(families_joined) #205


# get a list of already investigated families 
families_good <- families_joined %>%
  filter(!is.na(BLFamilyLatin))
nrow(families_good) # 163 families

# get a list of newly added families that need to be investigated
families_investigate <- families_joined %>%
  filter(is.na(BLFamilyLatin))
View(families_investigate) # 42 new families

#there are a few additional families, so let's look through Birds of the World and add their information in via editing csv! 
write.csv(families_investigate, here("Notes", "new_families_to_investigate.csv"))

# Emma looked these 42 families up on BoW
# import results of that process
families_investigated <- read.csv(here("Notes", "new_families_INVESTIGATED.csv"), header=T)

# we need to combine these 42 with the 163 from above
colnames(families_investigated)
colnames(families_good)

# first simplify both to keep them clean
f1_simple <- families_investigated %>% select(Family_Sci, Family_English, Coastal, Notes)
f2_simple <- families_good %>% select(Family_Sci, Family_English, Coastal, Notes)
f_all <- bind_rows(f1_simple, f2_simple)

# few last clean up steps
# convert Accipitridae, Cathartidae, and Falconiidae to be Coastal = Yes
f_all$Coastal[f_all$Notes == "INVESTIGATE"] <- "Yes" # change them to be coastal
f_all$Notes[f_all$Notes == "INVESTIGATE"] <- "" # remove text that says "INVESTIGATE"
# make all families that were not identified as Coastal have "No" in Coastal column
families_FINAL <- f_all %>% mutate(Coastal = if_else(Coastal == "Yes", Coastal, "No"))

write.csv(families_FINAL, here("Notes", "families_all_coastal.csv"))

##################################################################
#################################################################

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


#########SARAH- note: we will be marking Accipitriduae, Cathartidae, and Falconidae as "Yes" for coastal families. 


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

saveRDS(Coastal_Round_1_f, here("Outputs", "Coastal_Round_1.rds"))

####################################Round 2###################################
####
##
#in this round, we will sort through the common names of species that were marked as "Yes" for Urban Tolerance 
Round_1_yes <- Coastal_Round_1_f %>% 
  filter(Coastal == "Yes")

View(Round_1_yes)
nrow(Round_1_yes)
#825!  

write.csv(Round_1_yes, here("Notes", "Round_1_yes.csv"))

#now, look through the common names for these species... searching for key words that indicate NON-coastal habitats: "freshwater", "alpine", "upland", "lake", "river", 
#"mountain", "prairie", "highland", "forest", "desert" 

###### FINAL ROUND TO FIND ANY REMAINING COASTAL SPECIES #######
###### Before proceeding, do one last pass to identify any additional coastal species ######
# to do this, we will use the diet info columns that came from Wilman et al. 2014 (Elton traits) 

head(AllBirds)
# we need to search the original species list of 4434 species

# import elton traits (Wilman et al. 2014)
elton <- read.csv(here("Data", "elton.csv"), header=T)
head(elton)
colnames(elton) # look at column names of join 8 to identify columns that could be useful

# these ones seem useful:
unique(elton$Diet.Vfish) # percentage of diet that is fish. Filter to retain any species with > 0
unique(elton$Diet.5Cat) # Omnivore, VertFishScav, Invertebrate all seem potentially relevant
unique(elton$ForStrat.watbelowsurf) # percentage of time spent foraging below surf. Filter to retain any species with >0
unique(elton$ForStrat.wataroundsurf) # percentage of time spent foraging around surf. Filter to retain any species with >0
unique(elton$PelagicSpecialist) # Pelagic seabirds. Has values 0 or 1. Filter to retain species listed as 1


# join elton traits and AllBirds
AllBirds_elton <- elton %>%
  rename(Species_Jetz = Scientific) %>%
  select(Species_Jetz, Diet.Vfish, Diet.5Cat, ForStrat.watbelowsurf, # retain columns identified as useful above
         ForStrat.wataroundsurf, PelagicSpecialist) %>%
  left_join(AllBirds, ., by = "Species_Jetz")


# the diet category (Diet.5Cat) is the most vague and will likely pull many birds that are not coastal
# so, will implement an initial filter that the species must be one of the 3 diet categories AND...
# they must also fit one of the other criteria in the second filtering step
# note that the second set of filter requirements are OR statements 
# this means a bird should be retained if they have Diet.Vfish > 0 OR if they do any of their foraging below surf etc
# finally, we use a third filter to find all birds that meet the above requirements BUT are NOT currently classified as coastal
# we will want to investigate these species to see if they warrant inclusion 

coastaldiet <- AllBirds_elton %>% 
  filter(Diet.5Cat %in% c("Omnivore", "VertFishScav", "Invertebrate")) %>% # the first pass is to keep species with diet classifications of Omnivore, Invertebrate or VertFishScav 
  filter(Diet.Vfish > 0 |   # Now apply a second series of filtering requirements. First, keep any species with some fish in diet OR
           ForStrat.watbelowsurf > 0 | # keep any species that do any of their foraging below surf OR
           ForStrat.wataroundsurf >0 | # keep any species that do any of their foraging around surf OR
           PelagicSpecialist == 1)  # keep any species that are classified as Pelagic Specialists

nrow(coastaldiet) # 569 species
View(coastaldiet)

#### EMMA - 
# for the next step, you will want to take the data frame that contains all the coastal species you identified in all the previous steps
# you will need to use anti_join with the coastaldiet object to find which species were found by the diet classification method that are new
# eventually (after checking that all the species look reasonable) you'll want to drop all the elton trait columns and then use bind_rows to add these new species to the existing coastal list
coastal_new <- 