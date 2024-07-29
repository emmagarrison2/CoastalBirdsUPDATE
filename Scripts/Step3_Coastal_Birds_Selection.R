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

# extract a data frame that contains a list of family names
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
#"No" = Non-coastal family 


#We used Birds of the World (2022) to investigate each bird family represented by this dataset. If the family page mentioned use of
#coastal habitats (such as shorelines, beaches, estuaries, mangroves, etc.) then the family (and all representative species) was marked
#as "Yes" (coastal). If the Birds of the World family page did not mention use of coastal habitats or resources, then the all representative species of 
#said family were marked as "No" (not-coastal).

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



##################################################################
#################################################################


####### For all species from a "Coastal" (Yes) family, mark them as a coastal species. 
#For AllBirds --> column = Family_eBird
 
#For List_of_Families --> column = Family_Sci

#read in the csv that contains all families and their coastal status, as determined via investigation on Birds of the World (2022). 

List_of_Families <- read.csv(here("Notes", "families_all_coastal.csv"))
View(List_of_Families)
colnames(List_of_Families) 
colnames(AllBirds)
View(AllBirds)

#reformat AllBirds Family_eBird column so that it is compatable to join with List_of_Families 

AllBirds$Family_Sci <- word(AllBirds$Family_eBird, 1)
View(AllBirds)


#join them together 
Coastal_Round_1 <- left_join(AllBirds, List_of_Families, by = "Family_Sci")
View(Coastal_Round_1)
nrow(Coastal_Round_1)#4433 - correct number of rows in AllBirds, which was the left part of left_join 
 
#save rds of Coastal_Round_1, for quick recall 

saveRDS(Coastal_Round_1, here("Outputs", "Coastal_Round_1.rds"))


####################################Round 2###################################
####
##
#in this round, we will sort through the common names of species that were marked as "Yes" for Urban Tolerance 


Round_1_yes <- Coastal_Round_1 %>% 
  filter(Coastal == "Yes")

View(Round_1_yes)
nrow(Round_1_yes)
#823!  

write.csv(Round_1_yes, here("Notes", "Round_1_yes.csv"))

#now, look through the common names for these species... searching for key words that indicate NON-coastal habitats: "freshwater", "alpine", "upland", "lake", "river", 
#"mountain", "prairie", "highland", "forest", "desert" 


#read in edited csv 

edits_for_round_2 <- read.csv(here("Notes", "Round_1_yes_edited.csv"))

colnames(edits_for_round_2)


Coastal_Round_2 <- Coastal_Round_1 %>%
  left_join(edits_for_round_2, by = "Species_eBird", suffix = c("", ".r2")) %>%
  mutate(
    Coastal = ifelse(!is.na(Coastal.r2), Coastal.r2, Coastal), 
    Notes = ifelse(!is.na(Notes.r2), Notes.r2, Notes)
  ) %>%
  select(-ends_with(".r2"), -"X", -"X.1")

#View Coastal_Round_2 for double-checking 
View(Coastal_Round_2)


#looks good! 


####################################Round 2###################################
####
##
#in this round, we will sort through the common names of species that were marked as "Yes" for Urban Tolerance 


Round_1_no <- Coastal_Round_1 %>% 
  filter(Coastal == "No")

View(Round_1_no)
nrow(Round_1_no)
#3610!  

write.csv(Round_1_no, here("Notes", "Round_1_no.csv"))

#now, edit this Round_1_no.csv file (containing all species from families marked as "No" in Round 1) -> 
#search through common names for coastal-identifier words: "Coast", "Coastal", Sea", "Tide", "Intertidal", "Beach", "Mangrove",
#"Ocean", "Barnacle", "Kelp", "Fish", "Crab", "Bay", "Cove", "Cape", "Estuary", "Lagoon", and "Reef" 
#for all flagged species, look on Birds of the World (2022) species page for mentions of coastal habitat/resource use. 









###### FINAL ROUND TO FIND ANY REMAINING COASTAL SPECIES #######
###### Before proceeding, do one last pass to identify any additional coastal species ######
# to do this, we will use the diet info columns that came from Wilman et al. 2014 (Elton traits) 

head(AllBirds)
# we need to search the original species list of 4433 species

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