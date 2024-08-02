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

########################## Round 1 ###############################################

### Using Birds of the World (Cornell Lab or Ornithology), we classified each family represented by AllBirds dataframe as either 
#"Coastal" or "Not-Coastal", based on whether the family habitat page mentioned coastal areas. 
#The column "Coastal" categorizes these families. 
#"Yes" = Coastal family 
#"No" = Non-coastal family 


#We used Birds of the World (2022) to investigate each bird family represented by this dataset. If the family page mentioned use of
#coastal habitats (such as shorelines, beaches, estuaries, mangroves, etc.) then the family (and all representative species) was marked
#as "Yes" (coastal). If the Birds of the World family page did not mention use of coastal habitats or resources, then the all representative species of 
#said family were marked as "No" (not-coastal).

####### For all species from a "Coastal" (Yes) family, mark them as a coastal species. 
#For AllBirds --> column = Family_eBird
 
#For List_of_Families --> column = Family_Sci

#read in the csv that contains all families and their coastal status, as determined via investigation on Birds of the World (2022). 
List_of_Families <- read.csv(here("Notes", "families_all_coastal.csv")) %>% select(-X)
View(List_of_Families)
colnames(List_of_Families) 
colnames(AllBirds)
View(AllBirds)

#reformat AllBirds Family_eBird column so that it is compatable to join with List_of_Families 

AllBirds$Family_Sci <- word(AllBirds$Family_eBird, 1)
head(AllBirds)


#join them together 
Coastal_Round_1 <- left_join(AllBirds, List_of_Families, by = "Family_Sci")
View(Coastal_Round_1)
nrow(Coastal_Round_1) #4433 - correct number of rows in AllBirds, which was the left part of left_join 
 
#save rds of Coastal_Round_1, for quick recall 

saveRDS(Coastal_Round_1, here("Outputs", "Coastal_Round_1.rds"))


#################################### Round 2 ############################################
####
##
#in this round, we will sort through the common names of species that were marked as "Yes" for Coastal 

Round_1_yes <- Coastal_Round_1 %>% 
  filter(Coastal == "Yes")

View(Round_1_yes)
nrow(Round_1_yes)
#823!  

###
# generate a list of descriptors in the common names of the 4434 species in AllBirds
# we will manually filter this list to identify words in the common names that relate to habitat and/or diet

twowordnames <- AllBirds %>%
  select(CommonName_eBird) %>%
  mutate(space = str_count(CommonName_eBird, " ")) %>% # identify how many words each species has in its names using number of spaces
  filter(space == 1) %>% # keep only species with names that contain two words (a single space)
  separate_wider_delim(CommonName_eBird, delim = " ", names = c("descriptor1", "species")) 
# some of these still have hyphenated terms in the species column that could be descriptive

threewordnames <- AllBirds %>%
  select(CommonName_eBird) %>%
  mutate(space = str_count(CommonName_eBird, " ")) %>% # identify how many words each species has in its names using number of spaces
  filter(space == 2) %>% # keep only species with names that contain two words (two spaces)
  separate_wider_delim(CommonName_eBird, delim = " ", names = c("descriptor1", "descriptor2", "species") )
# some of these still have hyphenated terms in the species column that could be descriptive

hyphennames <- data.frame(species = c(twowordnames$species, threewordnames$species)) %>%
  mutate(hyphen = str_count(species, "-")) %>% # identify species with hyphenated names
  filter(hyphen == 1) %>% # keep only species with names that contain a hyphen
  separate_wider_delim(species, delim = "-", names = c("descriptor1", "species") )

# put all the descriptive words from common names together into a single data frame and remove duplicates
alldescriptors <- data.frame( # make a data frame
  descriptors = sort( # arrange in alphabetical order
    c(twowordnames$descriptor1, threewordnames$descriptor1, threewordnames$descriptor2, hyphennames$descriptor1))) %>% 
  distinct() # remove duplicated descriptors

# export as csv 
# manually mark all terms that are associated with:
# a habitat (but not a specific place)
# a diet
# an Ocean or Sea (e.g., Pacific)
# or indicated a species association or interaction between the bird and another species (e.g., plant, insect, mammal)
write.csv(alldescriptors, here("Notes", "descriptors.csv"))

# import the marked up list
# use column coastal_Y.N.M to sort
descriptors_marked <- read.csv(here("Notes", "descriptors_marked.csv"), header=T)
head(descriptors_marked)
View(descriptors_marked)

# first focus on non-coastal terms
noncoastal_terms <- descriptors_marked %>% 
  filter(coastal_Y.N.M =="N") %>% # keep all terms that are clearly non-coastal
  select(descriptors) # select the column of descriptive terms

# convert to a character vector 
NCterms_vec <- as.vector(noncoastal_terms[,1])
class(NCterms_vec)

# convert complete list of common names from Round_1_yes to a character vector 
head(Round_1_yes)
commname_vec <- as.vector(Round_1_yes[,5]) # common names are in 5th column

# now, look through commname_vec which contains the species in Round_1_yes using the terms stored in NCterms_vec
# we are searching for key words that indicate NON-coastal habitats, diets or associations with other species that are non-coastal
# e.g., "freshwater", "alpine", "upland", "lake", "river", 
#"mountain", "prairie", "highland", "forest", "desert" 

R1_noncoastal <- map(NCterms_vec, str_detect, string = commname_vec) %>%
  reduce(`|`) %>% 
  magrittr::extract(commname_vec, .) %>%
  tibble()

colnames(R1_noncoastal) <- "CommonName_eBird"

print(R1_noncoastal, n=Inf)

# export data frame, look up species and mark ones to keep and remove
write.csv(R1_noncoastal, here("Notes", "Round_1_noncoastal.csv"))

# manually look up species and determine whether they should be marked as coastal or not
# read in edited csv
# use it to update coastal species list

# read in edited csv 
R1_noncoastal_edited <- read.csv(here("Notes", "Round_1_noncoastal_edited.csv"), header=T) %>%
  select(-X)

head(R1_noncoastal_edited)
colnames(R1_noncoastal_edited)
nrow(R1_noncoastal_edited)

# combine with Coastal_Round_1 to make the updates
# this takes a few steps

# remind ourselves how many species are in the data
nrow(Coastal_Round_1) # 4433

# get all the species that were not modified in any way in the past step (any species not in R1_noncoastal_edited)
Coastal_Round_1_edit1 <- anti_join(Coastal_Round_1, R1_noncoastal_edited, by="CommonName_eBird")
# note: you need to specify to use CommonName_eBird for this join for it to work correctly

# does the number of rows in edit1 equal the number of rows in Round_1 minus the number of rows in noncoastal_edited?
# this should be "TRUE"
nrow(Coastal_Round_1_edit1) == nrow(Coastal_Round_1) - nrow(R1_noncoastal_edited)

# update the coastal classification for all species in R1_noncoastal_edited
Coastal_Round_1_edit2 <- Coastal_Round_1 %>% 
  select(-Coastal, -Notes) %>% # removes these columns because they will be updated with new versions in the next step
  left_join(R1_noncoastal_edited, ., by="CommonName_eBird")
nrow(Coastal_Round_1_edit2) # should be the same number as in R1_noncoastal_edited

# bind everything back together
Coastal_Round_2 <- bind_rows(Coastal_Round_1_edit1, Coastal_Round_1_edit2)  

# make sure all the species are still present
nrow(Coastal_Round_2) == nrow(Coastal_Round_1)

# View Coastal_Round_2 for double-checking 
View(Coastal_Round_2)

# looks good! 


#################################### Round 3 ###################################
####
##
#in this round, we will sort through the common names of species that were marked as "Yes" for Urban Tolerance 

Round_1_no <- Coastal_Round_1 %>% 
  filter(Coastal == "No")

View(Round_1_no)
nrow(Round_1_no)
#3610!  

# # # # # # # # # # # # # # # # # # 

# now, develop list of coastal descriptor terms to search the Common Names of species marked as not coastal

coastal_terms <- descriptors_marked %>% 
  filter(coastal_Y.N.M %in% c("M", "Y")) %>% # keep all terms that are clearly coastal (Y) and possibly coastal (M)
  select(descriptors) # select the column of descriptive terms

# convert to a character vector 
C.terms_vec <- as.vector(coastal_terms[,1])
class(C.terms_vec)

# convert complete list of common names from Round_1_yes to a character vector 
head(Round_1_no)
commname_vec_2 <- as.vector(Round_1_no[,5]) # common names are in 5th column

# now, edit this Round_1_no.csv file (containing all species from families marked as "No" in Round 1) -> 
#search through common names for coastal-identifier words stored in CTerms_vec (i.e. "Coastal", "Sea", "Tide", "Beach", "Mangrove", etc.) 
#for all flagged species, look on Birds of the World (2022) species page for mentions of coastal habitat/resource use. 

R1_coastal <- map(C.terms_vec, str_detect, string = commname_vec_2) %>%
  reduce(`|`) %>% 
  magrittr::extract(commname_vec_2, .) %>%
  tibble()

R1_coastal
nrow(R1_coastal) # 221

View(R1_coastal)
# export data frame, look up species and mark ones to keep and remove

write.csv(R1_coastal, here("Notes", "Round_1_Coastal.csv"))

# # # # # # # # # # # # # # # # # # 


# upload edited csv 
# SEE THE CODE FOR LINES 188 to 220 FOR HOW TO BEST JOIN THE EDITS
# YOU SHOULD BE ABLE TO MIMIC THAT PROCESS HERE
# YOU'LL WANT TO MAKE THE JOIN TO Coastal_Round_2


R1_coastal_edited <- read.csv(here("Notes", "Round_1_Coastal_edited.csv")) %>% select(-X, -X.1)

#View to visually check 
View(R1_coastal_edited)
colnames(R1_coastal_edited)
nrow(R1_coastal_edited)



# combine with Coastal_Round_2 to make the updates
# this takes a few steps

# remind ourselves how many species are in the data
nrow(Coastal_Round_2) # 4433

# get all the species that were not modified in any way in the past step (any species not in R1_coastal_edited)
Coastal_Round_2_edit1 <- anti_join(Coastal_Round_2, R1_coastal_edited, by="CommonName_eBird")
# note: you need to specify to use CommonName_eBird for this join for it to work correctly

# does the number of rows in edit1 equal the number of rows in Round_2 minus the number of rows in coastal_edited?
# this should be "TRUE"
nrow(Coastal_Round_2_edit1) == nrow(Coastal_Round_2) - nrow(R1_coastal_edited)
#TRUE
nrow(R1_coastal_edited)

# update the coastal classification for all species in R1_coastal_edited
Coastal_Round_2_edit2 <- Coastal_Round_2 %>% 
  select(-Coastal, -Notes) %>% # removes these columns because they will be updated with new versions in the next step
  left_join(R1_coastal_edited, ., by="CommonName_eBird")
nrow(Coastal_Round_2_edit2) # should be the same number as in R1_coastal_edited

# bind everything back together
Coastal_Round_3 <- bind_rows(Coastal_Round_2_edit1, Coastal_Round_2_edit2)  

# make sure all the species are still present
nrow(Coastal_Round_3) == nrow(Coastal_Round_2)

# View Coastal_Round_2 for double-checking 
View(Coastal_Round_3)

#looks good! 



########################### Final Step ###########################################################
# one last pass to identify any additional coastal species or species that are not actually coastal 
# we will use the diet info columns that came from Wilman et al. 2014 (Elton traits) 

head(AllBirds)
# we need to search the original species list of 4433 species

# import elton traits (Wilman et al. 2014)
elton <- read.csv(here("Data", "elton.csv"), header=T)
head(elton)
colnames(elton) # look at column names of join 8 to identify columns that could be useful

# these ones seem useful:
unique(elton$Diet.Vfish) # percentage of diet that is fish. Filter to retain any species with > 0
unique(elton$ForStrat.watbelowsurf) # percentage of time spent foraging below surf. Filter to retain any species with >0
unique(elton$ForStrat.wataroundsurf) # percentage of time spent foraging around surf. Filter to retain any species with >0
unique(elton$PelagicSpecialist) # Pelagic seabirds. Has values 0 or 1. Filter to retain species listed as 1


# join elton traits and AllBirds
AllBirds_elton <- elton %>%
  rename(Species_Jetz = Scientific) %>%
  select(Species_Jetz, Diet.Vfish, Diet.5Cat, ForStrat.watbelowsurf, # retain columns identified as useful above
         ForStrat.wataroundsurf, PelagicSpecialist) %>%
  left_join(AllBirds, ., by = "Species_Jetz")


# note that the filter requirements are OR statements 
# this means a bird should be retained if they have Diet.Vfish > 0 OR 
# if they do any of their foraging below surf OR
# if they forage around water OR
# if they are a pelagic specialist
# finally, we use a third filter to find all birds that meet the above requirements BUT are NOT currently classified as coastal
# we will want to investigate these species to see if they warrant inclusion 

coastaldiet <- AllBirds_elton %>% 
  filter(Diet.Vfish > 0 |   # Now apply a of filtering requirements. First, keep any species with some fish in diet OR
           ForStrat.watbelowsurf > 0 | # keep any species that do any of their foraging below surf OR
           ForStrat.wataroundsurf >0 | # keep any species that do any of their foraging around surf OR
           PelagicSpecialist == 1)  # keep any species that are classified as Pelagic Specialists

nrow(coastaldiet) # 688 species
View(coastaldiet)

# Are any birds in the coastaldiet list currently marked as non-coastal?
head(Coastal_Round_3)
R3_No <- Coastal_Round_3 %>% 
  filter(Coastal == "No")
R3_noncoastal <- inner_join(R3_No, coastaldiet)

nrow(R3_noncoastal)
#58 species --> these are species tht have been marked as "No" over the past 3 rounds, but their diet traits suggest a coastal association. 

# Export this list, look up all species and check whether they should be coastal
write.csv(R3_noncoastal, here("Notes", "Round_3_noncoastal.csv"))

#######################################

# Are any birds that are marked as coastal that don't have a diet with fish or forage in/around water?
# NOTE: species can still be coastal and not eat fish or foraging in/around water
# so this step may not result in the removal of any species from the coastal list
# we are using this as a final check to flag any species that may need a second look
R3_Yes <- Coastal_Round_3 %>% 
  filter(Coastal == "Yes") 

R3_coastal <- anti_join(R3_Yes, coastaldiet) %>%
  filter(Notes=="") %>% # only keep species where there is no Note. Any species with a Note has already been investigated
  arrange(Family_eBird)

nrow(R3_coastal)
#170
# many of these appear to be coastal
# we will export the list, double check their classification as coastal and reimport the list with any needed changes

# Export
write.csv(R3_coastal, here("Notes", "Round_3_coastal.csv"))






#edit files Round_3_coastal.csv and Round_3_noncoastal.csv by checking all species on BOTW 



# Import edited files
R3_coastal_edited <- read.csv(here("Notes", "Round_3_coastal_edited2.csv"), header = T) %>% select(-X)

R3_noncoastal_edited <- read.csv(here("Notes", "Round_3_noncoastal_edited2.csv"), header = T) %>% select(-X)


# combine them into one object that contains all possible changes
R3_edits <- bind_rows(R3_coastal_edited, R3_noncoastal_edited)

# get all the species that were not modified in any way in the past step 
Coastal_Round_4_edit1 <- anti_join(Coastal_Round_3, R3_edits, by="CommonName_eBird")
# note: you need to specify to use CommonName_eBird for this join for it to work correctly

# does the number of rows in edit1 equal the number of rows in Round_3 minus the number of rows in R3_edits ?
# this should be "TRUE"
nrow(Coastal_Round_4_edit1) == nrow(Coastal_Round_3) - nrow(R3_edits)
#TRUE

colnames(R3_edits)
colnames(Coastal_Round_3)
# update the coastal classification for all species in R3_edits
Coastal_Round_4_edit2 <- Coastal_Round_3 %>% 
  select(-Coastal, -Notes) %>% # removes these columns because they will be updated with new versions in the next step
  left_join(R3_edits, ., by="CommonName_eBird")

nrow(Coastal_Round_4_edit2) # 228 - should be the same number as in R3_edits
nrow(R3_edits) #228!


# bind everything back together
Coastal_Final <- bind_rows(Coastal_Round_4_edit1, Coastal_Round_4_edit2)  

# make sure all the species are still present
nrow(Coastal_Round_3) == nrow(Coastal_Final)
#TRUE 

# View Coastal_Final for double-checking 
View(Coastal_Final)

#looks good! 

#how many coastal species do we have? 

nrow(Coastal_Final)
colnames(Coastal_Final)
#4433, just like it should be 

Coastal_Birds <- Coastal_Final %>% 
  filter(Coastal == "Yes")

nrow(Coastal_Birds)
#826 birds! Wow! 


#View to confirm  

View(Coastal_Birds)


#save as csv and rds in data folder! 

write.csv(Coastal_Birds, here("Data", "Coastal_Birds_List.csv"))
saveRDS(Coastal_Birds, here("Data", "Coastal_Birds_list.rds"))
