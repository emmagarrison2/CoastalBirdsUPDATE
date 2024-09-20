####################The objective of this script is to join all of our current indexes
#(MUTI, UAI, and UN) --> after we successfully join and get all scientific names to line up, 
#then we can run this file through our "retracing_coastal_species_sorting.r" script! 


library(ape)
library(geiger)
library(nlme)
library(effects)
library(ggplot2)
library(ggeffects)
library(sjPlot)
library(dplyr)
library(here)

if(!require(here)){
  install.packages("here")
  require(here)
}


here::here()

#install.packages("gitcreds")
library(gitcreds)

UAI <- readRDS (here("Outputs", "UAI_eBirdJetzBirdLife_final.rds"))
View(UAI)
colnames(UAI)
nrow(UAI)
#4347
min(UAI$aveUAI) #0
max(UAI$aveUAI) #3.972

#UAI Naming Scheme Methods --> weird combo of eBird and other naming schemes --> Sarah Investigated on UAINamestoSpecies.R

#MUTI Naming Scheme Methods --> seems to use eBird scheme --> uses eBird 6 letter codes, and AOS species checklist for taxonomic information. 
####There may be a few differences, but for the most part, eBird/Clements naming scheme and AOS agree. 

#UN Naming Scheme Methods --> ?? unstated naming scheme for bird species. 


#Plan --> try to join UAI and MUTI first... then try with UN (might be more difficult?)

#Load in MUTI! 

MUTI <- read.csv (here("Data", "Fanelli_Urban_Tolerance.csv"), header = T)
colnames(MUTI)
View(MUTI)
MUTI.r <- MUTI %>% rename(MUTIscore = PC1_scores) 
MUTI.f <- MUTI.r %>% select(MUTIscore, scientific_name, common_name, family) %>%
  rename(SciName_MUTI = scientific_name, CommonName_MUTI=common_name, Family_MUTI = family)

colnames(MUTI.f)
View(MUTI.f)
nrow(MUTI.f)
#432 rows in MUTI (species)
min(MUTI.f$MUTIscore) #-2.0492
max(MUTI.f$MUTIscore) #5.461



#############################################################
##########let's take a look at Jetz and see how many we get 

MUTI_jetz <- MUTI.f %>% mutate(Species_Jetz = SciName_MUTI)
colnames(MUTI_jetz)


JetzMUTI_UAI <- inner_join(UAI, MUTI_jetz, by="Species_Jetz") 
head(JetzMUTI_UAI)
nrow(JetzMUTI_UAI)
#330.... 

nrow(MUTI_jetz)- nrow(JetzMUTI_UAI) 
#102 species off! prob not Jetz


#############################################################
##########let's take a look at BirdLife and see how many we get 

MUTI_BL <- MUTI.f %>% mutate(Species_BirdLife = SciName_MUTI)
colnames(MUTI_BL)


BirdLifeMUTI_UAI <- inner_join(UAI, MUTI_BL, by="Species_BirdLife") 
nrow(BirdLifeMUTI_UAI)
#404.... 

nrow(MUTI_BL)- nrow(BirdLifeMUTI_UAI) 
#28 species off! This is a possibility, but eBird may do better


#############################################################
##########let's take a look at eBird and see how many we get 

MUTI_eBird <- MUTI.f %>% mutate(Species_eBird = SciName_MUTI)
colnames(MUTI_eBird)

#test for how many species match UAI and MUTI via eBird/Clements naming scheme 

eBirdMUTI_UAI <- inner_join(UAI, MUTI_eBird, by="Species_eBird") 
nrow(eBirdMUTI_UAI)
#421.... 
colnames(eBirdMUTI_UAI)

nrow(MUTI_eBird)- nrow(eBirdMUTI_UAI) 
#11 species off! this is the most cohesive naming scheme for MUTI and UAI!

#let's left_join UAI and MUTI_eBird! this will add MUTI scores to all 421 species with UAI scientific name equivalent 
colnames(UAI)
colnames(MUTI.f)
UAI_and_MUTI <- MUTI_eBird %>% 
  right_join(., UAI, by="Species_eBird") %>% # use right join to keep all rows in UAI and only those in MUTI_eBird that match UAI
  rename(CommonName_UAI = CommonName) # rename column of common names from UAI, so it is clear where this column came from
      
nrow(UAI_and_MUTI) # 4347
View(UAI_and_MUTI)
colnames(UAI_and_MUTI)

# let's find out what 11 species are missing  and then manually join them! 

nomatch_eBird_MUTI <- anti_join (MUTI_eBird, eBirdMUTI_UAI, by = "Species_eBird") %>%
  select(-Species_eBird) # drop this column as these names do not match with eBird names
nrow(nomatch_eBird_MUTI)
View(nomatch_eBird_MUTI)

#manual inspection 
#8 of the 11 not in UAI 
#3 of the 11 in UAI but not under eBird names, under Jetz names 
###### what to do about these three species? They are Brandt's Cormorant, Double-crested cormorant, and Pelagic Cormorant, so they ARE Going to be coastal. It would be nice to keep them. 

# since they are Jetz Taxonomy names, we will perform a join between nomatch_eBird_MUTI and UAI using Species_Jetz 
#clean up column names for 
colnames(UAI_and_MUTI)
colnames(nomatch_eBird_MUTI)

nomatch_eBird_MUTI.r <- nomatch_eBird_MUTI %>% 
  mutate(Species_Jetz = SciName_MUTI) %>%
  select (MUTIscore, SciName_MUTI, CommonName_MUTI, Species_Jetz, Family_MUTI) 

colnames(UAI)
colnames(UAI_and_MUTI)

UAI_and_MUTI_finalfew <- UAI %>% rename(CommonName_UAI = CommonName) %>%
  full_join(., nomatch_eBird_MUTI.r, by="Species_Jetz") %>% # use full join to match cormorants and keep all species from MUTI that are not in UAI
  filter(!is.na(MUTIscore)) %>% # keep only the rows where species have a MUTI score
  mutate(Species_Jetz = if_else(!is.na(MUTIscore) & is.na(aveUAI), NA, Species_Jetz)) %>% # remove the Species_Jetz names for species only in MUTI as we are not sure these actually align with Jetz taxonomy
  select(MUTIscore, SciName_MUTI, CommonName_MUTI, Family_MUTI, Species_eBird, Species_Jetz, Family_Jetz, Order_Jetz,
         CommonName_UAI, aveUAI, Order_eBird, Species_BirdLife) # reorder columns to match UAI_and_MUTI
  
View(UAI_and_MUTI_finalfew)
# here are the 11 species. The 3 cormorants joined to UAI using Jetz names.
# the 8 remaining species that are not found in UAI are also included but with no UAI score
# these 8 species are missing names from the 3 naming schemes
# attempt to fix this

UAI_and_MUTI_neednames <- UAI_and_MUTI_finalfew %>%
  filter(is.na(Species_eBird))

# import name conversion sheet
namesconvert <- read.csv(here("Data", "BirdNamesConversion.csv"), header=T) %>%
  rename(Species_BirdLife = Species1_BirdLife, Species_eBird = Species2_eBird, Species_Jetz = Species3_BirdTree) %>%
  select(-Avibase.ID)
  
head(namesconvert)

# try using BirdLife names
UAI_and_MUTI_namefix1 <- UAI_and_MUTI_neednames %>% 
  mutate(Species_BirdLife = SciName_MUTI) %>%
  select(-Species_Jetz, -Species_eBird) %>%
  left_join(., namesconvert, by = "Species_BirdLife")

# try using eBird names
UAI_and_MUTI_namefix2 <- UAI_and_MUTI_neednames %>% 
  mutate(Species_eBird = SciName_MUTI) %>%
  select(-Species_Jetz, -Species_BirdLife) %>%
  left_join(., namesconvert, by = "Species_eBird")

# this got all but one species (Mexican duck, Anas diazi)
# check if this is present anywhere in the names conversion data frame?
colSums(namesconvert =='Anas diazi') # no
# there is no point trying to use Jetz names to add this
# a google search of this species describes it as a lake/river duck, so not coastal
# proceeding without it

UAI_and_MUTI_namefix1b <- UAI_and_MUTI_namefix1 %>% filter(!is.na(Species_eBird)) 
UAI_and_MUTI_namefix2b <- UAI_and_MUTI_namefix2 %>% filter(!is.na(Species_BirdLife)) 
UAI_and_MUTI_namesfixed <-  bind_rows(UAI_and_MUTI_namefix1b, UAI_and_MUTI_namefix2b) %>% 
  distinct()
UAI_and_MUTI_namesfixed
# 7 species with fixed names

# we need to bind everything together
# one last issue to resolve - the 3 cormorants are in both UAI_and_MUTI_finalfew and UAI_and_MUTI
# we only want to keep the version in UAI_and_MUTI_finalfew as this has the MUTI and UAI score for these species
# use anti_join to resolve this issue and then bind the two data frames together

# check these two data frames have the same column names
# we need to do this before we use bind_rows in the next step
colnames(UAI_and_MUTI_finalfew)
colnames(UAI_and_MUTI)
colnames(UAI_and_MUTI_namesfixed)

# first use anti_join and then use bind_rows
UAI_and_MUTI_cormorants <- UAI_and_MUTI_finalfew %>%
  filter(!is.na(Species_eBird))# simply to 3 cormorant species

UAI_and_MUTI_all <- anti_join(UAI_and_MUTI, UAI_and_MUTI_cormorants, by = join_by(Species_Jetz, CommonName_UAI, aveUAI)) %>%
  bind_rows(., UAI_and_MUTI_cormorants, UAI_and_MUTI_namesfixed)
  
nrow(UAI_and_MUTI_all) # 4354. This is correct because we dropped one species (Mexican duck)
head(UAI_and_MUTI_all)
colnames(UAI_and_MUTI_all)

# perform some checks to make sure everything worked as expected

# how many rows/species have a MUTI score?
UAI_and_MUTI_all %>% filter(!is.na(MUTIscore)) %>% nrow()
# 431 which is 421 eBird species, the 3 Jetz cormorant species, and the 7 species from MUTI with no UAI match!

# how many rows/species have a UAI score?
UAI_and_MUTI_all %>% filter(!is.na(aveUAI)) %>% nrow() 
# 4347 species

# how many rows/species have both a UAI and MUTI score?
UAI_and_MUTI_all %>% filter(!is.na(aveUAI) & !is.na(MUTIscore)) %>% nrow() 
# 424 species

View(UAI_and_MUTI_all)

##################################################################
##################################################################
##################################################################
##################################################################
############# now it's time for UN Joining! ######################

UN <- read.csv (here("Data", "HuAndCardosoData.csv"), header = T)
head(UN)
UN <- UN %>% select(Species, Urban) %>% rename(Species_UN = Species)
colnames(UN)
nrow(UN) # 528 rows -- we want to keep as many of these as possible! 

###### Try joining with eBird ##################################
# it is best to try and match with eBird as all the species names are unique

eBird_UN_MUTI_UAI <- UN %>% 
  mutate(Species_eBird = Species_UN) %>%
  left_join(UAI_and_MUTI_all, ., by="Species_eBird") 

nrow(eBird_UN_MUTI_UAI) # good. All 4354 rows are present

# how many species matched?
eBird_UN_MUTI_UAI_match <- eBird_UN_MUTI_UAI %>% filter(!is.na(Species_UN)) 
nrow(eBird_UN_MUTI_UAI_match) # 337

# are they all unique? This should match the number output by the row of code above (337)
length(unique(eBird_UN_MUTI_UAI_match$Species_UN))

# how many species did not match?
nrow(UN)- nrow(eBird_UN_MUTI_UAI_match) 
# 191 species off

# get the unmatched UN species in a data frame
colnames(eBird_UN_MUTI_UAI_match)
colnames(UN)
UN_nomatch <- eBird_UN_MUTI_UAI_match %>%
  select(Species_UN, Urban) %>%
  anti_join(UN, .)

nrow(UN_nomatch) # 191
head(UN_nomatch)

###### Try matching the 191 remaining species with Jetz ##############

Jetz_UN_MUTI_UAI <- UN_nomatch %>% 
  mutate(Species_Jetz = Species_UN) %>%
  inner_join(UAI_and_MUTI_all, ., by="Species_Jetz")
View(Jetz_UN_MUTI_UAI)

# how many UAI/MUTI species matched with UN using Jetz ?
nrow(Jetz_UN_MUTI_UAI) # 116

# are they all unique?
length(unique(Jetz_UN_MUTI_UAI$Species_UN))
# 116-111 = 5
# not all are unique. This is because some UAI/MUTI species have the same scientific name using Jetz
# so certain species from UN are matching with more than one UAI/MUTI species

# look at the duplicated species
Jetz_dups <- Jetz_UN_MUTI_UAI %>% count(Species_UN) %>%
  filter(n>1) %>%
  left_join(., Jetz_UN_MUTI_UAI)

View(Jetz_dups)
# Rallus longirostris in UN is differentiated as Ridgway's Rail (Rallus obsoletus) and Clapper Rail (Rallus crepitans) in both MUTI and UAI
# this is the primary duplication that seems problematic from a coastal bird perspective

# which UN species still do not have a match with UAI?
UN_stillnomatch <- Jetz_UN_MUTI_UAI %>% distinct(Species_UN) %>%
  anti_join(UN_nomatch, .)
nrow(UN_stillnomatch) # 80 species

# tried to match using BirdLife but it did not help, so these steps were deleted
# going to assume these species are not present in UAI or MUTI data sets

# try adding BirdLife, eBird and Jetz names from names convert

head(namesconvert)
colnames(UN_stillnomatch)
nrow(UN_stillnomatch)

# try using Jetz names
UN_namefix <- UN_stillnomatch %>% 
  mutate(Species_Jetz = Species_UN) %>%
  left_join(., namesconvert, by = "Species_Jetz")

nrow(UN_namefix) # 85. some duplicates happened

UN_namedupes <- UN_namefix %>% count(Species_UN) %>%
  filter(n>1) %>%
  left_join(., UN_namefix)
View(UN_namesdupes)

# keep entries where Jetz and BirdLife name are the same
UN_namefix2 <- UN_namedupes %>%
  distinct() %>%
  mutate(clean = if_else(Species_Jetz == Species_BirdLife, "keep", "drop")) %>%
  filter(clean == "keep") %>%
  select(-clean, - n) 

# put together 75 species that were not duplicated, with the 5 resolved speices
UN_namefix_final <- UN_namefix %>% count(Species_UN) %>%
  filter(n==1) %>%
  left_join(., UN_namefix) %>%
  select(-n) %>%
  bind_rows(., UN_namefix2)

nrow(UN_namefix_final)

# how well did this work?
UN_namefix_final %>% filter(is.na(Species_BirdLife)) %>% nrow()
UN_namefix_final %>% filter(is.na(Species_eBird)) %>% nrow()
UN_namefix_final %>% filter(is.na(Species_Jetz)) %>% nrow()
# all are zero, so all names were filled in! 

##### Put matches together and build towards final data frame ######

# combine Jetz and eBird matches into one using bind_rows
# first confirm the colnames match and are in the same order
colnames(eBird_UN_MUTI_UAI_match)
colnames(Jetz_UN_MUTI_UAI)

UN_combine1 <- bind_rows(Jetz_UN_MUTI_UAI, eBird_UN_MUTI_UAI_match)
nrow(UN_combine1) # 453

# get all the species with MUTI and/or UAI with no UN match
UAI_and_MUTI_noUN <- UN_combine1 %>% select(-Species_UN, -Urban) %>%
  anti_join(UAI_and_MUTI_all, .) 
nrow(UAI_and_MUTI_noUN) # 3902

# add the Species_UN and Urban columns back to this data frame but make them NA
# we need to do this to combine these species with UN_combine2
UAI_and_MUTI_noUN[ , 'Species_UN'] = NA
UAI_and_MUTI_noUN[ , 'Urban'] = NA

# combine UN_combine1 (all species from UN with UAI or MUTI match) and all species in UAI/MUTI with no UN match 
colnames(UAI_and_MUTI_noUN)
colnames(UN_combine1)
UN_combine2 <- bind_rows(UN_combine1, UAI_and_MUTI_noUN)
nrow(UN_combine2) # 4354

# get 80 UN species with no matches that have names from all 3 naming schemes and add those
UAI_MUTI_UN_combined <- full_join(UN_combine2, UN_namefix_final)
nrow(UAI_MUTI_UN_combined) # should equal 4434 (4354 + 80)

# do all the remaining species have eBird, BirdLife and Jetz names?
# if so, these lines of code should produce zeros
UAI_MUTI_UN_combined %>% filter(is.na(Species_eBird)) %>% nrow()
UAI_MUTI_UN_combined %>% filter(is.na(Species_BirdLife)) %>% nrow()
UAI_MUTI_UN_combined %>% filter(is.na(Species_Jetz)) %>% nrow()

# There are currently multiple columns that contain family and order information
# the next script (Step 3) will involve reducing this species list to only coastal birds and one of the steps relies on the family 
# it would be helpful to have a single column for family and to use eBird families as those will best align with Birds of the World, which we will rely on to look up families
# we also need to add common names for UN species, particularly the ones that don't also have UAI or MUTI scores
# the UN data do not contain common names

# import eBird taxonomy
ebird_tax <- read.csv(here("Data","ebird_taxonomy_v2022.csv"), header=T) 

ebird <- ebird_tax %>%
  select(SCI_NAME, FAMILY) %>% rename(Species_eBird = SCI_NAME, Family_eBird = FAMILY)


colnames(UAI_MUTI_UN_combined)

UAI_MUTI_UN_families <- UAI_MUTI_UN_combined %>% 
  select(-Family_MUTI, - Order_Jetz, -Order_eBird, -Family_Jetz) %>% # remove some columns to clean things up
  left_join(., ebird)

# are any species missing a family from eBird?
UAI_MUTI_UN_families %>% filter(is.na(Family_eBird))

# does this species have a match in the eBird taxonomy using its common name (Yellow-throated scrubwren)
ebird_tax %>% filter(PRIMARY_COM_NAME =="Yellow-throated Scrubwren")
# it belongs in family Acanthizidae (Thornbills and Allies)
# manually make this fix
UAI_MUTI_UN_families$Family_eBird[UAI_MUTI_UN_families$Species_eBird == "Sericornis citreogularis"] <- "Acanthizidae (Thornbills and Allies)"

# confirm this worked by checking again that all species have a family and a common name
UAI_MUTI_UN_families %>% filter(is.na(Family_eBird))

# add common names from eBird taxonomy
UAI_MUTI_UN_common <- ebird_tax %>% select(PRIMARY_COM_NAME, SCI_NAME) %>%
  rename(CommonName_eBird = PRIMARY_COM_NAME, Species_eBird = SCI_NAME) %>%
  right_join(., UAI_MUTI_UN_families, by="Species_eBird")

# do all species have a common name?
UAI_MUTI_UN_common %>% filter(is.na(CommonName_eBird))

# manually fix this species
UAI_MUTI_UN_common$CommonName_eBird[UAI_MUTI_UN_common$Species_eBird == "Sericornis citreogularis"] <- "Yellow-throated Scrubwren"

# do all species have unique common names?
length(unique(UAI_MUTI_UN_common$CommonName_eBird))

commonname_dups <- UAI_MUTI_UN_common %>% count(CommonName_eBird) %>%
  filter(n>1) %>%
  left_join(., UAI_MUTI_UN_common)

# manually make a few fixes

# differentiate American and Northwestern Crows in the eBird common names
UAI_MUTI_UN_common$CommonName_eBird[UAI_MUTI_UN_common$CommonName_MUTI == "Northwestern Crow"] <- "Northwestern Crow"

# UN data has Little-Bronze Cuckoo and Gould's Bronze-Cuckoo so make this distinction in eBird common names
UAI_MUTI_UN_common$CommonName_eBird[UAI_MUTI_UN_common$Species_Jetz == "Chrysococcyx russatus"] <- "Gould's Bronze-Cuckoo"

# Lesser Whitethroat looks like it can be condensed into a single record with both UN and UAI
# these steps mostly accomplish this, but we still need to delete the row that no longer contains any urban scores
UAI_MUTI_UN_common$Species_UN[UAI_MUTI_UN_common$Species_Jetz == "Sylvia althaea"] <- "Sylvia curruca"
UAI_MUTI_UN_common$Urban[UAI_MUTI_UN_common$Species_Jetz == "Sylvia althaea"] <- "U"
UAI_MUTI_UN_common$Species_Jetz[UAI_MUTI_UN_common$Species_Jetz == "Sylvia althaea"] <- "Sylvia curruca"

# delete the row for Lesser Whitethroat that is no longer needed and make the final data frame
UAI_MUTI_UN_final <- UAI_MUTI_UN_common %>% slice(-3183)
nrow(UAI_MUTI_UN_final) # we now have 4433 rows or species

# do all species have unique common names now?
length(unique(UAI_MUTI_UN_final$CommonName_eBird)) # should be 4433

##### Perform some Final Checks ####

# how many species have UAI score?
UAI_MUTI_UN_final %>% filter(!is.na(aveUAI)) %>% nrow() # 4347

# how many species have MUTI score?
UAI_MUTI_UN_final %>% filter(!is.na(MUTIscore)) %>% nrow() # 431

# how many species have Urban score?
UAI_MUTI_UN_final %>% filter(!is.na(Urban)) %>% nrow() # 533
# how many are distinct species names (this should equal the number in the original data frame)?
UAI_MUTI_UN_final %>% filter(!is.na(Urban)) %>% distinct(Species_UN) %>% nrow() # 528

# how many species have UN but no MUTI or UAI?
UAI_MUTI_UN_final %>% filter(!is.na(Urban)) %>% filter(is.na(MUTIscore) & is.na(aveUAI)) %>% nrow() # 79

# how many species have MUTI but no UAI or UN?
UAI_MUTI_UN_final %>% filter(!is.na(MUTIscore)) %>% filter(is.na(Urban) & is.na(aveUAI)) %>% nrow() # 7


# last step: reorder columns to be in a more logical order
# removing common names from UAI and MUTI and keeping only eBird common names to make next steps less confusing

UAI_MUTI_UN_final <- UAI_MUTI_UN_final %>%
  rename(SciName_UN = Species_UN) %>%
  select(Species_Jetz, Species_eBird, Species_BirdLife, 
         Family_eBird, CommonName_eBird,
         aveUAI,
         SciName_MUTI, MUTIscore,
         SciName_UN, Urban)


#View for visual assessment 
View(UAI_MUTI_UN_final)

##### Export list of species and urban tolerance indices ######
saveRDS (UAI_MUTI_UN_final, here("Outputs", "UAI_MUTI_UN_final.rds"))
write.csv(UAI_MUTI_UN_final, here("Outputs", "UAI_MUTI_UN_final.csv"))
