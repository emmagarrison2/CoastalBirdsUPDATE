##########################################################
######### Coastal Birds Urban Tolerance Project ##########
##########################################################
# Step 1: UAI Names to Species
##### Author(s): SL Jennings 
##########################################################
### Objective: Align UAI Common Names with Scientific Names from eBird, Jetz and BirdLife 

# Load packages
library(here)
library(tidyverse)
library(fuzzyjoin)
library(stringi)

# Import UAI data from Neate Clegg et al. 2023
UAI <- read.csv(here("Data", "UAIs.csv"), header=T) 
head(UAI)
str(UAI)
# Species is an English common name stored as character with spaces in between words
# Both parts of the name are capitalized

# UAI data only has English common names. The taxonomy they are based on is not stated in the manuscript
# The authors are using eBird data, so potentially the common names in their data set could be most similar to the eBird taxonomy
# But the R script provided with the manuscript imports a csv file called "Traits.csv" and they use a column called "Jetz" to join trait info with UAI
# So it is not entirely clear which taxonomy the provided names best match

# Also, the UAI data has multiple rows per species. One for each city where that species occurs.
# For joining the data sets, we will simplify to keep just one row for each species 
# We will retain an average UAI for each species that represents the mean value across the cities where the species occurs

# Group by species. Calc average UAI for each species. 

UAISpp <- UAI %>% group_by(Species) %>% summarize(aveUAI=mean(UAI)) %>% 
  rename(CommonName = Species) # Rename column as CommonName to allow join with eBird taxonomy
nrow(UAISpp) # 4415 species
head(UAISpp)
summary(UAISpp)

# Import eBird taxonomy (downloaded from eBird website)
# this list contains common names and scientific names
ebird <- read.csv(here("Data","ebird_taxonomy_v2022.csv"), header=T)
# PRIMARY_COM_NAME is the English common name stored as a character with spaces between words. Each word capitalized
# SCI_NAME is the scientific name stored as a character with a space between words. Genus capitalized, species lower case

######## STEP 2: Use eBird taxonomy to add eBird scientific names to UAI data #########

# Simplify and rename eBird taxonomy to keep only columns of interest
eBirdNames <- ebird %>% rename(Species_eBird=SCI_NAME, CommonName = PRIMARY_COM_NAME, Order_eBird=ORDER1) %>% select(Species_eBird, CommonName, Order_eBird) %>% distinct() # retain only some columns
head(eBirdNames) # contain columns: Species_eBird, CommonName, Order_eBird
nrow(eBirdNames) # 16860 rows

# Join UAI and eBird using Common Names. This will add the eBird scientific names to the UAI data
UAI_eBird <- inner_join(UAISpp, eBirdNames, by="CommonName") 
nrow(UAI_eBird) # 4364, so pretty successful but not perfect

# How many species in UAI did not find a match using eBird common names?
nrow(UAISpp) - nrow(UAI_eBird)
# 51 species

# Identify the 51 species that did not match
UAIeBird_Missing <- anti_join(UAISpp, eBirdNames, by="CommonName")
View(UAIeBird_Missing) # who are they?

# Use a fuzzy join to see if we can find any just have slight spelling/formatting differences for the common names betwen UAI and eBird
UAIeBirdFuzzy <- stringdist_join(UAIeBird_Missing, eBirdNames, 
                                 by="CommonName",
                                 mode="inner",
                                 method="hamming", # counts the number of character differences that turns a into b
                                 max_dist=2, # can have up to 2 character differences
                                 distance_col="dist") 

# We need to manually inspect the resulting data frame to confirm these are indeed matches
View(UAIeBirdFuzzy) # all of these look like matches
nrow(UAIeBirdFuzzy) # this has added 11 more species

# extract 11 species with fuzzy match
UAIeBird_Found <- UAIeBirdFuzzy %>% rename(CommonName = CommonName.x) %>% select(Species_eBird, CommonName, Order_eBird, aveUAI)
# note: use CommonName.x as this is from UAI and is needed for the next step

# Get updated list of species missing match
UAIeBird_Missing2 <- anti_join(UAIeBird_Missing, UAIeBird_Found)
nrow(UAIeBird_Missing2) # 40

# export the list of 40 species that are missing a match
write.csv(UAIeBird_Missing2, here("Align_UAI", "UAIeBird_Missing2.csv"))

# used google to look up all common names from UAIeBird_Missing2 and tried to find match in eBird taxonomy
# this was often fairly easy as the eBird website often lists other common name for the species
# the city column in UAI was also used to verify some species based on their geographic range
# Some of the species were part of a group (subspecies that make up a group). These were dropped to be conservative
# Some other UAI common names did not provide an obvious/clear match to eBird names, so those were also dropped
UAIeBird_Missing_fix <- read.csv(here("Data", "UAIeBird_Missing2_fix.csv"), header=T)
head(UAIeBird_Missing_fix)

# remove rows with NA where we did not find a match
UAIeBird_fix <- UAIeBird_Missing_fix %>% filter(!is.na(Species_eBird_fix))
nrow(UAIeBird_fix) # 26 species were matched, 14 were not (removed)

### Put the "fixed" species from the fuzzy join and the manual search together with the species in UAI_eBird

# First, add the species from the fuzzy join to UAI_eBird
# Keep CommonName.y as Common Name, which  this is the entry in the eBird taxonomy
UAI_eBird_fuzzy <- UAIeBirdFuzzy %>% rename(CommonName = CommonName.y) %>% select(CommonName, aveUAI, Species_eBird, Order_eBird) %>%
  bind_rows(.,UAI_eBird)
nrow(UAI_eBird_fuzzy) # 4375

# Second, add the species from the manual search process
# Keep CommonName_eBird as Common Name
UAI_eBird_All <- UAIeBird_fix %>% dplyr::select(CommonName_eBird, aveUAI, Species_eBird_fix, Order_eBird) %>%
  rename(CommonName = CommonName_eBird, Species_eBird = Species_eBird_fix) %>%
  bind_rows(UAI_eBird_fuzzy)
nrow(UAI_eBird_All) # 4401 out of 4415 species

######## STEP 3: Use AVONET species list to add the Jetz/Bird Tree scientific names and Bird Life names to UAI data #########

# Import AVONET species list. This is from Tobias et al. (2022) Ecology Letters
# It contains scientific names for Jetz et al taxonomy, eBird taxonomy, and Bird Life
namesconvert <- read.csv(here("Data", "BirdNamesConversion.csv"), header=T)
head(namesconvert)
# Species2_eBird is the eBird species name stored as a character with spaces between words. Genus capitalized, species lower case
# Species3_BirdTree is the Jetz et al taxonomy species name stored as a character with spaces between words. Genus capitalized, species lower case

eBird_Jetz <- namesconvert %>%
  select(Species2_eBird, Species3_BirdTree) %>% # keep names for eBird and BirdTree
  rename(Species_eBird = Species2_eBird, Species_Jetz = Species3_BirdTree) %>%
  distinct() # remove any duplicate rows

nrow(eBird_Jetz) # 10773

UAI_eBird_Jetz <- UAI_eBird_All %>%
  left_join(., eBird_Jetz)
nrow(UAI_eBird_Jetz) # 4490
# should be 4401 if there is a perfect one-to-one match

length(unique(UAI_eBird_All$CommonName)) # check there are still 4401 unique common names

# which species don't have a match with Jetz?
UAI_JetzNA <- UAI_eBird_Jetz %>% filter(is.na(Species_Jetz))
nrow(UAI_JetzNA) # 53 species with no match
# probably can't easily resolve these
View(UAI_JetzNA)

# which species have multiple matches?
UAI_Jetz_multi <- UAI_eBird_Jetz %>%
  group_by(CommonName) %>%
  filter(n()>1) %>% ungroup()

nrow(UAI_Jetz_multi) # 165 rows
length(unique(UAI_Jetz_multi$CommonName)) # 76 species
View(UAI_Jetz_multi)

# get data frame with no NAs for Species_Jetz and only single matches between Jetz and eBird
UAI_eBird_Jetz_good <- UAI_eBird_Jetz %>% 
  filter(!is.na(Species_Jetz)) %>%
  anti_join(., UAI_Jetz_multi)
nrow(UAI_eBird_Jetz_good) # 4272

# For the 76 species with multiple matches, find the instances where eBird and Jetz names are identical and keep those, dropping the other name(s)
UAI_eBird_Jetz_match <- UAI_Jetz_multi %>%
  group_by(CommonName) %>%
  mutate(Match = as.numeric(if_else(Species_eBird == Species_Jetz, "1", "0"))) %>%
  filter(Match == 1) %>% select(-Match)

length(unique(UAI_eBird_Jetz_match$CommonName)) # resolved 66 species

# get the remaining species that did not have a matching name in eBird and Jetz
UAI_eBird_Jetz_NoMatch <- UAI_Jetz_multi %>%
  group_by(CommonName) %>%
  mutate(Match = as.numeric(if_else(Species_eBird == Species_Jetz, "1", "0"))) %>%
  summarize(NoMatch=sum(Match)) %>% filter(NoMatch==0) %>%
  left_join(UAI_Jetz_multi) %>% select(-NoMatch)

# write these species as csv and look up species to pick best match and manually modify
write.csv(UAI_eBird_Jetz_NoMatch, here("Data", "UAI_eBird_Jetz_NoMatch.csv"))

# read in the fixed names
UAI_eBird_Jetz_NoMatch_fix <- read.csv(here("Data", "UAI_eBird_Jetz_NoMatch_fix.csv"))

# Combine all fixed names
head(UAI_eBird_Jetz_good)
head(UAI_eBird_Jetz_NoMatch_fix)
head(UAI_eBird_Jetz_match)

UAI_eBird_Jetz_clean <- bind_rows(UAI_eBird_Jetz_good, UAI_eBird_Jetz_NoMatch_fix, UAI_eBird_Jetz_match)
nrow(UAI_eBird_Jetz_clean) # 4348

length(unique(UAI_eBird_Jetz_clean$Species_eBird)) # 4348
length(unique(UAI_eBird_Jetz_clean$Species_Jetz)) # 4209

#################### Add BirdLife Names ##########################

UAI_eBird_Jetz_clean

eBird_Jetz_BirdLife <- namesconvert %>%
  select(Species1_BirdLife, Species2_eBird, Species3_BirdTree) %>% # keep names for eBird and BirdTree
  rename(Species_BirdLife = Species1_BirdLife, Species_eBird = Species2_eBird, Species_Jetz = Species3_BirdTree) %>%
  distinct() # remove any duplicate rows

UAI_eBird_Jetz_BirdLife <- UAI_eBird_Jetz_clean %>%
  left_join(., eBird_Jetz_BirdLife) %>%
  filter(!is.na(Species_BirdLife))

nrow(UAI_eBird_Jetz_BirdLife) # 4662

# which species have multiple matches?
UAI_BirdLife_multi <- UAI_eBird_Jetz_BirdLife %>%
  group_by(CommonName) %>%
  filter(n()>1) %>% ungroup()

nrow(UAI_BirdLife_multi) # 563 rows
length(unique(UAI_BirdLife_multi$CommonName)) # 249 species
View(UAI_BirdLife_multi)

# find the instances where BirdLife, eBird and Jetz names are identical and keep those, dropping other name
UAI_BirdLife_match <- UAI_BirdLife_multi %>%
  group_by(CommonName) %>%
  mutate(Match = as.numeric(if_else(Species_BirdLife == Species_Jetz & Species_BirdLife == Species_eBird, "1", "0"))) %>%
  filter(Match == 1) %>% select(-Match)

length(unique(UAI_BirdLife_match$CommonName)) # resolved 201 species

# get the remaining species that did not have a matching names across the 3 taxonomies
UAI_BirdLife_NoMatch <- UAI_BirdLife_multi %>%
  group_by(CommonName) %>%
  mutate(Match = as.numeric(if_else(Species_BirdLife == Species_Jetz & Species_BirdLife == Species_eBird, "1", "0"))) %>%
  summarize(NoMatch=sum(Match)) %>% filter(NoMatch==0) %>%
  left_join(UAI_BirdLife_multi) %>% select(-NoMatch)

# try getting match with only eBird and BirdLife
UAI_BirdLife_eBird <- UAI_BirdLife_NoMatch %>%
  group_by(CommonName) %>%
  mutate(Match = as.numeric(if_else(Species_BirdLife == Species_eBird, "1", "0"))) %>%
  filter(Match==1) %>% select(-Match)

nrow(UAI_BirdLife_eBird) # 36 species

UAI_BirdLife_stillmultimatch <- bind_rows(UAI_BirdLife_match, UAI_BirdLife_eBird) %>% select(CommonName) %>%
  anti_join(UAI_BirdLife_NoMatch, .)

# Last attempts to resolve multiple matches
# Keep record if BirdLife matches Jetz
# Finally, if the species part of the scientific name is a match for Jetz and eBird but the genus name is different, keep that record
# Drop any other species

# records where species name for Jetz and BirdLife match
UAI_BirdLife_Jetz <- UAI_BirdLife_stillmultimatch %>%
  group_by(CommonName) %>%
  mutate(Match = as.numeric(if_else(Species_BirdLife == Species_Jetz, "1", "0"))) %>%
  filter(Match==1) %>% select(-Match)

# records where species part of genus-species matches for eBird and Bird Life
UAI_BirdLife_finalfew <- UAI_BirdLife_Jetz %>% select(CommonName) %>%
  anti_join(UAI_BirdLife_stillmultimatch, .) %>% # get remaining species not resolved by previous step
 separate_wider_delim(Species_eBird, delim=" ", names = c("genus_eBird", "spp_eBird"), cols_remove=F) %>%
  separate_wider_delim(Species_BirdLife, delim=" ", names = c("genus_BirdLife", "spp_BirdLife"), cols_remove=F) %>%
  mutate(Match = if_else(spp_eBird == spp_BirdLife, 1, 0)) %>% filter(Match==1) %>% 
  select(CommonName, aveUAI, Species_eBird, Order_eBird, Species_Jetz, Species_BirdLife)

########
# put everything together
head(UAI_BirdLife_match) # 201 species
head(UAI_BirdLife_eBird) # 36 species
head(UAI_BirdLife_Jetz) # 7 species
head(UAI_BirdLife_finalfew) # 4 species

# get data frame with no NAs for Species_Jetz and only single matches between Jetz and eBird
UAI_BirdLife_good <- UAI_eBird_Jetz_BirdLife %>% 
  filter(!is.na(Species_BirdLife)) %>%
  anti_join(., UAI_BirdLife_multi)

nrow(UAI_BirdLife_good) # 4099
head(UAI_BirdLife_good)

# Combine everything
UAI_eBirdJetzBirdLife_combine <- bind_rows(
  UAI_BirdLife_good,
  UAI_BirdLife_match,
  UAI_BirdLife_eBird,
  UAI_BirdLife_Jetz,
  UAI_BirdLife_finalfew
)

nrow(UAI_eBirdJetzBirdLife_combine) # 4347

length(unique(UAI_eBirdJetzBirdLife_combine$Species_BirdLife)) # 4312
length(unique(UAI_eBirdJetzBirdLife_combine$Species_eBird)) # 4347
length(unique(UAI_eBirdJetzBirdLife_combine$Species_Jetz)) # 4208
# we still have instances of where there is not a one-to-one match
# however, these may disappear as we refine this list to only include coastal birds
# we will check again once the list contains only coastal birds and decide how to resolve any species if needed

# finally, add the family and order
jetznames <- read.csv(here("Data", "JetzNames.csv"), header=T)
head(jetznames)

UAI_eBirdJetzBirdLife_final <- jetznames %>% rename(Species_Jetz = SciName, Family_Jetz = Family, Order_Jetz = Order) %>%
  right_join(., UAI_eBirdJetzBirdLife_combine)

# export as csv and rds
write.csv(UAI_eBirdJetzBirdLife_final, here("Outputs", "UAI_eBirdJetzBirdLife_final.csv"))
saveRDS(UAI_eBirdJetzBirdLife_final, here("Outputs", "UAI_eBirdJetzBirdLife_final.rds"))

        