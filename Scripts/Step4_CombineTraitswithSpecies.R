#The objective of this script is to join Coastal_Birds_list.rds to trait data. 
#We will create a joined dataframe for each category of predictor trait variables, as follows... 
###### SENSORY TRAITS 
###### DIET TRAITS 
###### NESTING TRAITS
###### LIFE HISTORY TRAITS 
###### SEXUAL SELECTION TRAITS 
###### SOCIAL TRAITS 


# Load packages
library(here)
library(tidyverse)
library(stringi)


# read in Coastal Bird List 

Coastal_Birds_list <- readRDS (here("Data", "Coastal_Birds_list.rds"))
View(Coastal_Birds_list)
nrow(Coastal_Birds_list)
#807 total coastal species 

#let's see how many birds have UAI scores, MUTI scores, and UN scores 

Coastal_birds_UAI <- Coastal_Birds_list %>% 
  filter(!is.na(aveUAI))
nrow(Coastal_birds_UAI)
#798 species 

Coastal_birds_MUTI <- Coastal_Birds_list %>% 
  filter(!is.na(MUTIscore))
nrow(Coastal_birds_MUTI)
#130 species 

Coastal_birds_UN <- Coastal_Birds_list %>% 
  filter(!is.na(Urban))
nrow(Coastal_birds_UN)
#129 species 

################# JOIN BODY MASS #################
# body mass from avonet (Tobias et al. 2022 Ecol Letters)
# this join between coastal species and body mass will be used as the starting df to add other predictor traits 

AVONET <- read.csv(here("Data", "avonet.csv"))
head(AVONET)

head(Coastal_Birds_list)

# we are looking to join body mass (Mass (g)) with Coastal_Birds_list based on scientific name for species. 
# Avonet uses Bird Life names, so join via Coastal_Birds_list column Species_BirdLife

AVONET.2 <- AVONET %>% 
  rename(Species_BirdLife = Species1) 

AVONET.3 <- AVONET.2 %>% 
  select(Species_BirdLife, Mass)

colnames(AVONET.3)
head(AVONET.3)

Coastal_Bodymass <- left_join(Coastal_Birds_list, AVONET.3, by="Species_BirdLife")
head(Coastal_Bodymass)
nrow(Coastal_Bodymass) #807, looks good. 

# Body mass in birds is typically very skewed. Check that for coastal birds
hist(Coastal_Bodymass$Mass)

# look at distribution of log transformed body mass
hist(log(Coastal_Bodymass$Mass)) # much better

# add column for log transformed body mass
Coastal_Bodymass <- Coastal_Bodymass %>%
  mutate(Mass_log = log(Mass))

#check to see if we have NAs for body mass (Mass)
Mass_test <- Coastal_Bodymass %>% 
  filter(!is.na(Mass))
nrow(Mass_test)
#all 807 species have a body mass value from Avonet 

#save Coastal_Bodymass as .rds for easy retrieval 
saveRDS(Coastal_Bodymass, here("Outputs", "Coastal_Species_w_Mass.rds"))

##################################################################
################## JOIN SENSORY TRAITS ###########################
##################################################################

# two sensory traits to join: peak vocal frequency and dim light vision (C.T. ratio)

# # # # # # Peak Vocal Frequency # # # # # #

# data sources: Hu and Cardoso 2009 and Mikula at al. 2021

# import Hu and Cardoso frequency measurements
HU_CARDOSO <- read.csv(here("Data", "HuandCardosoData.csv"), header=T)

head(HU_CARDOSO)
# column Domi contains the dominant/peak frequency in Khz
# should join to SciName_UN column

# import Mikula at el song frequency measurements for song birds
MIKULA <- read.csv(here("Data", "mikula_peakfreq.csv"), header=T)
head(MIKULA)
# column peak_frequency contains the dominant/peak frequency in Hz
# uses Bird Tree names

# steps: first join Hu and Cardoso
freqHC <- HU_CARDOSO %>% 
  rename(SciName_UN = Species,) %>% # rename column to allow for join
  mutate(peak_freq_HC = Domi*1000) %>% # convert peak freq measurements to Hz to match Mikula et al.
  select(-Mini, -Body, -Veget, -Domi, -Urban) %>% # remove some columns
  left_join(Coastal_Bodymass, ., by="SciName_UN" )

# now join Mikula measurements 
# combine data from both sources into a single column
Coastal_freq <- MIKULA %>%
  mutate(Species_Jetz = str_replace(genus_species, "_", " ")) %>% # change underscore to space and rename column to allow join
  rename(peak_freq_Mik = peak_frequency) %>%
  select(Species_Jetz, peak_freq_Mik) %>% # select only the columns we need to keep
  left_join(freqHC, ., by="Species_Jetz") %>%
  mutate(peak_freq = case_when( # make new colum "peak_freq" that combines the two sources of freq measurements
    peak_freq_HC > 0 & peak_freq_Mik > 0 ~ peak_freq_Mik, # if both sources have values, keep freq value from Mikula et al
    peak_freq_HC > 0 & is.na(peak_freq_Mik) ~ peak_freq_HC, # if only Hu and Cardoso has freq value, keep it
    is.na(peak_freq_HC) & peak_freq_Mik >  0 ~ peak_freq_Mik # if only Hu and Cardoso has freq value, keep it
  )) %>%
  select(-peak_freq_Mik, -peak_freq_HC) # remove these columns as a final step to keep peak_freq as only column with vocal frequency

# how many species have peak freq measurements?
# for UAI?
Coastal_freq %>% filter(!is.na(peak_freq)) %>% filter(!is.na(aveUAI)) %>% nrow()
# for MUTI?
Coastal_freq %>% filter(!is.na(peak_freq)) %>% filter(!is.na(MUTIscore)) %>% nrow()
# for UN?
Coastal_freq %>% filter(!is.na(peak_freq)) %>% filter(!is.na(Urban)) %>% nrow()

# # # # # # Dim Light Vision # # # # # #

# data sources: multiple sources compiled by Francis lab at Cal Poly

# import eye morphometrics
EYE <- read.csv(here("Data", "eye.csv"), header=T)
head(EYE)
# we need column "C.T"

# remove rows/species with no C.T values and simplify data frame in preparation to join
eye_CT <- EYE %>% filter(C.T > 0) %>% 
  select(Scientific, C.T) %>% # keep only necessary columns
  rename(Species_Jetz = Scientific) # rename column with scientific names to enable join with coastal birds

# join eye data to Coastal_freq
Coastal_Sensory <- left_join(Coastal_freq, eye_CT, by = "Species_Jetz")

# confirm we still have all the coastal species
nrow(Coastal_Bodymass) == nrow(Coastal_Sensory) # should be "TRUE"

# how many species have dim light vision values?
# for UAI?
Coastal_Sensory %>% filter(!is.na(C.T)) %>% filter(!is.na(aveUAI)) %>% nrow()
# for MUTI?
Coastal_Sensory %>% filter(!is.na(C.T)) %>% filter(!is.na(MUTIscore)) %>% nrow()
# for UN?
Coastal_Sensory %>% filter(!is.na(C.T)) %>% filter(!is.na(Urban)) %>% nrow()

# export Coastal_Sensory as rds to be used in models
saveRDS(Coastal_Sensory, here("Outputs", "Coastal_Species_Sensory.rds"))

##################################################################
######################### JOIN DIET TRAITS #######################
##################################################################
# diet traits: % invertebrates, % vertebrates, % plant/seed, % fruit/nectar
# from Wilman et al. 2014 (elton traits)

# import elton traits

ELTON_DIET <- read.csv(here("Data", "elton.csv"))

head(ELTON_DIET)


# Since Elton Traits use BirdTree/Jetz taxonomy, let's join by Species_Jetz 

#change name of Species column in ELTON_DIET to Species_Jetz 

ELTON_DIET2 <- ELTON_DIET %>% 
  rename(Species_Jetz = Scientific)

ELTON_DIET3 <- ELTON_DIET2 %>% 
  select(Species_Jetz, Diet.Inv, Diet.Vend, Diet.Vect, Diet.Vfish, Diet.Vunk, Diet.Scav, Diet.Fruit, Diet.Nect, Diet.Seed, Diet.PlantO)

colnames(ELTON_DIET3)
head(ELTON_DIET3)

# join diet traits with species
Coastal_Diet_Traits <- left_join(Coastal_Bodymass, ELTON_DIET3)

# view to check in on it 
View(Coastal_Diet_Traits)
nrow(Coastal_Diet_Traits) #807, as it should be 


#check to see if we have NAs for Diet traits (test using Diet.Inv)

Diet_test <- Coastal_Diet_Traits %>% 
  filter(!is.na(Diet.Inv))
nrow(Diet_test)
#all 807 species have diet traits from Wilman et al. 2014 (elton traits)

# Next, we will simplify the Diet Traits into fewer categories

######################################################################################
# Start binning our categories for diet! 4 new bins --> vert, invert, plant/seed, and fruit/nectar (which we will test, but may have to drop)
# First let's combine diet of vertebrate endotherms, vertbrate ectotherms, vertebrate unknown, vertebrate fish, and scavenging together to make a general "Diet Vertebrate" bin 

Coastal_Diet_Traits$Diet.Vend <- as.numeric(Coastal_Diet_Traits$Diet.Vend)
Coastal_Diet_Traits$Diet.Vect <- as.numeric(Coastal_Diet_Traits$Diet.Vect)
Coastal_Diet_Traits$Diet.Vunk <- as.numeric(Coastal_Diet_Traits$Diet.Vunk)
Coastal_Diet_Traits$Diet.Vfish <- as.numeric(Coastal_Diet_Traits$Diet.Vfish)
Coastal_Diet_Traits$Diet.Scav <- as.numeric(Coastal_Diet_Traits$Diet.Scav)

Coastal_Diet_Traits <- Coastal_Diet_Traits %>%
  rowwise() %>%
  mutate(Diet.Vert = Diet.Vend + Diet.Vect + Diet.Vunk + Diet.Vfish + Diet.Scav)
head(Coastal_Diet_Traits)

# remove the columns that we used to make Diet.Vert

Coastal_Diet_Traits2 <- Coastal_Diet_Traits %>% select (-Diet.Vend, -Diet.Vect, -Diet.Vunk, -Diet.Vfish, - Diet.Scav)
colnames(Coastal_Diet_Traits2)

#
# combine fruit and nectar into a single column "Diet Fruit / Nectar" or Diet.FN
# 

Coastal_Diet_Traits2$Diet.Fruit <- as.numeric(Coastal_Diet_Traits2$Diet.Fruit)
Coastal_Diet_Traits2$Diet.Nect <- as.numeric(Coastal_Diet_Traits2$Diet.Nect)

Coastal_Diet_Traits2 <- Coastal_Diet_Traits2 %>%
  rowwise() %>%
  mutate(Diet.FN = Diet.Fruit + Diet.Nect)

head(Coastal_Diet_Traits2)


# remove columns that we used to make Diet.FN

Coastal_Diet_Traits3 <- Coastal_Diet_Traits2 %>% select (-Diet.Fruit, -Diet.Nect)
colnames(Coastal_Diet_Traits3)

#
# combine Plant Other and Seed into a new column Diet.PS (Diet Plant/Seed)
#

Coastal_Diet_Traits3$Diet.PlantO <- as.numeric(Coastal_Diet_Traits3$Diet.PlantO)
Coastal_Diet_Traits3$Diet.Seed <- as.numeric(Coastal_Diet_Traits3$Diet.Seed)

Coastal_Diet_Traits3 <- Coastal_Diet_Traits3 %>%
  rowwise() %>%
  mutate(Diet.PS =  Diet.PlantO + Diet.Seed)

#remove columns that we used to make Diet.FN

Coastal_Diet_Traits4 <- Coastal_Diet_Traits3 %>% select (-Diet.PlantO, -Diet.Seed)
colnames(Coastal_Diet_Traits4)

# how many species have diet values?
# using Diet.Inv for this test but it should be the same across the 4 diet categories
# for UAI?
Coastal_Diet_Traits4 %>% filter(!is.na(Diet.Inv)) %>% filter(!is.na(aveUAI)) %>% nrow()
# for MUTI?
Coastal_Diet_Traits4 %>% filter(!is.na(Diet.Inv)) %>% filter(!is.na(MUTIscore)) %>% nrow()
# for UN?
Coastal_Diet_Traits4 %>% filter(!is.na(Diet.Inv)) %>% filter(!is.na(Urban)) %>% nrow()

# save Coastal_Diet_Traits as .rds for easy retrieval 
saveRDS(Coastal_Diet_Traits4, here("Outputs", "Coastal_Species_Diet.rds"))


##################################################################
################# LIFE HISTORY TRAITS ############################
##################################################################

# life history traits include: longevity, clutch size, brood value, developmental mode
# sources: longevity from Bird et al. 2020, clutch size from Myhrvold et al. 2015,
# sources continued: brood value calculated using Bird et al. and Myhrvold et al.
# sources continued:developmental mode from Delhey et al. 2023 (originally from Ning et al. 2016)

# # # # # # Longevity # # # # # #

# import longevity data from Bird et al.
BIRD <- read.csv(here("Data", "longevity.csv"), header=T)
head(BIRD)

# prepare trait data frame for join
longevity <- BIRD %>%
  rename(Species_BirdLife = Scientific.name, longevity = Maximum.longevity) %>%
  select(Species_BirdLife, longevity)

# join longevity to Coastal_Bodymass
Coastal_Longevity <- left_join(Coastal_Bodymass, longevity, by = "Species_BirdLife")

# confirm all species are still present
nrow(Coastal_Bodymass) == nrow(Coastal_Longevity)

# how many species have longevity values?
# for UAI?
Coastal_Longevity %>% filter(!is.na(longevity)) %>% filter(!is.na(aveUAI)) %>% nrow()
# for MUTI?
Coastal_Longevity %>% filter(!is.na(longevity)) %>% filter(!is.na(MUTIscore)) %>% nrow()
# for UN?
Coastal_Longevity %>% filter(!is.na(longevity)) %>% filter(!is.na(Urban)) %>% nrow()

# # # # # # Clutch Size # # # # # #

# import clutch data from Myhrvold et al. 2015
AMNIOTE <- read.csv(here("Data", "amniote.csv"), header= T) # Myhrvold et al. 2015 (egg/clutch traits)
head(AMNIOTE)

# modify the clutch data frame slightly to make join easier
clutch <- AMNIOTE %>% 
  filter(class=="Aves") %>% # filter data to retain only the info for birds
  mutate(scientific_name = paste(genus, species, sep=" ")) %>% # there is no column with full sci name (genus and species), so creating one
  mutate(across(where(is.numeric), ~na_if(., -999))) %>% # this spreadsheet uses -999 if there is missing data. Changing missing values to NA 
  rename(clutch_size = litter_or_clutch_size_n, clutches_per_y = litters_or_clutches_per_y) %>% # simplify column names slightly
  select(scientific_name, clutch_size, clutches_per_y) # select columns we want to keep

head(clutch)

# It is not obvious which naming scheme this trait databases uses
# first, try joining to using BirdLife scientific names
# joining to Coastal_Longevity
Coastal_clutch_BL <- clutch %>% 
  rename(Species_BirdLife = scientific_name) %>% # rename columns to allow join
  left_join(Coastal_Longevity, ., by="Species_BirdLife") 

# are any species missing clutch size?
Coastal_clutchsize_NA <- Coastal_clutch_BL %>%
  filter(is.na(clutch_size))

nrow(Coastal_clutchsize_NA) # how many are missing?

# are any species missing number of clutches per yr?
Coastal_clutchyr_NA <- Coastal_clutch_BL %>%
  filter(is.na(clutches_per_y))

nrow(Coastal_clutchyr_NA) # how many are missing?
# as expected more species are missing litters or clutches per year. this is the case in the original data that contains all the species

# get a list of species that are missing both clutch size and clutches per y
# if they have one or the other, it means the species name matched using Bird Life, and the NA is due to no value being present in the Myhrvold data
# but if they are missing both, then it could be due to mismatch in names
Coastal_clutch_missing <- Coastal_clutch_BL %>%
  mutate(status = if_else(is.na(clutch_size) & is.na(clutches_per_y), "missing", "okay")) %>% # label species that are missing both, all others are "okay"
  filter(status == "missing") %>% # keep only rows where the status is "missing"
  dplyr::select(-status, -clutch_size, -clutches_per_y) # drop status column as it is not needed for future steps
# also drop clutch size and clutches per year to not mess up join in next step
nrow(Coastal_clutch_missing) # how many are missing?

# we also need a data frame of all the species from UAI_clutch_BL that did succesfully match
Coastal_clutch_okay <- Coastal_clutch_BL %>%
  mutate(status = if_else(is.na(clutch_size) & is.na(clutches_per_y), "missing", "okay")) %>% # label species that are missing both, all others are "okay"
  filter(status == "okay") %>% # keep only rows where the status is "okay"
  dplyr::select(-status) # drop status column as it is not needed for future steps

# try joining species that are missing clutch traits using Jetz names
Coastal_clutch_Jetz <- clutch %>% 
  rename(Species_Jetz = scientific_name) %>% 
  dplyr::select(Species_Jetz, clutch_size, clutches_per_y) %>%
  right_join(., Coastal_clutch_missing, by="Species_Jetz") 

# how many species were added for clutch size?
Coastal_clutch_Jetz %>% filter(! is.na(clutch_size)) %>% nrow()

# how many species were added for number of clutches per y?
Coastal_clutch_Jetz %>% filter(! is.na(clutches_per_y)) %>% nrow()

# combine Coastal_clutch_Jetz and Coastal_clutch_okay (from Bird Life join) to get the "max" list
# check data frames contain the same columns before binding them together
colnames(Coastal_clutch_Jetz)
colnames(Coastal_clutch_okay)

Coastal_Clutch <- bind_rows(Coastal_clutch_okay, Coastal_clutch_Jetz)

# check this worked. Should equal "TRUE" if all species are present
nrow(Coastal_Bodymass) == nrow(Coastal_Clutch)


# how many species have clutch size values?
# for UAI?
Coastal_Clutch %>% filter(!is.na(clutch_size)) %>% filter(!is.na(aveUAI)) %>% nrow()
# for MUTI?
Coastal_Clutch %>% filter(!is.na(clutch_size)) %>% filter(!is.na(MUTIscore)) %>% nrow()
# for UN?
Coastal_Clutch %>% filter(!is.na(clutch_size)) %>% filter(!is.na(Urban)) %>% nrow()

# # # # # # Brood Value # # # # # #

#  calculate brood value using longevity from Bird et al. and clutches_per_y from Mryhvold et al

# check the variables (clutches per year and longevity) that will be used for the calculation
clutch.yr <- Coastal_Clutch %>% filter(!is.na(clutches_per_y)) %>% select(clutches_per_y)
range(clutch.yr$clutches_per_y) # 21 seems high. Will inspect
hist(clutch.yr$clutches_per_y, breaks=50)
# clutches per yr for Australian Brushturkey and Killdeer seem inaccurate. Check values on Birds of the World

# Australian Brushturkey: a female lays 15-27 eggs in a season
# the eggs are often can be spread out over separate nests that are tended by different males (they are polyandrous) 
# If each nest where she lays is considered a separate brood then this high value is accurate
# leaving as is for now

# Killdeer: clutch size is typically 4 eggs. Birds of the world mentions one brood reared per season
# value in data frame needs to be changed

Coastal_Clutch$clutches_per_y[Coastal_Clutch$CommonName_eBird == "Killdeer"] <- "1.00"

max.longevity <- Coastal_Clutch %>% filter(!is.na(longevity)) %>% select(longevity)
range(max.longevity$longevity) # seems reasonable as there are long-lived seabirds in the data
hist(max.longevity$longevity)

# perform calculation to get brood value
# building off of Coastal_Clutch to do this
Coastal_BroodValue <- Coastal_Clutch %>%
  mutate_at(vars(longevity, clutches_per_y), as.numeric) %>%
  rowwise() %>%
  mutate(brood_value = log(1/(longevity*clutches_per_y)))

# how many species have clutch size values?
# for UAI?
Coastal_BroodValue %>% filter(!is.na(brood_value)) %>% filter(!is.na(aveUAI)) %>% nrow()
# for MUTI?
Coastal_BroodValue %>% filter(!is.na(brood_value)) %>% filter(!is.na(MUTIscore)) %>% nrow()
# for UN?
Coastal_BroodValue %>% filter(!is.na(brood_value)) %>% filter(!is.na(Urban)) %>% nrow()

# # # # # # Developmental Mode # # # # # #

# import trait data from Delhey et al. 2023
DELHEY <- read.csv(here("Data", "Delhey_2023_DS7.csv"), header=T)
head(DELHEY)

# simplify to keep only developmental mode and scientific names
DELHEY_developmental <- DELHEY %>%
  mutate(Species_Jetz = str_replace(phylo, "_", " ")) %>% # rename sci names and remove underscore
  select(Species_Jetz, developmental_mode) %>%
  distinct() # there are duplicate rows because many species have males and females. Remove them

# join developmental mode to all the other life history traits
Coastal_LifeHistory <- left_join(Coastal_BroodValue, DELHEY_developmental, by = "Species_Jetz")

# check we still have all the species
nrow(Coastal_LifeHistory) == nrow(Coastal_Bodymass)

# how many species have developmental mode?
# for UAI? 
# Print numbers for Precocial = 0 and Altricial = 1
Coastal_LifeHistory %>% filter(!is.na(developmental_mode)) %>% filter(!is.na(aveUAI)) %>% group_by(developmental_mode) %>% count()
# for MUTI? 
# Print numbers for Precocial = 0 and Altricial = 1
Coastal_LifeHistory %>% filter(!is.na(developmental_mode)) %>% filter(!is.na(MUTIscore)) %>% group_by(developmental_mode) %>% count()
# for UN? 
# Print numbers for Precocial = 0 and Altricial = 1
Coastal_LifeHistory %>% filter(!is.na(developmental_mode)) %>% filter(!is.na(Urban)) %>% group_by(developmental_mode) %>% count()

# save joined Life History traits and Coastal Species as an .rds file 
saveRDS(Coastal_LifeHistory, here("Outputs", "Coastal_Species_LifeHistory.rds"))


##################################################################
######################## NEST TRAITS #############################
##################################################################

# nest traits include: nest site (low, high), nest strategy (open or enclosed), nest safety
# data sources include Chia et al. 2023 (nest site and nest structure), Delhey et al. 2023 (nest safety)


# import nest site and strategy info from Chia et al 2023
# uses the BirdLife taxonomy

CHIA_NEST <- read.csv(here("Data", "nests.csv"))

head(CHIA_NEST)
colnames(CHIA_NEST)

# rename column to facilitate join
# Chia et al. uses Bird Life names
nests_names <- CHIA_NEST %>%
  rename(Species_BirdLife = Scientific_name) %>%
  select(Species_BirdLife, Parasite:NestStr_second_cavity) # keep only necessary columns

# join nest traits with Coastal_Bodymass
Coastal_Nest_Chia <- left_join(Coastal_Bodymass, nests_names)
# note: there are a few species where nest data is missing or is partially missing

# # # # # # NEST STRATEGY # # # # # #
# sort nest strategies into "open" and "enclosed" (NestStr = Nest Strategy)

# classify nests as Open and Enclosed
# we also need to identify species that fall into both Open and Enclosed nest types
# Open nests types are scrape, platform and cup
# Enclosed nest types are dome, tunnel, primary_cavity, second_cavity
Coastal_Nest_OpEnc <- Coastal_Nest_Chia %>%
  mutate(NestStr_Open = ifelse(NestStr_scrape > 0 | NestStr_platform > 0 | NestStr_cup > 0, 1, 0),
         NestStr_Enclosed = ifelse(NestStr_dome > 0 | NestStr_dome_tunnel > 0 | NestStr_primary_cavity > 0 | NestStr_second_cavity > 0, 1, 0),
         OpenAndEnclosed = ifelse( NestStr_Open == 1 & NestStr_Enclosed == 1, "both", "one"))

# which species and how many were assigned to both Open and Enclosed groups?
Coastal_Nest_StrBoth <- Coastal_Nest_OpEnc %>% filter(OpenAndEnclosed == "both")
nrow(Coastal_Nest_StrBoth)

# we will resolve these into a single category based on the nest type they predominantly use
# any that can't be resolved will be marked as NA
# this may occur if they use an mix of open and enclosed nests and there is no obvious preference described in Birds of World species page

# export for editing
write.csv(Coastal_Nest_StrBoth, here("Notes", "nest_structure_both.csv"))

# import the edited version
Coastal_Nest_Str_edited <- read.csv(here("Notes", "nest_structure_both_edited.csv"), header=T)
head(Coastal_Nest_Str_edited)
nrow(Coastal_Nest_Str_edited)

# aiming for a single column NestStr where 0 is enclosed and 1 is open

# modify the data frame with the recently edited species to have just a single column for Nest Strategy (NestStr)
Coastal_Nest_Str_1 <- Coastal_Nest_Str_edited %>%
  select(Species_Jetz:NestSite_termite_ant, NestStr)

# modify the full species list data frame (Coastal_Nest_OpEnc) to have just a single column for Nest Strategy (NestStr)
# remove the species that needed resolving in previous steps
Coastal_Nest_Str_2 <- Coastal_Nest_OpEnc %>%
  mutate(NestStr = case_when(
    NestStr_Open == 1 & NestStr_Enclosed == 0 ~ 1, # open nests get a value of 1
    NestStr_Open == 0 & NestStr_Enclosed == 1 ~ 0, # enclosed nests get a value of 0
    NestStr_Open == 1 & NestStr_Enclosed == 1 ~ NA # nests that needed to be manually resolved/edited get NA and will be fixed in next step
  )) %>%
  select(Species_Jetz:NestSite_termite_ant, NestStr) %>% # identify columns to keep
  filter(!is.na(NestStr)) # remove all the species that needed editing/resolving or that are missing from Chia et al. nest database

nrow(Coastal_Nest_Str_2)

# get the species that have resolved nest strategies and 10 species that are missing from nest database
# add the modified nest strategy for the species that were looked up
Coastal_Nest_Str_3 <- Coastal_Nest_OpEnc %>%
  mutate(NestStr = case_when(
    NestStr_Open == 1 & NestStr_Enclosed == 0 ~ 1, # open nests get a value of 1
    NestStr_Open == 0 & NestStr_Enclosed == 1 ~ 0, # enclosed nests get a value of 0
    NestStr_Open == 1 & NestStr_Enclosed == 1 ~ NA # nests that needed to be manually resolved/edited get NA and will be fixed in next step
  )) %>%
  select(Species_Jetz:NestSite_termite_ant, NestStr) %>% # identify columns to keep
  filter(is.na(NestStr)) %>%
  left_join(., Coastal_Nest_Str_1)

# combine Coastal_Nest_Str_3 and Coastal_Nest_Str_2 to get final data frame
Coastal_Nest_Str <- bind_rows(Coastal_Nest_Str_3, Coastal_Nest_Str_2)
nrow(Coastal_Nest_Str) == nrow(Coastal_Bodymass)

# how many species have nest strategy data ?
# for UAI? 
# Print numbers for Open and Enclosed (Open = 1, Enclosed = 0)
Coastal_Nest_Str %>% filter(!is.na(NestStr)) %>% filter(!is.na(aveUAI)) %>% group_by(NestStr) %>% count()
# for MUTI? 
# Print numbers for Open and Enclosed (Open = 1, Enclosed = 0)
Coastal_Nest_Str %>% filter(!is.na(NestStr)) %>% filter(!is.na(MUTIscore)) %>% group_by(NestStr) %>% count()
# for UN? 
# Print numbers for Open and Enclosed (Open = 1, Enclosed = 0)
Coastal_Nest_Str %>% filter(!is.na(NestStr)) %>% filter(!is.na(Urban)) %>% group_by(NestStr) %>% count()


# # # # # # NEST SITE # # # # # # 
# sort nest sites into the bins "low" (on or underground) and "high" (above the ground)
# low nest sites include ground, underground, waterbody
# high nest sites include tree, nontree, cliff_bank, termite_ant
Coastal_Nest_Site <-  Coastal_Nest_Str %>% 
  mutate(NestSite_Low = ifelse(NestSite_ground > 0 | NestSite_underground > 0 | NestSite_waterbody > 0, 1, 0),
         NestSite_High = ifelse(NestSite_tree > 0 | NestSite_nontree > 0 | NestSite_cliff_bank > 0 | NestSite_termite_ant > 0, 1, 0),
         LowAndHigh = ifelse( NestSite_Low == 1 & NestSite_High == 1, "both", "one"))

# which species and how many were assigned to both Low and High groups?
Coastal_Nest_SiteBoth <- Coastal_Nest_Site %>% filter(LowAndHigh == "both")
nrow(Coastal_Nest_SiteBoth) # many! 
# Many coastal species nest on the ground and on cliffs/banks.
# This phenomenon is likely driving such a high number of species falling into both nest site types
# This can't easily be resolved
# Instead, we will use two models one to examine Low species and one to examine High species

# remove all unnecessary columns to keep data frame clean and organized
Coastal_Nest_StrSite <- Coastal_Nest_Site %>% 
  select(Species_Jetz:Mass, NestStr, NestSite_High, NestSite_Low)

colnames(Coastal_Nest_StrSite)

# how many species have nest site high data ?
# for UAI? 
# Print numbers for both groups. A value of 1 is High nests.
Coastal_Nest_StrSite %>% filter(!is.na(NestSite_High)) %>% filter(!is.na(aveUAI)) %>% group_by(NestSite_High) %>% count()
# for MUTI? 
# Print numbers for both groups. A value of 1 is High nests.
Coastal_Nest_StrSite %>% filter(!is.na(NestSite_High)) %>% filter(!is.na(MUTIscore)) %>% group_by(NestSite_High) %>% count()
# for UN? 
# Print numbers for both groups. A value of 1 is High nests.
Coastal_Nest_StrSite %>% filter(!is.na(NestSite_High)) %>% filter(!is.na(Urban)) %>% group_by(NestSite_High) %>% count()


# how many species have nest site low data ?
# for UAI? 
# Print numbers for both groups. A value of 1 is for Low nests.
Coastal_Nest_StrSite %>% filter(!is.na(NestSite_Low)) %>% filter(!is.na(aveUAI)) %>% group_by(NestSite_Low) %>% count()
# for MUTI? 
# Print numbers for both groups. A value of 1 is for Low nests.
Coastal_Nest_StrSite %>% filter(!is.na(NestSite_Low)) %>% filter(!is.na(MUTIscore)) %>% group_by(NestSite_Low) %>% count()
# for UN? 
# Print numbers for both groups. A value of 1 is for Low nests.
Coastal_Nest_StrSite %>% filter(!is.na(NestSite_Low)) %>% filter(!is.na(Urban)) %>% group_by(NestSite_Low) %>% count()


# # # # # # Nest Safety # # # # # #

# this comes from Delhey et al. 2023 which was already imported above
head(DELHEY)

# simplify to keep only nest safety and scientific names
DELHEY_safety <- DELHEY %>%
  mutate(Species_Jetz = str_replace(phylo, "_", " ")) %>% # rename sci names and remove underscore
  select(Species_Jetz, nest.safety) %>%
  distinct() # there are duplicate rows because many species have males and females. Remove them

# join nest safety to Coastal_Nest_StrSite
Coastal_Nest_Traits <- left_join(Coastal_Nest_StrSite, DELHEY_safety, by = "Species_Jetz")

# confirm all the species are still present
nrow(Coastal_Nest_Traits) == nrow(Coastal_Bodymass)

# how many species have nest safety values?
# for UAI?
Coastal_Nest_Traits %>% filter(!is.na(nest.safety)) %>% filter(!is.na(aveUAI)) %>% nrow()
# for MUTI?
Coastal_Nest_Traits %>% filter(!is.na(nest.safety)) %>% filter(!is.na(MUTIscore)) %>% nrow()
# for UN?
Coastal_Nest_Traits %>% filter(!is.na(nest.safety)) %>% filter(!is.na(Urban)) %>% nrow()

# save joined Nest traits and Coastal Species as an .rds file 
saveRDS(Coastal_Nest_Traits, here("Outputs", "Coastal_Species_Nest.rds"))


#############################################################################
######################### JOIN SEXUAL SELECTION TRAITS #######################
#############################################################################

# sexual selection traits: sexual dimorphism of plumage brightness, sexual dimorphism of plumage hue, 
# intensity of sexual selection Males (SSM), and Intensity of sexual selection on Females (SSF)
# SSM and SSM from - Delhey et al. (2023)
# Plumage sexual dimorphism (brightness and hue) from - Dunn et al. 2015


# # # # # # Strength of Sexual Selection for Males and Females (SSM and SSF) # # # # # #

# modify Delhey et al (2023) to get sexual selection traits - includes SSM and SSF. 

head(DELHEY) # this was imported above

# refine to have only sexual selection traits
DELHEY_ss <- DELHEY %>%
  mutate(Species_Jetz = str_replace(phylo, "_", " ")) %>% # rename sci names and remove underscore
  select(Species_Jetz, sex.sel.m, sex.sel.f) %>%
  distinct() # there are duplicate rows because many species have males and females. Remove them


# Since Delhey et. al (2023) uses BirdTree/Jetz taxonomy, let's join by Species_Jetz 

Coastal_SS_Traits <- left_join(Coastal_Bodymass, DELHEY_ss, by="Species_Jetz")

#view to check in on it 
View(Coastal_SS_Traits)
nrow(Coastal_SS_Traits) #807, as it should be 

#check to see if we have NAs for sexual selection intensity
Coastal_SS_Traits %>% filter(is.na(sex.sel.m)) %>% nrow()
Coastal_SS_Traits %>% filter(is.na(sex.sel.f)) %>% nrow()

# how many species have values for sex.sel for males ?
# for UAI?
Coastal_SS_Traits %>% filter(!is.na(sex.sel.m)) %>% filter(!is.na(aveUAI)) %>% nrow()
# for MUTI?
Coastal_SS_Traits %>% filter(!is.na(sex.sel.m)) %>% filter(!is.na(MUTIscore)) %>% nrow()
# for UN?
Coastal_SS_Traits %>% filter(!is.na(sex.sel.m)) %>% filter(!is.na(Urban)) %>% nrow()

# how many species have values for sex.sel for females ?
# for UAI?
Coastal_SS_Traits %>% filter(!is.na(sex.sel.f)) %>% filter(!is.na(aveUAI)) %>% nrow()
# for MUTI?
Coastal_SS_Traits %>% filter(!is.na(sex.sel.f)) %>% filter(!is.na(MUTIscore)) %>% nrow()
# for UN?
Coastal_SS_Traits %>% filter(!is.na(sex.sel.f)) %>% filter(!is.na(Urban)) %>% nrow()

# # # # # # Dichromatism HUE and BRIGHTNESS # # # # # #

# import Dunn et. al (2015) traits 

DUNN <- read.csv(here("Data", "2015DunnDichromatism.csv"))

head(DUNN)
colnames(DUNN)

# SumDiffPC1 = difference between M and F in plumage brightness 
# SumDiffPC2 = difference between M and F in plumage hue



#Dunn et al. 2015 uses Jetz taxonomy, so we will join by Species_Jetz
# rename and reformat column to be Species_Jetz so that it will join with Coastal_Bodymass 
DUNN.2 <- DUNN %>%
  mutate(Species_Jetz = str_replace(Species, "_", " "))

# rename hue and brightness columns
DUNN.3 <- DUNN.2 %>% 
  rename(Dichrom_bright = sumDiffPC1, # Change SummDiffPC1 to Dichrom_bright
        Dichrom_hue = sumDiffPC2) # Change SUmmDiffPC2 to Dichrom_hue

#select and keep only relevant columns 
DUNN.4 <- DUNN.3 %>%
  select (Species_Jetz, Dichrom_bright, Dichrom_hue)
colnames(DUNN.4)

# time to join the edited DUNN.4 with Coastal_SS_Traits 
Coastal_SS_Traits.2 <- left_join(Coastal_SS_Traits, DUNN.4)

#view to check in on it 
View(Coastal_SS_Traits.2)
nrow(Coastal_SS_Traits.2) #807, as it should be 

#check to see if we have NAs for sexual dichromatism and hue 
Coastal_SS_Traits.2 %>% filter(is.na(Dichrom_bright)) %>% nrow()
Coastal_SS_Traits.2 %>% filter(is.na(Dichrom_hue)) %>% nrow()

# how many species have dichromatism values?
# for UAI?
Coastal_SS_Traits.2 %>% filter(!is.na(Dichrom_bright)) %>% filter(!is.na(aveUAI)) %>% nrow()
# for MUTI?
Coastal_SS_Traits.2 %>% filter(!is.na(Dichrom_bright)) %>% filter(!is.na(MUTIscore)) %>% nrow()
# for UN?
Coastal_SS_Traits.2 %>% filter(!is.na(Dichrom_bright)) %>% filter(!is.na(Urban)) %>% nrow()

# how many species have hue values?
# for UAI?
Coastal_SS_Traits.2 %>% filter(!is.na(Dichrom_hue)) %>% filter(!is.na(aveUAI)) %>% nrow()
# for MUTI?
Coastal_SS_Traits.2 %>% filter(!is.na(Dichrom_hue)) %>% filter(!is.na(MUTIscore)) %>% nrow()
# for UN?
Coastal_SS_Traits.2 %>% filter(!is.na(Dichrom_hue)) %>% filter(!is.na(Urban)) %>% nrow()


# save joined Sexual selection traits and Coastal Species as an .rds file 

saveRDS(Coastal_SS_Traits.2, here("Outputs", "Coastal_Species_SSelect.rds"))

####################################################################
######################### JOIN SOCIAL TRAITS #######################
####################################################################
# social traits: territoriality, cooperative breeding 

# these traits are from Delhey et al. (2023), and are compiled from other original sources. 
# # # # # Original source for Territoriality - Tobias et al. (2016)
# # # # # Original source for Cooperative Breeding _ Cockburn (2006)

# modify Delhey et al (2023) to get territoriality and cooperative breeding
head(DELHEY) # this was imported above

# refine to have only social traits
DELHEY_social <- DELHEY %>%
  mutate(Species_Jetz = str_replace(phylo, "_", " ")) %>% # rename sci names and remove underscore
  select(Species_Jetz, territoriality, cooperative) %>%
  distinct() # there are duplicate rows because many species have males and females. Remove them


# Since Delhey et. al (2023) uses BirdTree/Jetz taxonomy, let's join by Species_Jetz 

Coastal_Social_Traits <- left_join(Coastal_Bodymass, DELHEY_social)

# view to check in on it 
View(Coastal_Social_Traits)
nrow(Coastal_Social_Traits) #807, as it should be 


# check to see if we have NAs for cooperative breeding and territoriality
Coastal_Social_Traits %>% filter(is.na(cooperative)) %>% nrow()
Coastal_Social_Traits %>% filter(is.na(territoriality)) %>% nrow()

# how many species have cooperative breeding scores?
# for UAI?
Coastal_Social_Traits %>% filter(!is.na(cooperative)) %>% filter(!is.na(aveUAI)) %>% group_by(cooperative) %>% count()
# for MUTI?
Coastal_Social_Traits %>% filter(!is.na(cooperative)) %>% filter(!is.na(MUTIscore)) %>% group_by(cooperative) %>% count()
# for UN?
Coastal_Social_Traits %>% filter(!is.na(cooperative)) %>% filter(!is.na(Urban)) %>% group_by(cooperative) %>% count()

# how many species have territoriality scores?
# for UAI?
Coastal_Social_Traits %>% filter(!is.na(territoriality)) %>% filter(!is.na(aveUAI)) %>% group_by(territoriality) %>% count()
# for MUTI?
Coastal_Social_Traits %>% filter(!is.na(territoriality)) %>% filter(!is.na(MUTIscore)) %>% group_by(territoriality) %>% count()
# for UN?
Coastal_Social_Traits %>% filter(!is.na(territoriality)) %>% filter(!is.na(Urban)) %>% group_by(territoriality) %>% count()


# save joined Social traits and Coastal Species as an .rds file 

saveRDS(Coastal_Social_Traits, here("Outputs", "Coastal_Species_Social.rds"))

