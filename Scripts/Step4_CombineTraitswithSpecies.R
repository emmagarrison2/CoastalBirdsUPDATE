#The objective of this script is to join Coastal_Birds_list.rds to trait data. 
#We will create a joined dataframe for each category of predictor trait variables, as follows... 
######SENSORY TRAITS 
######DIET TRAITS - done
######NESTING TRAITS - done 
######LIFE HISTORY TRAITS 
######SEXUAL SELECTION TRAITS 
######SOCIAL TRAITS 



# Load packages
library(here)
library(tidyverse)
if(!require(fuzzyjoin)){
  install.packages("fuzzyjoin")
  require(fuzzyjoin)
}
library(fuzzyjoin)
library(stringi)



#read in Coastal Bird List 

Coastal_Birds_list <- readRDS (here("Data", "Coastal_Birds_list.rds"))
View(Coastal_Birds_list)
nrow(Coastal_Birds_list)
#827 total coastal species 

#let's see how many birds have UAI scores, MUTI scores, and UN scores 

Coastal_birds_UAI <- Coastal_Birds_list %>% 
  filter(!is.na(aveUAI))
nrow(Coastal_birds_UAI)
#818 species 

Coastal_birds_MUTI <- Coastal_Birds_list %>% 
  filter(!is.na(MUTIscore))
nrow(Coastal_birds_MUTI)
#133 species 

Coastal_birds_UN <- Coastal_Birds_list %>% 
  filter(!is.na(Urban))
nrow(Coastal_birds_UN)
#132 species 

################# JOIN BODY MASS #################
# body mass from avonet (Tobias et al. 2022 Ecol Letters)

######
#####
## SARAH --- 
## WE USED ELTON.CSV FOR BODYMASS.VALUE LAST TIME WE DID THE JOINS... WHY ARE WE USING AVONET // AND IS AVONET BODY MASS (MASS) BETTER? 
######
######

# this join between coastal species and body mass will be used as the starting df to add other predictor traits 

AVONET <- read.csv(here("Data", "avonet.csv"))
head(AVONET)

head(Coastal_Birds_list)

#we are looking to join body mass (Mass (g)) with Coastal_Birds_list based on scientific name for species. 
#Avonet uses Bird Life names, so join via Coastal_Birds_list column Species_BirdLife

AVONET.2 <- AVONET %>% 
  rename(Species_BirdLife = Species1) 

AVONET.3 <- AVONET.2 %>% 
  select(Species_BirdLife, Mass)

colnames(AVONET.3)
head(AVONET.3)

Coastal_Bodymass <- left_join(Coastal_Birds_list, AVONET.3)
head(Coastal_Bodymass)
nrow(Coastal_Bodymass) #827, looks good. 

#check to see if we have NAs for body mass (Mass)

Mass_test <- Coastal_Bodymass %>% 
  filter(!is.na(Mass))
nrow(Mass_test)
#all 827 species have a body mass value from Avonet 

#save Coastal_Bodymass as .rds for easy retrieval 

saveRDS(Coastal_Bodymass, here("Outputs", "Coastal_Species_w_Mass.rds"))



######################### JOIN SENSORY TRAITS ####################### 
# sensory traits: dim light vision and dominant vocal frequency
## Why ARE WE USING HU & CARDOSO 2009 FOR VOCAL FREQUENCY?

####### Dim light vision (eye CT ratio) ######
# uses BirdTree/Jetz taxonomy

####### Dominant vocal frequency ######
# Hu and Cardoso 2009


######################### JOIN DIET TRAITS #######################
# diet traits: % invertebrates, % vertebrates, % plant/seed, % fruit/nectar
# from Wilman et al. 2014 (elton traits)

#read in Coastal_Bodymass  

Coastal_Bodymass <- readRDS(here("Outputs", "Coastal_Species_w_Mass.rds"))

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

Coastal_Diet_Traits <- left_join(Coastal_Bodymass, ELTON_DIET3)

#view to check in on it 
View(Coastal_Diet_Traits)
nrow(Coastal_Diet_Traits) #827, as it should be 


#check to see if we have NAs for Diet traits (test using Diet.Inv)

Diet_test <- Coastal_Diet_Traits %>% 
  filter(!is.na(Diet.Inv))
nrow(Diet_test)
#all 827 species have diet traits from Wilman et al. 2014 (elton traits)

#Simplify the Diet Traits 

######################################################################################
######################################################################################
#####################time start binning our categories for diet! 4 new bins --> vert, invert, plant/seed, and fruit/nectar (which we will test, but may have to drop)
#First let's combine diet of vertebrate endotherms, vertbrate ectotherms, vertebrate unknown, vertebrate fish, and scavenging together to make a genreal "Diet Vertebrate" bin 

Coastal_Diet_Traits$Diet.Vend <- as.numeric(Coastal_Diet_Traits$Diet.Vend)
Coastal_Diet_Traits$Diet.Vect <- as.numeric(Coastal_Diet_Traits$Diet.Vect)
Coastal_Diet_Traits$Diet.Vunk <- as.numeric(Coastal_Diet_Traits$Diet.Vunk)
Coastal_Diet_Traits$Diet.Vfish <- as.numeric(Coastal_Diet_Traits$Diet.Vfish)
Coastal_Diet_Traits$Diet.Scav <- as.numeric(Coastal_Diet_Traits$Diet.Scav)

Coastal_Diet_Traits$Diet.Vert <- Coastal_Diet_Traits$Diet.Vend + Coastal_Diet_Traits$Diet.Vect + Coastal_Diet_Traits$Diet.Vunk + Coastal_Diet_Traits$Diet.Vfish + Coastal_Diet_Traits$Diet.Scav
head(Coastal_Diet_Traits)

#remove columns that we used to make Diet.Vert

Coastal_Diet_Traits2 <- Coastal_Diet_Traits %>% select (-Diet.Vend, -Diet.Vect, -Diet.Vunk, -Diet.Vfish, - Diet.Scav)
colnames(Coastal_Diet_Traits2)

#
#same thing for combining fruit and nectar combined to "Diet Fruit / Nectar"
# 

Coastal_Diet_Traits2$Diet.Fruit <- as.numeric(Coastal_Diet_Traits2$Diet.Fruit)
Coastal_Diet_Traits2$Diet.Nect <- as.numeric(Coastal_Diet_Traits2$Diet.Nect)

Coastal_Diet_Traits2$Diet.FN <- Coastal_Diet_Traits2$Diet.Fruit + Coastal_Diet_Traits2$Diet.Nect 
head(Coastal_Diet_Traits2$Diet.Fruit)
head(Coastal_Diet_Traits2$Diet.Nect)
head(Coastal_Diet_Traits2$Diet.FN)


#remove columns that we used to make Diet.FN

Coastal_Diet_Traits3 <- Coastal_Diet_Traits2 %>% select (-Diet.Fruit, -Diet.Nect)
colnames(Coastal_Diet_Traits3)

#
#same thing for combining fruit and nectar combined to "Diet Fruit / Nectar"
#

Coastal_Diet_Traits3$Diet.PlantO <- as.numeric(Coastal_Diet_Traits3$Diet.PlantO)
Coastal_Diet_Traits3$Diet.Seed <- as.numeric(Coastal_Diet_Traits3$Diet.Seed)

Coastal_Diet_Traits3$Diet.PS <- Coastal_Diet_Traits3$Diet.PlantO + Coastal_Diet_Traits3$Diet.Seed 
head(Coastal_Diet_Traits3$Diet.PlantO)
head(Coastal_Diet_Traits3$Diet.Seed)
head(Coastal_Diet_Traits3$Diet.PS)


#remove columns that we used to make Diet.FN

Coastal_Diet_Traits4 <- Coastal_Diet_Traits3 %>% select (-Diet.PlantO, -Diet.Seed)
colnames(Coastal_Diet_Traits4)


#save Coastal_Diet_Traits as .rds for easy retrieval 

saveRDS(Coastal_Diet_Traits4, here("Outputs", "Coastal_Species_Diet.rds"))




######################### JOIN LIFE HISTORY TRAITS #######################
# life history traits include: clutch size, longevity, brood value, developmental mode

####### Clutch Size #######
# from Myhrvold
# may need to use both BirdLife and BirdTree/Jetz

####### Longevity #######
# from Bird et al. 2020
# uses BirdLife

####### Brood Value #######
# calculate most species manually using traits in Bird and Myhrvold
# fix one sketchy species
# fill in additional from Sayol et al. 2020. I think this works with Jetz but need to check


####### Developmental Mode #######
# from Wang and Kimball 2016
# which taxonomy?




######################### JOIN NEST TRAITS #######################
# nest traits include: nest site (low, high), nest structure (open or enclosed), and nest safety

#read in Coastal_Bodymass  

Coastal_Bodymass <- readRDS(here("Outputs", "Coastal_Species_w_Mass.rds"))

# import nest site and structure info from Chia et al 2023
# uses the BirdLife taxonomy

CHIA_NEST <- read.csv(here("Data", "nests.csv"))

head(CHIA_NEST)

nests_names <- CHIA_NEST %>%
  rename(Species_BirdLife = Scientific_name)

Coastal_Nest_Traits <- left_join(Coastal_Bodymass, nests_names)


##### simplify the nest traits

# # # # NEST STRATEGY # # # # 
#ok now lets sort our nest types in the bins "open" and "enclosed" (NestStr = Nest Strategy)

#open 
Coastal_Nest_Traits <- mutate(Coastal_Nest_Traits, NestStr_Open = ifelse(NestStr_scrape > 0 | NestStr_platform > 0 | NestStr_cup > 0, 1, 0))
summary(Coastal_Nest_Traits$NestStr_Open)
colnames(Coastal_Nest_Traits)

#enclosed
Coastal_Nest_Traits <- mutate(Coastal_Nest_Traits, NestStr_Enclosed = ifelse(NestStr_dome > 0 | NestStr_dome_tunnel > 0 | NestStr_primary_cavity > 0 | NestStr_second_cavity > 0, 1, 0))
summary(Coastal_Nest_Traits$NestStr_Enclosed)
colnames(Coastal_Nest_Traits)


# # # # NEST SITE # # # #  
#ok now lets sort our nest sites in the bins "low" (on or underground) and "high" (above the ground)

#low 
Coastal_Nest_Traits <- mutate(Coastal_Nest_Traits, NestSite_low = ifelse(NestSite_ground > 0 | NestSite_underground > 0 | NestSite_waterbody > 0, 1, 0))
summary(Coastal_Nest_Traits$NestSite_low)
colnames(Coastal_Nest_Traits)

#high
Coastal_Nest_Traits <- mutate(Coastal_Nest_Traits, NestSite_high = ifelse(NestSite_tree > 0 | NestSite_nontree > 0 | NestSite_cliff_bank > 0 | NestSite_termite_ant > 0, 1, 0))
summary(Coastal_Nest_Traits$NestSite_high)
colnames(Coastal_Nest_Traits)



#remove all unnecessary columns 

Coastal_Nest_Traits2 <- Coastal_Nest_Traits %>% 
  select(-Seq_HBWBLv5, -Order, -Family, -Common_name, -SISRecID, -AvibaseID, -Parasite, -Mound, -NestSite_ground, -NestSite_tree, -NestSite_nontree, 
         -NestSite_cliff_bank, -NestSite_underground, -NestSite_waterbody, -NestSite_termite_ant, 
         -NestStr_scrape, -NestStr_platform, -NestStr_cup, -NestStr_dome, -NestStr_dome_tunnel, -NestStr_primary_cavity, -NestStr_second_cavity, -NestAtt_basal, 
         -NestAtt_forked, -NestAtt_lateral, -NestAtt_pensile, -Ref_site, -Ref_str, -Ref_att, -Identifier, -Retrieval_time)

colnames(Coastal_Nest_Traits2)

######
######Nest Safety from Delhey et al. (2023)

DELHEY <- read.csv(here("Data", "Delhey_2023_DS7.csv"))

head(DELHEY)

#rename the column ph
DELHEY.2 <- DELHEY %>%
  rename(Species_Jetz = phylo)

#reformat the DELHEY.2$Species_Jetz column 

DELHEY.3 <- DELHEY.2 %>% 
  mutate(Species_Jetz = str_replace_all(Species_Jetz, "_", " "),  
         Species_Jetz = str_to_sentence(Species_Jetz))
head(DELHEY.3)

#reduce the columns down to only relevant ones 
colnames(DELHEY.3)

#INCLUDED COOPERATIVE BREEDING AND TERRITORIALITY FOR NOW 
DELHEY.4 <- DELHEY.3 %>% 
  dplyr::select(Species_Jetz, sex.sel.m, sex.sel.f, nest.safety, developmental_mode, cooperative, territoriality) %>% 
  distinct(Species_Jetz, .keep_all = TRUE)

View(DELHEY.4)

saveRDS(DELHEY.4, here("OUtputs", "Delhey_refined_traits.rds"))


Coastal_Nest_Traits3 <- left_join(Coastal_Nest_Traits2, DELHEY.4)

View(Coastal_Nest_Traits3)
nrow(Coastal_Nest_Traits3) #827, as it should be 

#let's check to see how many Coastal species have some of these categorical variables... such as Developmental mode! 
Delhey_test <- Coastal_Nest_Traits3 %>% 
  filter(!is.na(developmental_mode))
nrow(Delhey_test)
#795

#save joined Nest traits and Coastal Species as an .rds file 

saveRDS(Coastal_Nest_Traits3, here("Outputs", "Coastal_Species_Nest.rds"))






######################### JOIN SEXUAL SELECTION TRAITS #######################
# sexual selection traits: sexual dimorphism of plumage brightness, sexual dimorphism of plumage hue, 
# intensity of sexual selection Males (SSM), and Intensity of sexual selection on Females (SSF)

# SSM and SSM from - Delhey et al. (2014)
# Plumage sexual dimorphism (brightness and hue) from - Dunn et al. 2015


#read in Coastal_Bodymass  

Coastal_Bodymass <- readRDS(here("Outputs", "Coastal_Species_w_Mass.rds"))

############## SSM and SSF

# import Delhey et al (2023) traits - includes SSM and SSF. Previously edited due to use of this data source for other trait categories. 

DELHEY_ss <- readRDS(here("Outputs", "Delhey_refined_traits.rds"))

head(DELHEY_ss)


# Since Delhey et. al (2023) uses BirdTree/Jetz taxonomy, let's join by Species_Jetz 

Coastal_SS_Traits <- left_join(Coastal_Bodymass, DELHEY_ss)

#view to check in on it 
View(Coastal_SS_Traits)
nrow(Coastal_SS_Traits) #827, as it should be 


#check to see if we have NAs for sexual selection intensity (test using sex.sel.m)

SSM_test <- Coastal_SS_Traits %>% 
  filter(!is.na(sex.sel.m))
nrow(SSM_test)
#795 coastal species have a value for sexual selection intensity on males (SSM) from Delhey et al. 2023 


############## Dichromatism HUE and BRIGHTNESS 

# import Dunn et. al (2015) traits 

DUNN <- read.csv(here("Data", "2015DunnDichromatism.csv"))

head(DUNN)

# SumDiffPC1 = difference between M and F in plumage  brightness 
# SumDiffPC2 = difference between M and F in plumage hue

colnames(DUNN)

#Dunn et al. 2015 uses JETZ taxonomy, so we will join by Species_BirdLife

#Change 

DUNN.2 <- DUNN %>%
  rename(Species_Jetz = Species) 

# reformat Species_BirdLife so that it will join with Coastal_Bodymass 

DUNN.3 <- DUNN.2 %>% 
  mutate(Species_Jetz = str_replace_all(Species_Jetz, "_", " "),  
         Species_Jetz = str_to_sentence(Species_Jetz))

head(DUNN.3)

# Change SummDiffPC1 to Dichrom_bright 

DUNN.4 <- DUNN.3 %>% 
  rename(Dichrom_bright = sumDiffPC1)

# Change SUmmDiffPC2 to Dichrom_hue

DUNN.5 <- DUNN.4 %>% 
  rename(Dichrom_hue = sumDiffPC2)

colnames(DUNN.5)

#select and keep only relevant columns 

DUNN.6 <- DUNN.5 %>%
  select (Species_Jetz, Dichrom_bright, Dichrom_hue)
colnames(DUNN.6)

#time to join the edited DUNN.6 with Coastal_SS_Traits 

Coastal_SS_Traits.2 <- left_join(Coastal_SS_Traits, DUNN.6)

#view to check in on it 
View(Coastal_SS_Traits.2)
nrow(Coastal_SS_Traits.2) #827, as it should be 


#check to see if we have NAs for sexual dichromatism (test using Dichrom_bright)

dichrom_test <- Coastal_SS_Traits.2 %>% 
  filter(!is.na(Dichrom_bright))
nrow(dichrom_test)
#208 coastal species have a value for sexual dichromatism of plumage brightness (Dichrom_bright) from Dunn et al. 2015


#save joined Sexual selection traits and Coastal Species as an .rds file 

saveRDS(Coastal_SS_Traits.2, here("Outputs", "Coastal_Species_SSelect.rds"))




######################### JOIN SOCIAL TRAITS #######################
# social traits: territoriality, cooperative breeding 

# these traits are from Delhey et al. (2014), and are compiled from other original sources. 
# # # # #Original source for Territoriality - Tobias et al. (2016)
# # # # #Original source for Cooperative Breeding _ Cockburn (2006)


#read in Coastal_Bodymass  

Coastal_Bodymass <- readRDS(here("Outputs", "Coastal_Species_w_Mass.rds"))


# import Delhey et al (2023) traits - includes territoriality and cooperative breeding.
#Previously edited due to use of this data source for other trait categories. 

DELHEY_social <- readRDS(here("Outputs", "Delhey_refined_traits.rds"))

head(DELHEY_social)


# Since Delhey et. al (2023) uses BirdTree/Jetz taxonomy, let's join by Species_Jetz 

Coastal_Social_Traits <- left_join(Coastal_Bodymass, DELHEY_social)

#view to check in on it 
View(Coastal_Social_Traits)
nrow(Coastal_Social_Traits) #827, as it should be 


#check to see if we have NAs for sexual selection intensity (test using sex.sel.m)

coop_test <- Coastal_Social_Traits %>% 
  filter(!is.na(cooperative))
nrow(coop_test)
#795 coastal species have a value for territoriality from Delhey et al. (2023) 


#save joined Social traits and Coastal Species as an .rds file 

saveRDS(Coastal_Social_Traits, here("Outputs", "Coastal_Species_Social.rds"))

