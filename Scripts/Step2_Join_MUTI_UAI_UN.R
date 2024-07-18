####################The objective of this script is to join all of our current indexes
#(MUTI, UAI, and UN) --> after we successfully join and get all scientific names to line up, 
#then we can run this file through our "retracing_coastal_species_sorting.r" script! 

#testing 
#testing 

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

#UAI Naming Scheme Methods --> weird combo of eBird and other naming schemes --> Sarah Investigated on UAINamestoSpecies.R

#MUTI Naming Scheme Methods --> seems to use eBird scheme --> uses eBird 6 letter codes, and AOS species checklist for taxonomic information. 
####There may be a few differences, but for the most part, eBird/Clements naming scheme and AOS agree. 

#UN Naming Scheme Methods --> ?? unstated naming scheme for bird species. 


#Plan --> try to join UAI and MUTI first... then try with UN (might be more difficult?)

#Load in MUTI! 

MUTI <- read.csv (here("Data", "Fanelli_Urban_Tolerance.csv"), header = T)
colnames(MUTI)
MUTI.r <- MUTI %>% rename(MUTIscore = PC1_scores) 
MUTI.f <- MUTI.r %>% select(MUTIscore, scientific_name, common_name)
colnames(MUTI.f)
View(MUTI.f)
nrow(MUTI.f)
#432 rows in MUTI (species)




#############################################################
##########let's take a look at Jetz and see how many we get 

MUTI.f$Species_Jetz <- MUTI.f$scientific_name
colnames(MUTI.f)


JetzMUTI_UAI <- inner_join(UAI, MUTI.f, by="Species_Jetz") 
head(JetzMUTI_UAI)
nrow(JetzMUTI_UAI)
#330.... 

nrow(MUTI.f)- nrow(JetzMUTI_UAI) 
#102 species off! prob not Jetz


#############################################################
##########let's take a look at BirdLife and see how many we get 

MUTI.f$Species_BirdLife <- MUTI.f$scientific_name
colnames(MUTI.f)


BirdLifeMUTI_UAI <- inner_join(UAI, MUTI.f, by="Species_BirdLife") 
nrow(BirdLifeMUTI_UAI)
#404.... 

nrow(MUTI.f)- nrow(BirdLifeMUTI_UAI) 
#28 species off! This is a possibility, but eBird may do better




#############################################################
##########let's take a look at eBird and see how many we get 

MUTI.f$Species_eBird <- MUTI.f$scientific_name
colnames(MUTI.f)

#test for how many species match UAI and MUTI via eBird/Clements naming scheme 

eBirdMUTI_UAI <- inner_join(UAI, MUTI.f, by="Species_eBird") 
nrow(eBirdMUTI_UAI)
#421.... 
colnames(eBirdMUTI_UAI)

nrow(MUTI.f)- nrow(eBirdMUTI_UAI) 
#11 species off! this is the most cohesive naming scheme for MUTI and UAI!

#let's left_join UAI and MUTI.f! this will add MUTI scores to all 421 species with UAI scientific name equivalent 
colnames(UAI)
colnames(MUTI.f)
UAI_and_MUTI <- MUTI.f %>% select(-Species_Jetz, -Species_BirdLife) %>% # drop these columns to avoid duplication when joining
  rename(CommonName_MUTI = common_name, SciName_MUTI = scientific_name) %>% # make columns labels more clearly attributable to either MUTI or UAI
  right_join(., UAI, by="Species_eBird") %>% # use right join to keep all rows in UAI and only those in MUTI.f that match UAI
  rename(CommonName_UAI = CommonName) # also rename column of common names from UAI, so it is clear where this column came from
  
nrow(UAI_and_MUTI) # 4347
View(UAI_and_MUTI)
colnames(UAI_and_MUTI)


# let's find out what 11 species are missing  and then manually join them! 

nomatch_eBird_MUTI <- anti_join (MUTI.f, eBirdMUTI_UAI, by = "Species_eBird")
nrow(nomatch_eBird_MUTI)
View(nomatch_eBird_MUTI)

#manual inspection 
#8 of the 11 not in UAI 
#3 of the 11 in UAI but not under eBird names, under Jetz names 
###### what to do about these two species? They are Brandt's Cormorant, Double-crested cormorant, and Pelagic Cormorant, so they ARE Going to be coastal. It would be nice to keep them. 

# since they are Jetz Taxonomy names... we need to do a join between nomatch_eBird_MUTI and UAI using Species_Jetz 
colnames(UAI_and_MUTI)
colnames(nomatch_eBird_MUTI)
nomatch_eBird_MUTI.r <- nomatch_eBird_MUTI %>% 
  rename(SciName_MUTI = scientific_name, CommonName_MUTI = common_name) %>%
  select (MUTIscore, SciName_MUTI, CommonName_MUTI, Species_Jetz) 
colnames(nomatch_eBird_MUTI.r)


# Emma - I simplified this section to do the same thing in fewer steps (SLJ)
# this keeps all the MUTI species, even the ones that do not have UAI scores as we may decide some are coastal 
# e.g., I think Northwestern crows sometimes forage on fish/intertidal inverts and they are in MUTI but not UAI

colnames(UAI)
colnames(UAI_and_MUTI)
UAI_and_MUTI_finalfew <- UAI %>% rename(CommonName_UAI = CommonName) %>%
  full_join(., nomatch_eBird_MUTI.r, by="Species_Jetz") %>% # use full join to match cormorants and keep all species from MUTI that are not in UAI
  filter(!is.na(MUTIscore)) %>% # keep only the rows where species have a MUTI score
  select(MUTIscore, SciName_MUTI, CommonName_MUTI, Species_eBird, Species_Jetz, Family_Jetz, Order_Jetz,
         CommonName_UAI, aveUAI, Order_eBird, Species_BirdLife) # reorder columns to match UAI_and_MUTI
  
View(UAI_and_MUTI_finalfew)
# here are the 11 species. The 3 cormorants joined to UAI using Jetz names.
# the 8 remaining species that are not found in UAI are also included but with no UAI score
# we need to bind these rows to the UAI_and_MUTI data frame to put everything together

# one last issue to resolve - the 3 cormorants are in both UAI_and_MUTI_finalfew and UAI_and_MUTI
# we only want to keep the version in UAI_and_MUTI_finalfew as this has the MUTI and UAI score for these species
# use anti_join to resolve this issue and then bind the two data frames together

UAI_and_MUTI_all <- anti_join(UAI_and_MUTI, UAI_and_MUTI_finalfew, by = join_by(Species_Jetz,CommonName_UAI, aveUAI)) %>%
  bind_rows(., UAI_and_MUTI_finalfew)
  
nrow(UAI_and_MUTI_all) # 4355
head(UAI_and_MUTI_all)

# perform some checks to make sure everything worked as expected

# how many rows/species have a MUTI score?
UAI_and_MUTI_all %>% filter(!is.na(MUTIscore)) %>% nrow()
# 432 which is 421 eBird species, the 3 Jetz cormorant species, and the 8 species from MUTI with no UAI match!

# how many rows/species have a UAI score?
UAI_and_MUTI_all %>% filter(!is.na(aveUAI)) %>% nrow() 
# 4347 species

# how many rows/species have both a UAI and MUTI score?
UAI_and_MUTI_all %>% filter(!is.na(aveUAI) & !is.na(MUTIscore)) %>% nrow() 
# 424 species

##################################################################
##################################################################
##################################################################
##################################################################
############# now it's time for UN Joining! ######################

UN <- read.csv (here("Data", "HuAndCardosoData.csv"), header = T)
View(UN)
UN <- UN %>% select(Species, Urban)
colnames(UN)
nrow(UN) # 528 rows -- we want to keep as many of these as possible! 

#create corresponding column names with UAI df 
UN$Species_eBird <- UN$Species
UN$Species_BirdLife <- UN$Species
UN$Species_Jetz <- UN$Species
colnames(UN)

###### Try with eBird ##################################
eBird_UN_MUTI_UAI <- UN %>% select(-Species_BirdLife, -Species_Jetz) %>%
  inner_join(UAI_and_MUTI_all, UN, by="Species_eBird") 
nrow(eBird_UN_MUTI_UAI)
#337 species --> not the best.... 

nrow(UN)- nrow(eBird_UN_MUTI_UAI) 
#191 species off! not looking like eBird


###### Try with BirdLife ##################################
BirdLife_UN_MUTI_UAI <- UN %>% select(-Species_eBird, -Species_Jetz) %>%
  inner_join(UAI_and_MUTI_all, ., by="Species_BirdLife") 
nrow(BirdLife_UN_MUTI_UAI)
#357 species --> not the best.... 

nrow(UN)- nrow(BirdLife_UN_MUTI_UAI) 
#171 species off! not looking like BirdLife

###### Try with Jetz ##################################

Jetz_UN_MUTI_UAI <- UN %>% select(-Species_eBird, -Species_BirdLife) %>%
  inner_join(UAI_and_MUTI_all, ., by="Species_Jetz") 
nrow(Jetz_UN_MUTI_UAI)
View(Jetz_UN_MUTI_UAI)
#479 species  

nrow(UN)- nrow(Jetz_UN_MUTI_UAI) 
#49 species off! it's prob jetz! 


################################################################
#Since Jetz has the most species match... let's join those 479 species to the df with UAI and MUTI 

#join 

UN_UAI_MUTI <- UN %>% select(-Species_eBird, -Species_BirdLife) %>%
  left_join(UAI_and_MUTI_all, .) 
nrow(UN_UAI_MUTI)
View(UN_UAI_MUTI)

####check how many species got added via this join... it should be 479

UN_UAI_MUTI$Urban <- ifelse(UN_UAI_MUTI$Urban == "U", 1, 0)
View(UN_UAI_MUTI)
#0 = N, 1 = U

UN_Check <- UN_UAI_MUTI %>% filter(!is.na(Urban))
nrow(UN_Check)
#479 species... thank goodness! this is correct  

#anti-join to find 49 species with UN who did not match with UAI using Jetz names
colnames(UN)
colnames(Jetz_UN_MUTI_UAI)
nomatch_Jetz_UN <- Jetz_UN_MUTI_UAI %>% select(Species, Species_Jetz) %>%
  anti_join (UN, .)
nrow(nomatch_Jetz_UN)
View(nomatch_Jetz_UN)

#why is this 80? according to earlier it should be 49... 
##
##
##
##
##THIS IS THE ISSUE (nomatch_Jetz_UN) CHECK BACK IN HERE ^

nrow(UAI)
#4347 - also this is a commit test in GitHub 


#see if there are eBird matches within the anti-join DF 


#see if there are BirdLife matches within the anit-join DF 