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
#28 species off! prob not Jetz




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
UAI_and_MUTI <- left_join(UAI, MUTI.f, by="Species_eBird") 
nrow(UAI_and_MUTI)
View(UAI_and_MUTI)
colnames(UAI_and_MUTI)


#let's find out what 11 species are missing  and then manually join them! 

nomatch_eBird_MUTI <- anti_join (MUTI.f, eBirdMUTI_UAI, by = "Species_eBird")
nrow(nomatch_eBird_MUTI)
View(nomatch_eBird_MUTI)
nomatch_eBird_MUTI$Species_Jetz <- nomatch_eBird_MUTI$scientific_name
colnames(nomatch_eBird_MUTI)
nomatch_eBird_MUTI$Species_Jetz.x <- nomatch_eBird_MUTI$Species_Jetz

#manual inspection 
#8 of the 11 not in UAI 
#3 of the 11 in UAI but not under eBird names, under Jetz names 
######what to do about these two species? They are Brandt's Cormorant, Double-crested cormorant, and Pelagic Cormorant, so they ARE Going to be coastal. It would be nice to keep them. 

#since they are Jetz Taxonomy names... what if we do a join between nomatch_eBird_MUTI and UAI using Species_Jetz? 
colnames(UAI_and_MUTI)
colnames(nomatch_eBird_MUTI)
nomatch_eBird_MUTI.r <- nomatch_eBird_MUTI %>% 
  select (MUTIscore, Species_Jetz) 
nomatch_eBird_MUTI.r$Species_Jetz.x <- nomatch_eBird_MUTI.r$Species_Jetz
colnames(nomatch_eBird_MUTI.r)


#FROM CHAT 
#this join maneuver works! 

UAI_and_MUTI <- UAI_and_MUTI %>%
  left_join(nomatch_eBird_MUTI %>% select(Species_Jetz.x, MUTIscore), by = "Species_Jetz.x") %>%
  mutate(MUTIscore = coalesce(MUTIscore.x, MUTIscore.y)) %>%
  select(-MUTIscore.x, -MUTIscore.y)
colnames(UAI_and_MUTI)

UAI_and_MUTI2 <- left_join(UAI_and_MUTI, nomatch_eBird_MUTI.r, by="Species_Jetz.x") 
nrow(UAI_and_MUTI2)
View(UAI_and_MUTI2)
colnames(UAI_and_MUTI2)
#MUTIscore.x is the column we want to look at here... we don't care about MUTIscore.y 

#lets select to retain only MUTIscore.x 

UAI_and_MUTI3 <- UAI_and_MUTI2 %>% select (Species_Jetz.x, Family_Jetz, Order_Jetz, CommonName, aveUAI, Species_eBird, Species_eBird, Order_eBird, Species_BirdLife.x, scientific_name, common_name, MUTIscore.x)
View(UAI_and_MUTI3)

#manually checked for Pelagic Cormorant, Double-Crested Cormorant, and Brandt's cormorant, and those have been successfully joined!
#now every MUTI that has a corresponding UAI species is joined! Yay! 


#QUESTION - there are 8 species from MUTI that are not represented by UAI... currently not joined... 
#should they be? 

MUTI_Check <- UAI_and_MUTI3 %>% filter(!is.na(MUTIscore.x))
nrow(MUTI_Check)
#424! which is 421 eBird species and the 3 Jetz cormorant species! 


##################################################################
##################################################################
##################################################################
##################################################################
#############now it's time for UN Joining! ######################

UN <- read.csv (here("Data", "HuAndCardosoData.csv"), header = T)
View(UN)
UN <- UN %>% select(Species, Urban)
colnames(UN)
nrow(UN) #528 rows -- we want to keep as many of these as possible! 

#create corresponding column names with UAI df 
UN$Species_eBird <- UN$Species
UN$Species_BirdLife <- UN$Species
UN$Species_Jetz <- UN$Species
colnames(UN)

######Try with eBird##################################
eBird_UN_MUTI_UAI <- inner_join(UAI, UN, by="Species_eBird") 
nrow(eBird_UN_MUTI_UAI)
#337 species --> not the best.... 

nrow(UN)- nrow(eBird_UN_MUTI_UAI) 
#191 species off! not looking like eBird


######Try with BirdLife##################################
BirdLife_UN_MUTI_UAI <- inner_join(UAI, UN, by="Species_BirdLife") 
nrow(BirdLife_UN_MUTI_UAI)
#357 species --> not the best.... 

nrow(UN)- nrow(BirdLife_UN_MUTI_UAI) 
#171 species off! not looking like BirdLife

######Try with Jetz##################################

Jetz_UN_MUTI_UAI <- inner_join(UAI, UN, by="Species_Jetz") 
nrow(Jetz_UN_MUTI_UAI)
View(Jetz_UN_MUTI_UAI)
#479 species  

nrow(UN)- nrow(Jetz_UN_MUTI_UAI) 
#49 species off! it's prob jetz! 

colnames(Jetz_UN_MUTI_UAI)
Jetz_UN_MUTI_UAI$Species_Jetz.x <- Jetz_UN_MUTI_UAI$Species_Jetz
colnames(Jetz_UN_MUTI_UAI)


################################################################
#Since Jetz has the most species match... let's join those 479 species to the df with UAI and MUTI 

#join 
colnames(UN_UAI_MUTI)
UN$Species_Jetz.x <- UN$Species_Jetz
UN_UAI_MUTI <- left_join(UAI_and_MUTI3, UN, by="Species_Jetz.x") 
nrow(UN_UAI_MUTI)
View(UN_UAI_MUTI)

####check how many species got added via this join... it should be 479

UN_UAI_MUTI$Urban <- ifelse(UN_UAI_MUTI$Urban == "U", 1, 0)
View(UN_UAI_MUTI)
#0 = N, 1 = U

UN_Check <- UN_UAI_MUTI %>% filter(!is.na(Urban))
nrow(UN_Check)
#479 species... thank goodness! this is correct  

#anti-join 
colnames(UN)
colnames(Jetz_UN_MUTI_UAI)
nomatch_Jetz_UN <- anti_join (UN, Jetz_UN_MUTI_UAI, by = "Species_Jetz.x")
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