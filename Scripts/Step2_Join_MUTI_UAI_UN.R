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
UN <- UN %>% select(Species, Urban) %>% rename(Species_UN = Species)
colnames(UN)
nrow(UN) # 528 rows -- we want to keep as many of these as possible! 

#create corresponding column names with UAI df 
UN$Species_eBird <- UN$Species_UN
UN$Species_BirdLife <- UN$Species_UN
UN$Species_Jetz <- UN$Species_UN
colnames(UN)

###### Try with eBird ##################################
# it is best to try and match with eBird as all the species names are unique
eBird_UN_MUTI_UAI <- UN %>% select(-Species_BirdLife, -Species_Jetz) %>%
  left_join(UAI_and_MUTI_all, ., by="Species_eBird") 
nrow(eBird_UN_MUTI_UAI)

# how many species matched?
eBird_UN_MUTI_UAI_match <- eBird_UN_MUTI_UAI %>% filter(!is.na(Species_UN)) 
nrow(eBird_UN_MUTI_UAI_match) # 337
# are they all unique? This should match the number output by the row of code above
length(unique(eBird_UN_MUTI_UAI_match$Species_UN))

# how many species did not match?
nrow(UN)- nrow(eBird_UN_MUTI_UAI_match) 
#191 species off

# get the unmatched UN species in a data frame
colnames(eBird_UN_MUTI_UAI_match)
colnames(UN)
UN_nomatch <- eBird_UN_MUTI_UAI_match %>%
  select(Species_UN, Urban) %>%
  anti_join(UN, .)
nrow(UN_nomatch) # 191

###### Try matching the 191 remaining species with Jetz ##############

Jetz_UN_MUTI_UAI <- UN_nomatch %>% select(-Species_eBird, -Species_BirdLife) %>%
  inner_join(UAI_and_MUTI_all, ., by="Species_Jetz")
View(Jetz_UN_MUTI_UAI)

# how many UAI/MUTI species matched with UN using Jetz ?
nrow(Jetz_UN_MUTI_UAI)
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
# if we look at the UN paper can we determine if the authors were studying one of these two species, or possibly both clumped under one name?


# which UN species still do not have a match with UAI?
UN_stillnomatch <- Jetz_UN_MUTI_UAI %>% distinct(Species_UN) %>%
  anti_join(UN_nomatch, .)
nrow(UN_stillnomatch) # 80 species
# we could try to match using BirdLife. I tried this but it did not help, so I deleted these steps


##### Put matches together and build towards final data frame ######

# combine Jetz and eBird matches into one using bind_rows
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
nrow(UN_combine2) # 4355

# get 80 UN species with no matches and add those
UAI_MUTI_UN_final <- full_join(UN_combine2, UN_stillnomatch)
  
nrow(UAI_MUTI_UN_final) # should equal 4355 + 80 

##### Perform some Final Checks ####

# how many species have UAI score?
UAI_MUTI_UN_final %>% filter(!is.na(aveUAI)) %>% nrow() # 4347

# how many species have MUTI score?
UAI_MUTI_UN_final %>% filter(!is.na(MUTIscore)) %>% nrow() # 432

# how many species have Urban score?
UAI_MUTI_UN_final %>% filter(!is.na(Urban)) %>% nrow() # 533
# how many are distinct species names (this should equal the number in the original data frame)?
UAI_MUTI_UN_final %>% filter(!is.na(Urban)) %>% distinct(Species_UN) %>% nrow() # 528

# how many species have UN but no MUTI or UAI?
UAI_MUTI_UN_final %>% filter(!is.na(Urban)) %>% filter(is.na(MUTIscore) & is.na(aveUAI)) %>% nrow() # 80

# how many species have MUTI but no UAI or UN?
UAI_MUTI_UN_final %>% filter(!is.na(MUTIscore)) %>% filter(is.na(Urban) & is.na(aveUAI)) %>% nrow() # 8

##### Export list of species and urban tolerance indices ######
saveRDS (UAI_MUTI_UN_final, here("Outputs", "UAI_MUTI_UN_final.rds"))
