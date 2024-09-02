# The objective of this script is to find Phylogenetic Signal (lambda) for each 
# continuous predictor trait, and numeric ordinal predictor traits. 
# (Unable to use picante to find phylo signal for binomial predictor traits). 
# This will require importing joined rds files of Coastal sp, index values, and predictor trait groups.

# This value is important for... ?? ensuring that it is appropriate to use phylo linear model, 
# instead of regular linear model 

#load/install packages 

library(here)
library(tidyverse)
library(ape)
library(geiger)
library(picante)
library(phytools)


######################## Phylogenetic Signal (lambda) for SENSORY TRAITS ######################## 

#import Coastal_Species_Sensory.rds 

Sensory_Traits <- readRDS(here("Outputs", "Coastal_Species_Sensory.rds"))
head(Sensory_Traits)
colnames(Sensory_Traits)

#should I be going off of... Species_Jetz, Species_BirdLife, or Species_eBird???? prob Jetz because I join it to Jetz tree.... 

#do some reformatting and re-naming... 
Sensory_Traits <- Sensory_Traits %>%
  rename(Species = Species_Jetz)
colnames(Sensory_Traits)
#nice 

#now, reformat Species column so that it is in Aaaa_aaaa format! 
Sensory_Traits <- Sensory_Traits %>%
  mutate(Species = str_replace(Species, " ", "_"))
head(Sensory_Traits)

### we need to trim Jetz tree that we imported above to get a tree that only contains coastal species

# create data frame with traits and species as the rownames
Sensoryspp <- Sensory_Traits %>%
  column_to_rownames(var="Species")

### load phylogenetic tree
jetztree <- read.tree(here("Data", "Jetz_ConsensusPhy.tre")) 
jetztree$tip.label # look at species name formatting for tree tips to confirm it matches formatting we are using - yes!

# check that the tree is ultrametric (do all the tree tips line up?)
is.ultrametric(jetztree)
#true! 


# check that all the species in our data are also present in the tree
check_sensory <- name.check(jetztree, Sensoryspp) # there should be many more species in the tree than our data
check_sensory

# trim the tree to match the data by dropping all the extra species identified in the previous step
jetztree_sensory <-drop.tip(jetztree, check_sensory$tree_not_data)

# check again whether tree has same species as data. Should say "OK"
name.check(jetztree_sensory, Sensoryspp)
#OK 

##########
colnames(Sensoryspp)

## for each trait (C.T and peak_freq) we need to create a vector of values with species as rownames

      #CT 

# look at distribution
hist(Sensoryspp$C.T) # not skewed. no transformation needed

# get vector of C.T values with species names as rownames
CT_vect <- setNames(Sensoryspp$C.T, rownames(Sensoryspp))

# drop species with NA values
CT_vect[!is.na(CT_vect)]

# get measure of lambda for CT
# will give message about dropping species. This is okay. There are species in our tree that do not have CT values and the function is removing those
CT_lambda2 <-phylosig(tree = jetztree_sensory, CT_vect, method="lambda", test=T) # inputs are trimmed phylogenetic tree and vector of body mass
CT_lambda2 # lambda = 0.562 , p << 0.001
plot.phylosig(CT_lambda2) # plot the likelihood surface for lambda



      #peak_freq

# look at distribution
hist(Sensoryspp$peak_freq) # not skewed. no transformation needed

# get vector of peak frequency values with species names as rownames
pf_vect <- setNames(Sensoryspp$peak_freq, rownames(Sensoryspp))

# drop species with NA values
pf_vect[!is.na(pf_vect)]

# get measure of lambda for peak frequency 
# will give message about dropping species. This is okay. There are species in our tree that do not have CT values and the function is removing those
pf_lambda2 <-phylosig(tree = jetztree_sensory, pf_vect, method="lambda", test=T) # inputs are trimmed phylogenetic tree and vector of body mass
pf_lambda2 # lambda = 0.936 , p << 0.001
plot.phylosig(pf_lambda2) # plot the likelihood surface for lambda


######################## Phylogenetic Signal (lambda) for DIET TRAITS ######################## 


#import Coastal_Species_Sensory.rds 

Diet_Traits <- readRDS(here("Outputs", "Coastal_Species_Diet.rds"))
head(Diet_Traits)
colnames(Diet_Traits)

#should I be going off of... Species_Jetz, Species_BirdLife, or Species_eBird???? prob Jetz because I join it to Jetz tree.... 

#do some reformatting and re-naming... 
Diet_Traits <- Diet_Traits %>%
  rename(Species = Species_Jetz)
colnames(Diet_Traits)
#nice 

#now, reformat Species column so that it is in Aaaa_aaaa format! 
Diet_Traits <- Diet_Traits %>%
  mutate(Species = str_replace(Species, " ", "_"))
View(Diet_Traits)

### we need to trim Jetz tree that we imported above to get a tree that only contains coastal species

# create data frame with traits and species as the rownames
Dietspp <- Diet_Traits %>%
  column_to_rownames(var="Species")

### load phylogenetic tree
jetztree <- read.tree(here("Data", "Jetz_ConsensusPhy.tre")) 
jetztree$tip.label # look at species name formatting for tree tips to confirm it matches formatting we are using - yes!

# check that the tree is ultrametric (do all the tree tips line up?)
is.ultrametric(jetztree)
#true! 


# check that all the species in our data are also present in the tree
check_diet <- name.check(jetztree, Dietspp) # there should be many more species in the tree than our data
check_diet

# trim the tree to match the data by dropping all the extra species identified in the previous step
jetztree_diet <-drop.tip(jetztree, check_diet$tree_not_data)

# check again whether tree has same species as data. Should say "OK"
name.check(jetztree_diet, Dietspp)
#OK 

##########
colnames(Dietspp)

     # % Diet Invert 

# look at distribution
hist(Dietspp$Diet.Inv) # not skewed. definitely not normal, but no log transformation needed

# get vector of C.T values with species names as rownames
invert_vect <- setNames(Dietspp$Diet.Inv, rownames(Dietspp))

# drop species with NA values
invert_vect[!is.na(invert_vect)]

# get measure of lambda for CT
# will give message about dropping species. This is okay. There are species in our tree that do not have CT values and the function is removing those
invert_lambda2 <-phylosig(tree = jetztree_diet, invert_vect, method="lambda", test=T) # inputs are trimmed phylogenetic tree and vector of body mass
invert_lambda2 # lambda = 0.845 , p << 0.001
plot.phylosig(invert_lambda2) # plot the likelihood surface for lambda


     # % Diet Vert 

# look at distribution
hist(Dietspp$Diet.Vert) # not skewed. definitely not normal, but no log transformation needed

# get vector of C.T values with species names as rownames
vert_vect <- setNames(Dietspp$Diet.Vert, rownames(Dietspp))

# drop species with NA values
vert_vect[!is.na(vert_vect)]

# get measure of lambda for CT
# will give message about dropping species. This is okay. There are species in our tree that do not have CT values and the function is removing those
vert_lambda2 <-phylosig(tree = jetztree_diet, vert_vect, method="lambda", test=T) # inputs are trimmed phylogenetic tree and vector of body mass
vert_lambda2 # lambda = 0.860 , p << 0.001
plot.phylosig(vert_lambda2) # plot the likelihood surface for lambda



     # % Diet Fruit/Nectar 


# look at distribution
hist(Dietspp$Diet.FN) # not skewed. definitely not normal, but no log transformation needed

# get vector of C.T values with species names as rownames
FN_vect <- setNames(Dietspp$Diet.FN, rownames(Dietspp))

# drop species with NA values
FN_vect[!is.na(FN_vect)]

# get measure of lambda for CT
# will give message about dropping species. This is okay. There are species in our tree that do not have CT values and the function is removing those
FN_lambda2 <-phylosig(tree = jetztree_diet, FN_vect, method="lambda", test=T) # inputs are trimmed phylogenetic tree and vector of body mass
FN_lambda2 # lambda = 0.912 , p << 0.001
plot.phylosig(FN_lambda2) # plot the likelihood surface for lambda


     # % Diet Plant/Seed

# look at distribution
hist(Dietspp$Diet.PS) # not skewed. definitely not normal, but no log transformation needed

# get vector of C.T values with species names as rownames
PS_vect <- setNames(Dietspp$Diet.PS, rownames(Dietspp))

# drop species with NA values
PS_vect[!is.na(PS_vect)]

# get measure of lambda for CT
# will give message about dropping species. This is okay. There are species in our tree that do not have CT values and the function is removing those
PS_lambda2 <-phylosig(tree = jetztree_diet, PS_vect, method="lambda", test=T) # inputs are trimmed phylogenetic tree and vector of body mass
PS_lambda2 # lambda = 0.828 , p << 0.001
plot.phylosig(PS_lambda2) # plot the likelihood surface for lambda



######################## Phylogenetic Signal (lambda) for NESTING TRAITS ######################## 

#import Coastal_Species_Sensory.rds 

Nest_Traits <- readRDS(here("Outputs", "Coastal_Species_Nest.rds"))
head(Nest_Traits)
colnames(Nest_Traits)

#should I be going off of... Species_Jetz, Species_BirdLife, or Species_eBird???? prob Jetz because I join it to Jetz tree.... 

#do some reformatting and re-naming... 
Nest_Traits <- Nest_Traits %>%
  rename(Species = Species_Jetz)
colnames(Nest_Traits)
#nice 

#now, reformat Species column so that it is in Aaaa_aaaa format! 
Nest_Traits <- Nest_Traits %>%
  mutate(Species = str_replace(Species, " ", "_"))
View(Nest_Traits)

### we need to trim Jetz tree that we imported above to get a tree that only contains coastal species

# create data frame with traits and species as the rownames
Nestspp <- Nest_Traits %>%
  column_to_rownames(var="Species")

### load phylogenetic tree
jetztree <- read.tree(here("Data", "Jetz_ConsensusPhy.tre")) 
jetztree$tip.label # look at species name formatting for tree tips to confirm it matches formatting we are using - yes!

# check that the tree is ultrametric (do all the tree tips line up?)
is.ultrametric(jetztree)
#true! 


# check that all the species in our data are also present in the tree
check_nest <- name.check(jetztree, Nestspp) # there should be many more species in the tree than our data
check_nest

# trim the tree to match the data by dropping all the extra species identified in the previous step
jetztree_nest <-drop.tip(jetztree, check_nest$tree_not_data)

# check again whether tree has same species as data. Should say "OK"
name.check(jetztree_nest, Nestspp)
#OK 


     # Nest safety (may be considered ordinal...)

# look at distribution
hist(Nestspp$nest.safety) # not skewed. definitely not normal, but no log transformation needed

# get vector of nest safety values with species names as rownames
safety_vect <- setNames(Nestspp$nest.safety, rownames(Nestspp))

# drop species with NA values
safety_vect[!is.na(safety_vect)]

# get measure of lambda for nest safety
# will give message about dropping species. This is okay. There are species in our tree that do not have CT values and the function is removing those
safety_lambda2 <-phylosig(tree = jetztree_nest, safety_vect, method="lambda", test=T) # inputs are trimmed phylogenetic tree and vector of body mass
safety_lambda2 # lambda = 0.876 , p << 0.001
plot.phylosig(safety_lambda2) # plot the likelihood surface for lambda



######################## Phylogenetic Signal (lambda) for LIFE HISTORY TRAITS ########################

#import Coastal_Species_Sensory.rds 

LifeHist_Traits <- readRDS(here("Outputs", "Coastal_Species_LifeHistory.rds"))
head(LifeHist_Traits)
colnames(LifeHist_Traits)

#should I be going off of... Species_Jetz, Species_BirdLife, or Species_eBird???? prob Jetz because I join it to Jetz tree.... 

#do some reformatting and re-naming... 
LifeHist_Traits <- LifeHist_Traits %>%
  rename(Species = Species_Jetz)
colnames(LifeHist_Traits)
#nice 

#now, reformat Species column so that it is in Aaaa_aaaa format! 
LifeHist_Traits <- LifeHist_Traits %>%
  mutate(Species = str_replace(Species, " ", "_"))
View(LifeHist_Traits)

### we need to trim Jetz tree that we imported above to get a tree that only contains coastal species

# create data frame with traits and species as the rownames
LifeHistspp <- LifeHist_Traits %>%
  column_to_rownames(var="Species")

### load phylogenetic tree
jetztree <- read.tree(here("Data", "Jetz_ConsensusPhy.tre")) 
jetztree$tip.label # look at species name formatting for tree tips to confirm it matches formatting we are using - yes!

# check that the tree is ultrametric (do all the tree tips line up?)
is.ultrametric(jetztree)
#true! 


# check that all the species in our data are also present in the tree
check_LifeHist <- name.check(jetztree, LifeHistspp) # there should be many more species in the tree than our data
check_LifeHist

# trim the tree to match the data by dropping all the extra species identified in the previous step
jetztree_LifeHist <-drop.tip(jetztree, check_LifeHist$tree_not_data)

# check again whether tree has same species as data. Should say "OK"
name.check(jetztree_LifeHist, LifeHistspp)
#OK 


     # Brood Value 

# look at distribution
hist(LifeHistspp$brood_value) # not skewed, looks normal dist. 

# get vector of nest safety values with species names as rownames
bv_vect <- setNames(LifeHistspp$brood_value, rownames(LifeHistspp))

# drop species with NA values
bv_vect[!is.na(bv_vect)]

# get measure of lambda for nest safety
# will give message about dropping species. This is okay. There are species in our tree that do not have CT values and the function is removing those
bv_lambda2 <-phylosig(tree = jetztree_LifeHist, bv_vect, method="lambda", test=T) # inputs are trimmed phylogenetic tree and vector of body mass
bv_lambda2 # lambda = 0.876 , p << 0.001
plot.phylosig(bv_lambda2) # plot the likelihood surface for lambda



     # Clutch size 

# look at distribution
hist(LifeHistspp$clutch_size) # not normal, but doesn't seem to need transformation... saving it as output for Sarah to look at with me... 

# get vector of nest safety values with species names as rownames
clutch_vect <- setNames(LifeHistspp$clutch_size, rownames(LifeHistspp))

# drop species with NA values
clutch_vect[!is.na(clutch_vect)]

# get measure of lambda for nest safety
# will give message about dropping species. This is okay. There are species in our tree that do not have CT values and the function is removing those
clutch_lambda2 <-phylosig(tree = jetztree_LifeHist, clutch_vect, method="lambda", test=T) # inputs are trimmed phylogenetic tree and vector of body mass
clutch_lambda2 # lambda = 0.876 , p << 0.001
plot.phylosig(clutch_lambda2) # plot the likelihood surface for lambda


     # Longevity 

# look at distribution
hist(LifeHistspp$longevity) #looks pretty normal, no need to transform 

# get vector of nest safety values with species names as rownames
longev_vect <- setNames(LifeHistspp$longevity, rownames(LifeHistspp))

# drop species with NA values
longev_vect[!is.na(longev_vect)]

# get measure of lambda for nest safety
# will give message about dropping species. This is okay. There are species in our tree that do not have CT values and the function is removing those
longev_lambda2 <-phylosig(tree = jetztree_LifeHist, longev_vect, method="lambda", test=T) # inputs are trimmed phylogenetic tree and vector of body mass
longev_lambda2 # lambda = 0.876 , p << 0.001
plot.phylosig(longev_lambda2) # plot the likelihood surface for lambda


     # Developmental mode = NOT continuous 

######################## Phylogenetic Signal (lambda) for SEXUAL SELECTION TRAITS ########################

#import Coastal_SSelect_Sensory.rds 

SS_Traits <- readRDS(here("Outputs", "Coastal_Species_SSelect.rds"))
head(SS_Traits)
colnames(SS_Traits)

#should I be going off of... Species_Jetz, Species_BirdLife, or Species_eBird???? prob Jetz because I join it to Jetz tree.... 

#do some reformatting and re-naming... 
SS_Traits <- SS_Traits %>%
  rename(Species = Species_Jetz)
colnames(SS_Traits)
#nice 

#now, reformat Species column so that it is in Aaaa_aaaa format! 
SS_Traits <- SS_Traits %>%
  mutate(Species = str_replace(Species, " ", "_"))
View(SS_Traits)

### we need to trim Jetz tree that we imported above to get a tree that only contains coastal species

# create data frame with traits and species as the rownames
SSelectspp <- SS_Traits %>%
  column_to_rownames(var="Species")

### load phylogenetic tree
jetztree <- read.tree(here("Data", "Jetz_ConsensusPhy.tre")) 
jetztree$tip.label # look at species name formatting for tree tips to confirm it matches formatting we are using - yes!

# check that the tree is ultrametric (do all the tree tips line up?)
is.ultrametric(jetztree)
#true! 


# check that all the species in our data are also present in the tree
check_SSelect <- name.check(jetztree, SSelectspp) # there should be many more species in the tree than our data
check_SSelect

# trim the tree to match the data by dropping all the extra species identified in the previous step
jetztree_SSelect <-drop.tip(jetztree, check_SSelect$tree_not_data)

# check again whether tree has same species as data. Should say "OK"
name.check(jetztree_SSelect, SSelectspp)
colnames(SSelectspp)
#OK 


     # Brightness Dimorphism 

# look at distribution
hist(SSelectspp$Dichrom_bright) #looks pretty normal, no need to transform 

# get vector of nest safety values with species names as rownames
brightness_vect <- setNames(SSelectspp$Dichrom_bright, rownames(SSelectspp))

# drop species with NA values
brightness_vect[!is.na(brightness_vect)]

# get measure of lambda for nest safety
# will give message about dropping species. This is okay. There are species in our tree that do not have CT values and the function is removing those
brightness_lambda2 <-phylosig(tree = jetztree_SSelect, brightness_vect, method="lambda", test=T) # inputs are trimmed phylogenetic tree and vector of body mass
brightness_lambda2 # lambda = ~0 , p = 1
plot.phylosig(brightness_lambda2) # plot the likelihood surface for lambda


     # Hue Dimorphism 

# look at distribution
hist(SSelectspp$Dichrom_hue) #looks pretty normal, no need to transform 

# get vector of nest safety values with species names as rownames
hue_vect <- setNames(SSelectspp$Dichrom_hue, rownames(SSelectspp))

# drop species with NA values
hue_vect[!is.na(hue_vect)]

# get measure of lambda for nest safety
# will give message about dropping species. This is okay. There are species in our tree that do not have CT values and the function is removing those
hue_lambda2 <-phylosig(tree = jetztree_SSelect, hue_vect, method="lambda", test=T) # inputs are trimmed phylogenetic tree and vector of body mass
hue_lambda2 # lambda = 0.0412 , p = 0.105
plot.phylosig(hue_lambda2) # plot the likelihood surface for lambda


     # Sex selection intensity MALES 

# look at distribution
hist(SSelectspp$sex.sel.m) #looks pretty normal, no need to transform 

# get vector of nest safety values with species names as rownames
SSmale_vect <- setNames(SSelectspp$sex.sel.m, rownames(SSelectspp))

# drop species with NA values
SSmale_vect[!is.na(SSmale_vect)]

# get measure of lambda for nest safety
# will give message about dropping species. This is okay. There are species in our tree that do not have CT values and the function is removing those
SSmale_lambda2 <-phylosig(tree = jetztree_SSelect, SSmale_vect, method="lambda", test=T) # inputs are trimmed phylogenetic tree and vector of body mass
SSmale_lambda2 # lambda = 0.746 , p << 0.001
plot.phylosig(SSmale_lambda2) # plot the likelihood surface for lambda


     # Sex selection intensity FEMALES


# look at distribution
hist(SSelectspp$sex.sel.f) #not normal looking, but no need to transform (most scores are 0)

# get vector of nest safety values with species names as rownames
SSfemale_vect <- setNames(SSelectspp$sex.sel.f, rownames(SSelectspp))

# drop species with NA values
SSfemale_vect[!is.na(SSfemale_vect)]

# get measure of lambda for nest safety
# will give message about dropping species. This is okay. There are species in our tree that do not have CT values and the function is removing those
SSfemale_lambda2 <-phylosig(tree = jetztree_SSelect, SSfemale_vect, method="lambda", test=T) # inputs are trimmed phylogenetic tree and vector of body mass
SSfemale_lambda2 # lambda = 0.938 , p << 0.001
plot.phylosig(SSfemale_lambda2) # plot the likelihood surface for lambda


######################## Phylogenetic Signal (lambda) for SOCIAL TRAITS ######################## 

     # no continuous traits 

######################## Phylogenetic Signal (lambda) for BODY MASS ######################## 

     # Log Body Mass 

#import Coastal_Species_w_Mass.rds 

Mass_Traits <- readRDS(here("Outputs", "Coastal_Species_w_Mass.rds"))
head(Mass_Traits)
colnames(Mass_Traits)

#should I be going off of... Species_Jetz, Species_BirdLife, or Species_eBird???? prob Jetz because I join it to Jetz tree.... 

#do some reformatting and re-naming... 
Mass_Traits <- Mass_Traits %>%
  rename(Species = Species_Jetz)
colnames(Mass_Traits)
#nice 

#now, reformat Species column so that it is in Aaaa_aaaa format! 
Mass_Traits <- Mass_Traits %>%
  mutate(Species = str_replace(Species, " ", "_"))
head(Mass_Traits)

### we need to trim Jetz tree that we imported above to get a tree that only contains coastal species

# create data frame with traits and species as the rownames
Mass.spp <- Mass_Traits %>%
  column_to_rownames(var="Species")

### load phylogenetic tree
jetztree <- read.tree(here("Data", "Jetz_ConsensusPhy.tre")) 
jetztree$tip.label # look at species name formatting for tree tips to confirm it matches formatting we are using - yes!

# check that the tree is ultrametric (do all the tree tips line up?)
is.ultrametric(jetztree)
#true! 


# check that all the species in our data are also present in the tree
check_Mass <- name.check(jetztree, Mass.spp) # there should be many more species in the tree than our data
check_Mass

# trim the tree to match the data by dropping all the extra species identified in the previous step
jetztree_Mass <-drop.tip(jetztree, check_Mass$tree_not_data)

# check again whether tree has same species as data. Should say "OK"
name.check(jetztree_Mass, Mass.spp)
colnames(Mass.spp)
#OK 



#(Log) Body Mass 

# look at distribution
hist(Mass.spp$Mass) # need to log transform 
hist(Mass.spp$Mass_log) # looks good 

# get vector of nest safety values with species names as rownames
mass_vect <- setNames(Mass.spp$Mass_log, rownames(Mass.spp))

# drop species with NA values
mass_vect[!is.na(mass_vect)]

# get measure of lambda for nest safety
# will give message about dropping species. This is okay. There are species in our tree that do not have CT values and the function is removing those
mass_lambda2 <-phylosig(tree = jetztree_Mass, mass_vect, method="lambda", test=T) # inputs are trimmed phylogenetic tree and vector of body mass
mass_lambda2 # lambda = 0.977 , p = 0
plot.phylosig(mass_lambda2) # plot the likelihood surface for lambda


