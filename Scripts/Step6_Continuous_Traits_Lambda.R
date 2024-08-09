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

# get vector of C.T values with species names as rownames
pf_vect <- setNames(Sensoryspp$peak_freq, rownames(Sensoryspp))

# drop species with NA values
pf_vect[!is.na(pf_vect)]

# get measure of lambda for CT
# will give message about dropping species. This is okay. There are species in our tree that do not have CT values and the function is removing those
pf_lambda2 <-phylosig(tree = jetztree_sensory, pf_vect, method="lambda", test=T) # inputs are trimmed phylogenetic tree and vector of body mass
pf_lambda2 # lambda = 0.936 , p << 0.001
plot.phylosig(pf_lambda2) # plot the likelihood surface for lambda


######################## Phylogenetic Signal (lambda) for DIET TRAITS ######################## 

######################## Phylogenetic Signal (lambda) for NESTING TRAITS ######################## 

######################## Phylogenetic Signal (lambda) for LIFE HISTORY TRAITS ########################

######################## Phylogenetic Signal (lambda) for SEXUAL SELECTION TRAITS ########################

######################## Phylogenetic Signal (lambda) for SOCIAL TRAITS ######################## 

######################## Phylogenetic Signal (lambda) for BODY MASS ######################## 

#import Coastal_Species_w_Mass.rds