##### the objective of this script is to find the phylogenetic distance (PD) and species richness (SR) 
# of each urban tolerance index (UAI, MUTI, and UN)

#install/load appropriate packages 

library(here)
library(tidyverse)
library(ape)
library(geiger)
if(!require(picante)){
  install.packages("picante")
  require(picante)
}
library(picante)
library(phytools)

citation("picante")
packageVersion("picante")



######################## Phylogenetic distance ########################
# in this section, we will use the pd function to get Phylogenetic distance (PD) and
# species richness (SR) for each Urban Tolerance Index.
# Additionally, we will count the number of distinct families in each urban tolerance index


# import coastal species list (these are joined with scores from all 3 indexes)
# trait data not needed, as PD is across the entire index 

Coastal <- read.csv(here("Data", "Coastal_Birds_List.csv"), header=T)
head(Coastal)
colnames(Coastal)
nrow(Coastal) #807

#Species_Jetz has a sci name for all 807 coastal species! 
###use this column to do PD measurements 
#but first, let's rename as "Species" 

Coastal <- Coastal %>%
  rename(Species = Species_Jetz)
colnames(Coastal)
#nice 

#now, reformat Species column so that it is in Aaaa_aaaa format! 
Coastal <- Coastal %>%
  mutate(Species = str_replace(Species, " ", "_"))
head(Coastal)

# get coastal UAI species
CoastalUAI <- Coastal %>% filter(!is.na(aveUAI))
nrow(CoastalUAI) # 798 species

# get coastal UN species
CoastalUN <- Coastal %>% filter(!is.na(Urban))
nrow(CoastalUN) # 129 species

# get coastal MUTI species
CoastalMUTI <- Coastal %>% filter(!is.na(MUTIscore))
nrow(CoastalMUTI) # 130 species

# import Jetz phylogeny

### load phylogenetic tree
jetztree <- read.tree(here("Data", "Jetz_ConsensusPhy.tre")) 
jetztree$tip.label # look at species name formatting for tree tips to confirm it matches formatting we are using - yes!

# check that the tree is ultrametric (do all the tree tips line up?)
is.ultrametric(jetztree)
#True! 

######
# reorganize the data
# the picante package which allows us to calculate phylogenetic distance (PD)
# they have the data organized a specific way and the following steps put our data in the same format
ForPD <- Coastal %>% dplyr::select(Species, aveUAI, Urban, MUTIscore) %>%
  rename(UAI = aveUAI, UN = Urban, MUTI = MUTIscore) %>%
  mutate(
    UAI = if_else(!is.na(UAI), 1, 0), # convert species with index to 1 and species without index to 0
    UN  = if_else(!is.na(UN), 1, 0),
    MUTI = if_else(!is.na(MUTI), 1, 0)) %>% 
  pivot_longer(!Species, names_to="UrbanIndex", values_to="Score") %>%
  pivot_wider(., names_from = Species, values_from = Score) %>%
  column_to_rownames(., var="UrbanIndex")

View(ForPD)



# prune the jetz phylogeny
coastal_jetz <- prune.sample(ForPD, jetztree)
coastal_jetz #807 tips 

# make sure the list of species in ForPD and the species in the pruned tree are in the same order
ForPD <- ForPD[, coastal_jetz$tip.label]


# get Faith's phylogenetic distance for each group of species (one for UAI, UTI, and UN)
PD <- pd(ForPD, coastal_jetz) # pd function from picante package
PD
# The pd function returns two values for each community, the PD and the species richness (SR)
# Faith’s PD is defined as the total branch length spanned by the tree for the species in the group

#interesting that UN is less phylogenetically diverse than MUTI, although they have almost = # of sp 

# trying out other phylogenetic species diversity metrics
if(!require(phyr)){
  install.packages("phyr")
  require(phyr)
}
library(phyr) # measures described in Helmus et al. 2007 Am Nat
# link to paper: https://www.journals.uchicago.edu/doi/10.1086/511334

# We can obtain psv and psr for our data (the authors also present an evenness measure that we can't calculate).
# These measures are explained in the paper as follows:

# The three metrics we present here are derived statistically by considering the 
# value of some unspecified neutral trait shared by all species in a community. 
# As this neutral trait evolves up a phylogenetic tree, 
# speciation occurs, and from this point forward, evolution proceeds independently along each phylogenetic lineage.
# Our metric of phylogenetic species variability (PSV) quantifies how phylogenetic relatedness 
# decreases the variance of this hypothetical unselected trait shared by all species in the community. 
# To calculate PSV, only information about the phylogenetic relatedness of species in a community is needed, 
# not information about any particular trait. 
# Nonetheless, framing this measure in the context of a hypothetical neutral trait gives a metric 
# that has not only an intuitive interpretation but also appealing statistical properties. 
# The second metric quantifies phylogenetic species richness (PSR) as the number of species in a community 
# multiplied by the community’s PSV. 
# This metric is directly comparable to the traditional metric of species richness but includes phylogenetic relatedness. 

# PSV = phylogenetic species variability
# Bound between zero and one
# approaches zero as relatedness of the species increases 
# therefore, higher values therefore reflect greater diversity
psv(ForPD, coastal_jetz, compute.var=F)
# we would report the PSVs value for each urban tolerance index
#UAI - 0.8044951 
#UN - 0.7606829 
#MUTI - 0.7970667 

# PSR = phylogenetic species richness
# Can take on any value, but is directly comparable to actual species richness
# if PSR is much lower than the actual richness, it would indicate a community where species are closely related (i.e. congenerics)
psr(ForPD, coastal_jetz, compute.var=F)
# PSR column gives the phylogenetic species richness and SR gives the actual species richness
#UAI - 641.98712 
#UN - 98.12809 
#MUTI - 103.61867 

#### how many bird families are there for each urban index?

# organize the data again to enable this
colnames(Coastal)

#make sure I know what family column to use 
Coastal_fam_test <- Coastal %>% filter (!is.na(Family_Sci))
nrow(Coastal_fam_test)

families <- Coastal %>% dplyr::select(Species, Family_Sci, aveUAI, Urban, MUTIscore) %>%
  rename(UAI = aveUAI, UN = Urban, MUTI = MUTIscore) %>%
  mutate(
    UAI = if_else(!is.na(UAI), 1, 0), # convert species with index to 1 and species without index to 0
    MUTI  = if_else(!is.na(MUTI), 1, 0),
    UN = if_else(!is.na(UN), 1, 0)) 

# now find the families associated with each urban tolerance index

#UAI
familiesUAI <- families %>% filter(UAI==1) %>% 
  distinct(Family_Sci) # this will print a list of families in UAI
nrow(familiesUAI) # gives number of families in UAI = 81

#UN
familiesUN <- families %>% filter(UN==1) %>% 
  distinct(Family_Sci) # this will print a list of families in UN
nrow(familiesUN) # gives number of families in UN = 24

#MUTI
familiesMUTI <- families %>% filter(MUTI==1) %>% 
  distinct(Family_Sci) # this will print a list of families in MUTI
nrow(familiesMUTI) # gives number of families in MUTI = 33


# combine the family counts with PD and Species Richness (SR)
PD$Family <- as.vector(c(nrow(familiesUAI), nrow(familiesUN), nrow(familiesMUTI)))
PD # we will want the info in this data frame in a table in the manuscript

#let's write this object (PD) as a csv, so that if things change in the future, we can 
#simply re-run this Step5 script with updated starting point. 
write.csv(PD, here("Notes", "PD_SR_Fam.csv"))

#done
