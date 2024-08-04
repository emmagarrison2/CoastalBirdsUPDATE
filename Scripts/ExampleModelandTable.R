library(nlme)
library(tidyverse)
library(here)
library(broom)
if(!require(broom.mixed)){
  install.packages("broom.mixed")
  require(broom.mixed)
}
library(broom.mixed)
if(!require(flextable)){
  install.packages("flextable")
  require(flextable)
}
library(flextable)
library(phytools)
library(ape)
library(geiger)
library(ggeffects)
if(!require(easystats)){
  install.packages("easystats")
  require(easystats)
}
library(easystats)

# read in trait and species data
Coastal_Sensory <- readRDS(here("Outputs", "Coastal_Species_Sensory.rds"))
View(Coastal_Sensory)
# confirm there are 807 species and 807 unique Jetz species names
nrow(Coastal_Sensory)
length(unique(Coastal_Sensory$Species_Jetz))

# add underscore to Jetz species names and move them to rownames
sensory_jetz <- Coastal_Sensory %>%
  mutate(phy_names = str_replace(Species_Jetz, " ", "_")) %>%
  column_to_rownames(., var = "phy_names")

# import tree
tree <- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))
tree$tip.label # look at species name formatting for tree tips to confirm it matches formatting we are using - yes!
# check that the tree is ultrametric (do all the tree tips line up?)
is.ultrametric(tree)


## EMMA - Making an example model. Please re-format this to match your method/approach
# look at dim light vision for UAI as an example model
# drop all rows/species that are missing C.T values
CT_UAI <- sensory_jetz %>% filter(!is.na(C.T)) %>% filter(!is.na(aveUAI))

# join tree with data
CT_UAI_tree  <- treedata(tree, CT_UAI, sort=T)
CT_UAI_phy <- CT_UAI_tree$phy # rename pruned phylogeny

# rename sorted data and change columns to numeric
CT_UAI_dat <- data.frame(CT_UAI_tree$data) %>% 
  mutate_at(c("aveUAI", "Mass_log", "C.T"), as.numeric) # make columns for model into numeric

# look at distribution of variables we will use in the model
hist(CT_UAI_dat$aveUAI) # UAI scores
hist(CT_UAI_dat$Mass_log) # log transformed body mass
hist(CT_UAI_dat$C.T) # C.T ratio

# run a phylogenetic gls model
CT_UAI_mod <- gls(aveUAI ~ C.T + Mass_log, data = CT_UAI_dat, 
    correlation = corPagel(0, phy=CT_UAI_phy, fixed=T), method = "ML") 

summary(CT_UAI_mod)

# make model summary into a tidy output
# use broom.mixed package to do this
CT_UAI_mod_tidy <- broom.mixed::tidy(CT_UAI_mod) %>%  
  mutate_if(is.numeric, round, 4) # round all the columns that are numeric to have 4 digits after the decimal
CT_UAI_mod_tidy

# at this stage you could write the tidy data frame from above as a csv file 
# write.csv(CT_UAI_mod_tidy, here("Results", "CT_UAI_model.csv")) 
# I would add a folder called Results to the R Project as a place to save these

# we can also use the flextable package on the tidy model summary to make a table that can be pasted into a document
# set defaults for these tables (put this code once at the beginning of the script. this will apply to every table you make)
set_flextable_defaults(
  font.family = "Arial",
    font.size = 12,
    digits = 4)

# use the flextable() function on the tidy model summary
# a table will show up in the Viewer window (on the right). This can be exported
flextable(CT_UAI_mod_tidy)

