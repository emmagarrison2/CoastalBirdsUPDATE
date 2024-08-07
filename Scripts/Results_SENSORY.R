######Results - Sensory Data Traits 

#load required packages 

library(nlme)
library(tidyverse)
library(here)
library(broom)
library(broom.mixed)
library(flextable)
library(phytools)
library(ape)
library(geiger)
library(ggeffects)
library(easystats)
if(!require(phylolm)){
  install.packages("phylolm")
  require(phylolm)
}
library(phylolm)

###################### Prep


# read in Sensory trait and species data

Coastal_Sensory <- readRDS(here("Outputs", "Coastal_Species_Sensory.rds"))

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

#For Urban (UN), we need to change Urban (U) = 1, and nonurban (N) = 0.

sensory_jetz$Urban <- ifelse(sensory_jetz$Urban == "U", 1, 0)
View(sensory_jetz)

################################### C.T. ################################### 




########## C.T. and UAI ###########

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
confint(CT_UAI_mod)

# make model summary into a tidy output
# use broom.mixed package to do this
CT_UAI_mod_tidy <- broom.mixed::tidy(CT_UAI_mod) %>%  
  mutate_if(is.numeric, round, 4) # round all the columns that are numeric to have 4 digits after the decimal
CT_UAI_mod_tidy


# we will combine this CT_UAI_mod_tidy tidy output with the other results... and this will become our sensory traits results table! 


########## C.T. and MUTI ###########
colnames(sensory_jetz)

# look at dim light vision for UAI as an example model
# drop all rows/species that are missing C.T values
CT_MUTI <- sensory_jetz %>% filter(!is.na(C.T)) %>% filter(!is.na(MUTIscore))

# join tree with data
CT_MUTI_tree  <- treedata(tree, CT_MUTI, sort=T)
CT_MUTI_phy <- CT_MUTI_tree$phy # rename pruned phylogeny

# rename sorted data (TREE) and change columns to numeric (DAT)
CT_MUTI_dat <- data.frame(CT_MUTI_tree$data) %>% 
  mutate_at(c("MUTIscore", "Mass_log", "C.T"), as.numeric) # make columns for model into numeric

# look at distribution of variables we will use in the model
hist(CT_MUTI_dat$MUTIscore) # MUTI scores - skewed right
hist(CT_MUTI_dat$Mass_log) # log transformed body mass - looks good 
hist(CT_MUTI_dat$C.T) # C.T ratio - looks good 

# run a phylogenetic gls model
CT_MUTI_mod <- gls(MUTIscore ~ C.T + Mass_log, data = CT_MUTI_dat, 
                  correlation = corPagel(0, phy=CT_MUTI_phy, fixed=T), method = "ML") 

summary(CT_MUTI_mod)
confint(CT_MUTI_mod)

# make model summary into a tidy output
# use broom.mixed package to do this
CT_MUTI_mod_tidy <- broom.mixed::tidy(CT_MUTI_mod) %>%  
  mutate_if(is.numeric, round, 3) # round all the columns that are numeric to have 4 digits after the decimal
CT_MUTI_mod_tidy


########## C.T. and UN ###########

colnames(sensory_jetz)

# look at dim light vision for UAI as an example model
# drop all rows/species that are missing C.T values
CT_UN <- sensory_jetz %>% filter(!is.na(C.T)) %>% filter(!is.na(Urban))

# join tree with data
CT_UN_tree  <- treedata(tree, CT_UN, sort=T)
CT_UN_phy <- CT_UN_tree$phy # rename pruned phylogeny

# rename sorted data (TREE) and change columns to numeric (DAT)
CT_UN_dat <- data.frame(CT_UN_tree$data) %>% 
  mutate_at(c("Urban", "Mass_log", "C.T"), as.numeric) # make columns for model into numeric

# look at distribution of variables we will use in the model
hist(CT_UN_dat$Urban) # MUTI scores - slightly skewed left
hist(CT_UN_dat$Mass_log) # log transformed body mass - looks okay 
hist(CT_UN_dat$C.T) # C.T ratio - looks okay 

# run a phylogenetic gls model
#CT_MUTI_mod <- gls(MUTIscore ~ C.T + Mass_log, data = CT_UN_dat, 
                  # correlation = corPagel(0, phy=CT_UN_phy, fixed=T), method = "ML") 

CT_UN_mod <- phylolm(Urban~ C.T + Mass_log, data= CT_UN_dat, phy= CT_UN_phy, model="lambda") 


summary(CT_UN_mod)
confint(CT_UN_mod)

#is this worth all the work when I can just hand-copy it in? 

# make model summary into a tidy outputs

# use broom.mixed package to do this

#have to make custom tidy function for phylolm 
#tidy_phylolm <- function(model) {
 # summary_model <- summary(model)
#  coefficients <- as.data.frame(summary_model$coefficients)
 # rownames_to_column(coefficients, var = "term") %>%
  #  rename(estimate = Estimate,
   #        std.error = `Std. Error`,
    #       statistic = `t value`,
     #      p.value = `Pr(>|t|)`)
#}

#CT_UN_mod_tidy <- tidy_phylolm(CT_UN_mod) %>%  
 # mutate_if(is.numeric, round, 4) # round all the columns that are numeric to have 4 digits after the decimal
#CT_UN_mod_tidy




################################### Peak Frequency ################################### 



########## Peak Frequency and UAI ###########
colnames(sensory_jetz)

# look at dim light vision for UAI as an example model
# drop all rows/species that are missing peak_freq values
pf_UAI <- sensory_jetz %>% filter(!is.na(peak_freq)) %>% filter(!is.na(aveUAI))

# join tree with data
pf_UAI_tree  <- treedata(tree, pf_UAI, sort=T)
pf_UAI_phy <- pf_UAI_tree$phy # rename pruned phylogeny

# rename sorted data and change columns to numeric
pf_UAI_dat <- data.frame(pf_UAI_tree$data) %>% 
  mutate_at(c("aveUAI", "Mass_log", "peak_freq"), as.numeric) # make columns for model into numeric

# look at distribution of variables we will use in the model
hist(pf_UAI_dat$aveUAI) # UAI scores - looks fine 
hist(pf_UAI_dat$Mass_log) # log transformed body mass
hist(pf_UAI_dat$peak_freq) # peak_freq ratio - bit skewed to the right, but looks fine 

# run a phylogenetic gls model
pf_UAI_mod <- gls(aveUAI ~ peak_freq + Mass_log, data = pf_UAI_dat, 
                  correlation = corPagel(0, phy=pf_UAI_phy, fixed=T), method = "ML") 

summary(pf_UAI_mod)
confint(pf_UAI_mod)

# make model summary into a tidy output
# use broom.mixed package to do this
pf_UAI_mod_tidy <- broom.mixed::tidy(pf_UAI_mod) %>%  
  mutate_if(is.numeric, round, 3) # round all the columns that are numeric to have 4 digits after the decimal
pf_UAI_mod_tidy




########## Peak Frequency and MUTI ###########

colnames(sensory_jetz)

# look at dim light vision for UAI as an example model
# drop all rows/species that are missing peak_freq values
pf_MUTI <- sensory_jetz %>% filter(!is.na(peak_freq)) %>% filter(!is.na(MUTIscore))

# join tree with data
pf_MUTI_tree  <- treedata(tree, pf_MUTI, sort=T)
pf_MUTI_phy <- pf_MUTI_tree$phy # rename pruned phylogeny

# rename sorted data and change columns to numeric
pf_MUTI_dat <- data.frame(pf_MUTI_tree$data) %>% 
  mutate_at(c("MUTIscore", "Mass_log", "peak_freq"), as.numeric) # make columns for model into numeric

# look at distribution of variables we will use in the model
hist(pf_MUTI_dat$MUTIscore) # MUTIscore scores - skewed right
hist(pf_MUTI_dat$Mass_log) # log transformed body mass - not necessarily normal but looks ok 
hist(pf_MUTI_dat$peak_freq) # peak_freq ratio - skewed right



# run a phylogenetic gls model
pf_MUTI_mod <- gls(MUTIscore ~ peak_freq + Mass_log, data = pf_MUTI_dat, 
                  correlation = corPagel(0, phy=pf_MUTI_phy, fixed=T), method = "ML") 

summary(pf_MUTI_mod)
confint(pf_MUTI_mod)

# make model summary into a tidy output
# use broom.mixed package to do this
pf_MUTI_mod_tidy <- broom.mixed::tidy(pf_MUTI_mod) %>%  
  mutate_if(is.numeric, round, 3) # round all the columns that are numeric to have 4 digits after the decimal
pf_MUTI_mod_tidy




########## Peak Frequency and UN ###########


colnames(sensory_jetz)

# look at dim light vision for UAI as an example model
# drop all rows/species that are missing C.T values
pf_UN <- sensory_jetz %>% filter(!is.na(peak_freq)) %>% filter(!is.na(Urban))

# join tree with data
pf_UN_tree  <- treedata(tree, pf_UN, sort=T)
pf_UN_phy <- pf_UN_tree$phy # rename pruned phylogeny

# rename sorted data (TREE) and change columns to numeric (DAT)
pf_UN_dat <- data.frame(pf_UN_tree$data) %>% 
  mutate_at(c("Urban", "Mass_log", "peak_freq"), as.numeric) # make columns for model into numeric

# look at distribution of variables we will use in the model
hist(pf_UN_dat$Urban) # MUTI scores - slightly skewed left
hist(pf_UN_dat$Mass_log) # log transformed body mass - looks okay 
hist(pf_UN_dat$peak_freq) # peak_freq ratio - looks okay 

# run a phylolm model

pf_UN_mod <- phylolm(Urban~ peak_freq + Mass_log, data= pf_UN_dat, phy= pf_UN_phy, model="lambda") 


summary(pf_UN_mod)
confint(pf_UN_mod)


############################################## Time to Combine ################


#combine TEST 

combined_tidy_test <- bind_rows(CT_UAI_mod_tidy, pf_UAI_mod_tidy)
combined_tidy_test

# looks good! will have to combine in order at end (i.e. CT_UAI, CT_MUTI, CT_UN, PF_UAI, PF_MUTI, PF_UN)


#combine CT_UAI, CT_MUTI, (no UN), PF_UAI, and PF_MUTI (no UN)

combined_Sensory_results <- bind_rows(CT_UAI_mod_tidy, CT_MUTI_mod_tidy, pf_UAI_mod_tidy, pf_MUTI_mod_tidy)
combined_Sensory_results
