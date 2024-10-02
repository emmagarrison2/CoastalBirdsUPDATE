##########################################################
######### Coastal Birds Urban Tolerance Project ##########
##########################################################
# Goal: Run Phylogenetic Logistic Trait Models for UN
# Authors: Sarah L. Jennings, Emma M. Garrison
##########################################################

# The objective of this script is to run UN models with phyloglm(), 
# this is a logistic regression, which is more fitting for binary response variables like UN 


# load required packages 
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
library(phylolm)
library(logistf)


# we will be boot strapping the final model for each trait
# therefore, we need 21 random numbers for set.seed() to produce consistent results 
set.seed(456)
seednums <- sample(100:999, 21, replace=F) # get 21 integers with 3 digits (between 100 and 999)
seednums


#######################################################
###################### BODY MASS ######################
#######################################################


###################### Prep

#load in "Coastal_Species_w_Mass.rds" - since this contains coastal species and all diet trait variables :) 

C_mass_dat <- readRDS(here("Outputs", "Coastal_Species_w_Mass.rds"))
str(C_mass_dat)

C_mass_dat2 <- C_mass_dat %>%
  mutate(Species_Jetz  = str_replace(Species_Jetz, " ", "_"))
str(C_mass_dat2)

C_mass_dat2$Urban <- ifelse(C_mass_dat2$Urban == "U", 1, 0)
View(C_mass_dat2)
colnames(C_mass_dat2)


######################## UN and body mass ##########################

# lets first simplify a NEW database by removing records where we don't have an UN / brood_value
UNDataUT <- C_mass_dat2 %>% filter(!is.na(Urban)) 
MassData3 <- UNDataUT %>% filter(!is.na(Mass_log)) 
length(MassData3$Mass_log)
#129 species with UN and Mass_log

colnames(MassData3)

###### add and pair tree

# add rownames to data
row.names(MassData3) <- MassData3$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Massphydat3 <- treedata(tree_out,MassData3, sort=T)

Massphy3 <- Massphydat3$phy
MassTraitDat3 <- as.data.frame(Massphydat3$data)

str(MassTraitDat3)
length(MassTraitDat3$Mass_log)
#129


### convert traits of interest to numeric

MassTraitDat3$Urban <- as.numeric(MassTraitDat3$Urban)
MassTraitDat3$Mass_log <- as.numeric(MassTraitDat3$Mass_log)



#lets run the model using Phyloglm

############################# 

# default method ="logistic_MPLE"

phyglm_UN_Mass <- phyloglm( Urban ~ Mass_log, 
                               data = MassTraitDat3, 
                               phy = Massphy3, 
                               boot = 1000)
summary(phyglm_UN_Mass)
# this fails to converge
# also produces warning about alpha reaching the upper bound

# scale and center mass and run model 
phyglm_UN_Mass_scale <- phyloglm( Urban ~ scale(Mass_log), 
                            data = MassTraitDat3, 
                            phy = Massphy3, 
                            boot = 1000) 

summary(phyglm_UN_Mass_scale)
# this also fails converge

# print AIC values for models with different upper bounds
# try intervals of 0.1 from 0 up to 4 for log.alpha.bound
# 4 is the default setting for log.alpha.bound
for (i in seq(0, 4, by = 0.1)) {
  print(phyloglm(Urban ~ scale(Mass_log), 
                 data = MassTraitDat3, 
                 phy = Massphy3,
                 log.alpha.bound = i)$aic)
}

# there is higher support for larger alpha values (later in the printed list) that correspond to low phylogenetic signal

# we can also look at logLik, which should show a similar pattern to AIC
for (i in seq(0, 4, by = 0.1)) {
  print(phyloglm(Urban ~ scale(Mass_log), 
                 data = MassTraitDat3, 
                 phy = Massphy3,
                 log.alpha.bound = i)$logLik)
}
# same as above. Lower log likelihood values for models with higher values of alpha/low phylo signal


# try model with scale(Mass_log) by fixing alpha at upper bounds
# upper bounds for alpha in our models is:
t <- phyglm_UN_Mass$mean.tip.height # set the mean tip height of the tree as a variable called t
t # t for our tree is 97.56
# use t to find the upper bound of alpha for our model
exp(4)/t # equals 0.5596322
exp(-4)/t # this is the lower bound
# Note: alpha can take on values from 0 to infinity, but around log.alpha.bound = 4 the model should produce results similar to a non-phylogenetic logistic model

# log.alpha.bound = 4 sets the upper bound at this value (0.5596332)
# we can also set start.alpha, the alpha value where the model begins its search
# we can't set start.alpha to exactly 0.5596332 nad have the run model (so it won't let us actually fix the value of alpha)
# we can give a start.alpha value very close to the upper bound (e.g., 0.557)
# this work around basically constrains the model's search area for the optimal alpha value within a very small range of values (similar to fixing it)
set.seed(seednums[1])
phyglm_UN_Mass_fix <- phyloglm(Urban ~ scale(Mass_log), 
                                  data = MassTraitDat3, 
                                  phy = Massphy3, 
                                  start.alpha = 0.557,
                                  log.alpha.bound = 4, 
                                  boot=1000)
summary(phyglm_UN_Mass_fix)


# save model
saveRDS(phyglm_UN_Mass_fix, here("Outputs", "phyglm_UN_Mass_fix.rds"))
# load model
phyglm_UN_Mass_fix <- readRDS(here("Outputs", "phyglm_UN_Mass_fix.rds"))

# as alpha is at upper bounds, we can also look at a regular logistic model
# use logistf package for this
# this runs as logistic regression with Firth's correction (a bias reduction method)
# Ives and Garland 2010 recommend log.alpha.bound = 4 as the limits (and this is default setting in phyloglm)
# they specify 4 because when the model reaches the upper bounds of this limit,
# the model estimates should be very similar to a logistic model using Firth's correction
glm_UN_Mass <- logistf(Urban ~ scale(Mass_log), data = MassTraitDat3)
summary(glm_UN_Mass)
summary(phyglm_UN_Mass_fix) # compare estimates
# very similar coefficients

# get values for results table for the fixed alpha model with scale(Mass_log)
summary(phyglm_UN_Mass_fix)
# use the bootstrapped CIs in the summary

# get alpha, t, and half life for the model
(phyglm_UN_Mass_fix$mean.tip.height) # this is t (mean tip height) of the tree
(alpha_mass <- phyglm_UN_Mass_fix$alpha) # this is alpha
(hl_mass <- log(2)/alpha_mass) # this is the half-life for the model
# compare the value for half life with the mean tip height of the tree
# compared to t, the half life is small -> therefore we conclude there is low phylogenetic signal


#######################################################
###################### SENSORY TRAITS #################
#######################################################


###################### Prep


#load in "CoastalUT_dat30April24.rds" - since this contains coastal species and all trait variables :) 

C_Sensory_dat <- readRDS(here("Outputs", "Coastal_Species_Sensory.rds"))
str(C_Sensory_dat)

C_Sensory_dat2 <- C_Sensory_dat %>%
  mutate(Species_Jetz  = str_replace(Species_Jetz, " ", "_"))
str(C_Sensory_dat2)

C_Sensory_dat2$Urban <- ifelse(C_Sensory_dat2$Urban == "U", 1, 0)
View(C_Sensory_dat2)
colnames(C_Sensory_dat2)


############################## C.T. and UN ####################################

# lets first simplify a NEW database by removing records where we don't have an UAI value from Neate-Clegg
URBANSensory <- C_Sensory_dat2 %>% filter(!is.na(Urban)) 
SensoryTraitData2 <- URBANSensory %>% filter(!is.na(C.T)) 
length(SensoryTraitData2$C.T)
#38 species with UAI and CT

###### add and pair tree

# add rownames to data
row.names(SensoryTraitData2) <- SensoryTraitData2$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Sensoryphydat2 <- treedata(tree_out,SensoryTraitData2, sort=T)

Sensoryphy2 <- Sensoryphydat2$phy
SensoryTraitDat2 <- as.data.frame(Sensoryphydat2$data)

str(SensoryTraitDat2)
length(SensoryTraitDat2$C.T)
#38

### convert traits of interest to numeric

SensoryTraitDat2$Urban <- as.numeric(SensoryTraitDat2$Urban)
SensoryTraitDat2$Mass_log <- as.numeric(SensoryTraitDat2$Mass_log)
SensoryTraitDat2$C.T <- as.numeric(SensoryTraitDat2$C.T)

# let's run a model using Phyloglm 


############################ Run Models

phyglm_UN_CT <- phyloglm( Urban ~ C.T + Mass_log, 
                              data = SensoryTraitDat2, 
                              phy = Sensoryphy2, 
                              boot = 1000) 
# this model fails to converge
summary(phyglm_UN_CT)

# scale and center the explanatory variables and run the model
set.seed(seednums[2])
phyglm_UN_CT_scale <- phyloglm(Urban ~ scale(C.T) + scale(Mass_log), 
                          data = SensoryTraitDat2, 
                          phy = Sensoryphy2, 
                          boot = 1000) 

# this version does converge
summary(phyglm_UN_CT_scale)

# save model
saveRDS(phyglm_UN_CT_scale, here("Outputs", "phyglm_UN_CT_scale.rds"))
# load model
phyglm_UN_CT_scale <- readRDS(here("Outputs", "phyglm_UN_CT_scale.rds"))

# get alpha, t, and half life for the model
(alpha_CT <- phyglm_UN_CT_scale$alpha) # alpha
(phyglm_UN_CT_scale$mean.tip.height) # t (aka mean tip height)
(hl_CT<- log(2)/alpha_CT) # half-life
# compared to t, this is a small Half-Life -> low phylogenetic signal

############################# Peak Frequency and UN #################################

# lets first simplify a NEW database by removing records where we don't have an UAI value from Neate-Clegg
UNSensory <- C_Sensory_dat2 %>% filter(!is.na(Urban)) 
SensoryTraitData5 <- UNSensory %>% filter(!is.na(peak_freq)) 
length(SensoryTraitData5$peak_freq)
#129 species with UAI and peak frequency

###### add and pair tree

# add rownames to data
row.names(SensoryTraitData5) <- SensoryTraitData5$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Sensoryphydat5 <- treedata(tree_out,SensoryTraitData5, sort=T)

Sensoryphy5 <- Sensoryphydat5$phy
SensoryTraitDat5 <- as.data.frame(Sensoryphydat5$data)

str(SensoryTraitDat5)
length(SensoryTraitDat5$peak_freq)
#129

### convert traits of interest to numeric

SensoryTraitDat5$Urban <- as.numeric(SensoryTraitDat5$Urban)
SensoryTraitDat5$Mass_log <- as.numeric(SensoryTraitDat5$Mass_log)
SensoryTraitDat5$peak_freq <- as.numeric(SensoryTraitDat5$peak_freq)

# let's run a model using Phyloglm 


##################### Run Models

phyglm_UN_pf <- phyloglm( Urban ~ peak_freq + Mass_log, 
                          data = SensoryTraitDat5, 
                          phy = Sensoryphy5, 
                          boot = 1000) 

summary(phyglm_UN_pf)

# run model with scaled and centered explanatory variables
set.seed(seednums[3])
phyglm_UN_pf_scale <- phyloglm( Urban ~ scale(peak_freq) + scale(Mass_log), 
                          data = SensoryTraitDat5, 
                          phy = Sensoryphy5, 
                          boot = 1000)
# warning: alpha reached upper bound
# but model does converge
summary(phyglm_UN_pf_scale)

# save model
saveRDS(phyglm_UN_pf_scale, here("Outputs", "phyglm_UN_pf_scale.rds"))
# load model
phyglm_UN_pf_scale <- readRDS(here("Outputs", "phyglm_UN_pf_scale.rds"))

# as alpha is at upper bound, also look at regular logistic model
glm_UN_pf_scale <- logistf(Urban ~ scale(peak_freq) + scale(Mass_log), 
                           data = SensoryTraitDat5)
summary(glm_UN_pf_scale)
# coefficient for peak freq is quite different, but we reach the same conclusions

# get alpha, t, and half life for the model
(phyglm_UN_pf_scale$mean.tip.height) # t
(alpha_pfreq <- phyglm_UN_pf_scale$alpha) # alpha
(hl_pfreq<- log(2)/alpha_pfreq) # half life
# low half life relative to mean tip height -> low phylogenetic signal


#######################################################
###################### DIET TRAITS #####################
#######################################################


###################### Prep


#load in "Coastal_Species_Diet.rds" - since this contains coastal species and all diet trait variables :) 

C_Diet_dat <- readRDS(here("Outputs", "Coastal_Species_Diet.rds"))
str(C_Diet_dat)

C_Diet_dat2 <- C_Diet_dat %>%
  mutate(Species_Jetz  = str_replace(Species_Jetz, " ", "_"))
str(C_Diet_dat2)

C_Diet_dat2$Urban <- ifelse(C_Diet_dat2$Urban == "U", 1, 0)
View(C_Diet_dat2)
colnames(C_Diet_dat2)

######################## UN and % Diet Invertebrates ##########################


# lets first simplify a NEW database by removing records where we don't have an UAI / %dietinv
UNDataUT <- C_Diet_dat2 %>% filter(!is.na(Urban)) 
DietData3 <- UNDataUT %>% filter(!is.na(Diet.Inv)) 
length(DietData3$Diet.Inv)
#129 species with UAI and CT

###### add and pair tree

DietData3 <- as.data.frame(DietData3)
# add rownames to data
row.names(DietData3) <- DietData3$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Dietphydat3 <- treedata(tree_out,DietData3, sort=T)

Dietphy3 <- Dietphydat3$phy
DietTraitDat3 <- as.data.frame(Dietphydat3$data)

str(DietTraitDat3)
length(DietTraitDat3$Diet.Inv)
#129

### convert traits of interest to numeric

DietTraitDat3$Urban <- as.numeric(DietTraitDat3$Urban)
DietTraitDat3$Mass_log <- as.numeric(DietTraitDat3$Mass_log)
DietTraitDat3$Diet.Inv <- as.numeric(DietTraitDat3$Diet.Inv)

# let's run a model using Phyloglm 


###################### Run Models

phyglm_UN_Invert <- phyloglm( Urban ~ Diet.Inv + Mass_log, 
                            data = DietTraitDat3, 
                            phy = Dietphy3, 
                            boot = 1000) 

#Warning message:
  #In phyloglm(Urban ~ Diet.Inv + Mass_log, data = DietTraitDat3, phy = Dietphy3,  :
                #the boundary of the linear predictor has been reached during the optimization procedure.
             # You can increase this bound by increasing 'btol'.

# run model with scaled/center explanatory variables
phyglm_UN_Invert_scale <- phyloglm( Urban ~ scale(Diet.Inv) + scale(Mass_log), 
                              data = DietTraitDat3, 
                              phy = Dietphy3, 
                              boot = 1000)

summary(phyglm_UN_Invert_scale) 
# this fails to converge
# we also get a warning that alpha reached upper bounds


# print AIC values for models with different upper bounds
# intervals of 0.1 from 0 up to 4 for log.alpha.bound
# I tried saving these to a vector but I am only able to get it to save AIC for models where it converged
# so a lot of values are lost and it's hard to know which AIC corresponds to which value of log.alpha.bound 
# this is not super sophisticated but it works for what we need it for
for (i in seq(0, 4, by = 0.1)) {
  print(phyloglm(Urban ~ scale(Diet.Inv) + scale(Mass_log), 
               data = DietTraitDat3, 
               phy = Dietphy3,
               log.alpha.bound = i)$aic)
}

# higher support for larger alpha values (later in the printed list) that correspond to low phylogenetic signal

# we can also look at logLik, which should show a similar pattern to AIC
for (i in seq(0, 4, by = 0.1)) {
  print(phyloglm(Urban ~ scale(Diet.Inv) + scale(Mass_log), 
                 data = DietTraitDat3, 
                 phy = Dietphy3,
                 log.alpha.bound = i)$logLik)
}
# same as above. Lower log likelihood values for models with higher values of alpha/lower phylo signal


# try to fix alpha at exp(4)/t (at the default upper bounds)
# this fails to converge
phyglm_UN_Invert_fix_4 <- phyloglm( Urban ~ scale(Diet.Inv) + scale(Mass_log), 
                                      data = DietTraitDat3, 
                                      phy = Dietphy3, 
                                      start.alpha = 0.55,
                                      log.alpha.bound = 4)

# try increasing upper bounds a small amount over 4
# both these models converge
set.seed(seednums[4])
phyglm_UN_Invert_fix_4.05 <- phyloglm( Urban ~ scale(Diet.Inv) + scale(Mass_log), 
                                    data = DietTraitDat3, 
                                    phy = Dietphy3, 
                                    start.alpha = 0.55,
                                    log.alpha.bound = 4.05, boot=1000)
summary(phyglm_UN_Invert_fix_4.05) 
phyglm_UN_Invert_fix_4.05$aic

phyglm_UN_Invert_fix_4.1 <- phyloglm( Urban ~ scale(Diet.Inv) + scale(Mass_log), 
                                       data = DietTraitDat3, 
                                       phy = Dietphy3, 
                                       start.alpha = 0.55,
                                       log.alpha.bound = 4.1, boot=1000)
summary(phyglm_UN_Invert_fix_4.1) 
phyglm_UN_Invert_fix_4.1$aic


# coefficients fairly stable across 4.05 and 4.1 models
# save model
saveRDS(phyglm_UN_Invert_fix_4.05, here("Outputs", "phyglm_UN_Invert_fix.rds"))
# load model
phyglm_UN_Invert_fix_4.05 <- readRDS(here("Outputs", "phyglm_UN_Invert_fix.rds"))

# also look at non-phylogenetic logistic model
glm_UN_Invert <- logistf(Urban ~ scale(Diet.Inv) + scale(Mass_log), 
        data = DietTraitDat3)
summary(glm_UN_Invert)
# we reach similar conclusions as two models above

# get alpha, t, and half life for the model
# using phyglm_UN_Invert_fix_4.05 as final model
(phyglm_UN_Invert_fix_4.05$mean.tip.height) # t
(alpha_invert <- phyglm_UN_Invert_fix_4.05$alpha) # alpha
(hl_invert <- log(2)/alpha_invert) # half life
# small half life -> low phylogenetic signal

######################## UN and % Diet Vertebrates ##########################


# lets first simplify a NEW database by removing records where we don't have an UAI / %dietinv
UNDataUT <- C_Diet_dat2 %>% filter(!is.na(Urban)) 
DietData6 <- UNDataUT %>% filter(!is.na(Diet.Vert)) 
length(DietData6$Diet.Vert)
#129 species with UAI and CT

###### add and pair tree

DietData6 <- as.data.frame(DietData6)
# add rownames to data
row.names(DietData6) <- DietData6$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Dietphydat6 <- treedata(tree_out,DietData6, sort=T)

Dietphy6 <- Dietphydat6$phy
DietTraitDat6 <- as.data.frame(Dietphydat6$data)

str(DietTraitDat6)
length(DietTraitDat6$Diet.Vert)
#129

### convert traits of interest to numeric

DietTraitDat6$Urban <- as.numeric(DietTraitDat6$Urban)
DietTraitDat6$Mass_log <- as.numeric(DietTraitDat6$Mass_log)
DietTraitDat6$Diet.Vert <- as.numeric(DietTraitDat6$Diet.Vert)


#lets run the model using Phyloglm!  


####################### Run Models


phyglm_UN_Vert <- phyloglm( Urban ~ Diet.Vert + Mass_log, 
                              data = DietTraitDat6, 
                              phy = Dietphy6, 
                              boot = 1000) 
summary(phyglm_UN_Vert)

# using scaled explanatory variables
phyglm_UN_Vert_scale <- phyloglm( Urban ~ scale(Diet.Vert) + scale(Mass_log), 
                            data = DietTraitDat6, 
                            phy = Dietphy6, 
                            boot = 1000)
summary(phyglm_UN_Vert_scale)
# this fails to converge

# print AIC values for models with different upper bounds
# intervals of 0.1 from 0 up to 4
for (i in seq(0, 4, by = 0.1)) {
  print(phyloglm(Urban ~ scale(Diet.Vert) + scale(Mass_log), 
                 data = DietTraitDat6, 
                 phy = Dietphy6,
                 log.alpha.bound = i)$aic)
}
# AIC values support models with larger values of alpha (low phylo signal)

for (i in seq(0, 4, by = 0.1)) {
  print(phyloglm(Urban ~ scale(Diet.Vert) + scale(Mass_log), 
                 data = DietTraitDat6, 
                 phy = Dietphy6,
                 log.alpha.bound = i)$logLik)
}
# we reach the same conclusion with log likelihood

# try to fix alpha at upper bounds
set.seed(seednums[5])
phyglm_UN_Vert_fix <- phyloglm( Urban ~ scale(Diet.Vert) + scale(Mass_log), 
                                     data = DietTraitDat6, 
                                     phy = Dietphy6, 
                                     start.alpha = 0.557,
                                     log.alpha.bound = 4, boot = 1000)
summary(phyglm_UN_Vert_fix)

# save model
saveRDS(phyglm_UN_Vert_fix, here("Outputs", "phyglm_UN_Vert_fix.rds"))
# load model
phyglm_UN_Vert_fix <- readRDS(here("Outputs", "phyglm_UN_Vert_fix.rds"))

# as alpha is at upper bound, look at a non-phylogenetic model
glm_UN_Vert <- logistf(Urban ~ scale(Diet.Vert) + scale(Mass_log), 
                       data = DietTraitDat6)
summary(glm_UN_Vert)
# very similar

# get alpha, t, and half life for the model
(phyglm_UN_Vert_fix$mean.tip.height) # t
(alpha_vert <- phyglm_UN_Vert_fix$alpha) # alpha
(hl_vert <- log(2)/alpha_vert) # half life
# compared to t, this is a small half life

######################## UN and % Diet Plant/Seed ##########################


# lets first simplify a NEW database by removing records where we don't have an UAI / % diet plant seed 
UNDataUT <- C_Diet_dat2 %>% filter(!is.na(Urban)) 
DietData9 <- UNDataUT %>% filter(!is.na(Diet.PS)) 
length(DietData9$Diet.PS)
#129 species with UAI and CT

###### add and pair tree

DietData9 <- as.data.frame(DietData9)
# add rownames to data
row.names(DietData9) <- DietData9$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Dietphydat9 <- treedata(tree_out,DietData9, sort=T)

Dietphy9 <- Dietphydat9$phy
DietTraitDat9 <- as.data.frame(Dietphydat9$data)

str(DietTraitDat9)
length(DietTraitDat9$Diet.PS)
#129

### convert traits of interest to numeric

DietTraitDat9$Urban <- as.numeric(DietTraitDat9$Urban)
DietTraitDat9$Mass_log <- as.numeric(DietTraitDat9$Mass_log)
DietTraitDat9$Diet.PS <- as.numeric(DietTraitDat9$Diet.PS)


#lets run the model using Phyloglm!  


######################## # Run Models


phyglm_UN_PS <- phyloglm( Urban ~ Diet.PS + Mass_log, 
                            data = DietTraitDat9, 
                            phy = Dietphy9, 
                            boot = 1000) 

summary(phyglm_UN_PS)

# try model with scaled and centered explanatory variables
set.seed(seednums[6])
phyglm_UN_PS_scale <- phyloglm( Urban ~ scale(Diet.PS) + scale(Mass_log), 
                          data = DietTraitDat9, 
                          phy = Dietphy9, 
                          boot = 1000) 

# this converges
summary(phyglm_UN_PS_scale)

# save model
saveRDS(phyglm_UN_PS_scale, here("Outputs", "phyglm_UN_PS_scale.rds"))
# load model
phyglm_UN_PS_scale <- readRDS(here("Outputs", "phyglm_UN_PS_scale.rds"))

# get alpha, t, and half life for the model
(phyglm_UN_PS_scale$mean.tip.height) # t
(alpha_PS <- phyglm_UN_PS_scale$alpha) # alpha
(hl_PS <- log(2)/alpha_PS) # half life
# small compared to t -> low phylogenetic signal

######################## UN and % Diet Fruit/Nectar ##########################


# lets first simplify a NEW database by removing records where we don't have an UAI / % diet plant seed 
UNDataUT <- C_Diet_dat2 %>% filter(!is.na(Urban)) 
DietData12 <- UNDataUT %>% filter(!is.na(Diet.FN)) 
length(DietData12$Diet.FN)
#129 species with UAI and CT

###### add and pair tree

DietData12 <- as.data.frame(DietData12)
# add rownames to data
row.names(DietData12) <- DietData12$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Dietphydat12 <- treedata(tree_out,DietData12, sort=T)

Dietphy12 <- Dietphydat12$phy
DietTraitDat12 <- as.data.frame(Dietphydat12$data)

str(DietTraitDat12)
length(DietTraitDat12$Diet.FN)
#129

### convert traits of interest to numeric

DietTraitDat12$Urban <- as.numeric(DietTraitDat12$Urban)
DietTraitDat12$Mass_log <- as.numeric(DietTraitDat12$Mass_log)
DietTraitDat12$Diet.FN <- as.numeric(DietTraitDat12$Diet.FN)




#lets run the model using Phyloglm!  


####################### Run Models


phyglm_UN_FN <- phyloglm( Urban ~ Diet.FN + Mass_log, 
                          data = DietTraitDat12, 
                          phy = Dietphy12, 
                          boot = 1000) 

# this model fails to converge
# multiple other warning messages


# run model with scaled and centered explanatory variables
set.seed(seednums[7])
phyglm_UN_FN_scale <- phyloglm( Urban ~ scale(Diet.FN) + scale(Mass_log), 
                          data = DietTraitDat12, 
                          phy = Dietphy12, 
                          boot = 1000)

# this model converges
summary(phyglm_UN_FN_scale)


# save model
saveRDS(phyglm_UN_FN_scale, here("Outputs", "phyglm_UN_FN_scale.rds"))
# load model
phyglm_UN_FN_scale <- readRDS(here("Outputs", "phyglm_UN_FN_scale.rds"))


# get alpha, t, and half life for the model
(phyglm_UN_FN_scale$mean.tip.height) # t
(alpha_FN <- phyglm_UN_FN_scale$alpha) # alpha
(hl_FN <- log(2)/alpha_FN) # half life
# compared to t, this is a small half life


#######################################################
###################### LIFE HISTORY TRAITS ############
#######################################################



###################### Prep


#load in "Coastal_Species_LifeHistory.rds" - since this contains coastal species and all diet trait variables :) 

C_LifeHist_dat <- readRDS(here("Outputs", "Coastal_Species_LifeHistory.rds"))
str(C_LifeHist_dat)

C_LifeHist_dat2 <- C_LifeHist_dat %>%
  mutate(Species_Jetz  = str_replace(Species_Jetz, " ", "_"))
str(C_LifeHist_dat2)

C_LifeHist_dat2$Urban <- ifelse(C_LifeHist_dat2$Urban == "U", 1, 0)
View(C_LifeHist_dat2)
colnames(C_LifeHist_dat2)



######################## UN and % brood value ##########################


# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UNDataUT <- C_LifeHist_dat2 %>% filter(!is.na(Urban)) 
LifehistData3 <- UNDataUT %>% filter(!is.na(brood_value)) 
length(LifehistData3$brood_value)
#102 species with UAI and brood_value

###### add and pair tree

# add rownames to data
LifehistData3 <- as.data.frame(LifehistData3)
row.names(LifehistData3) <- LifehistData3$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Lifehistphydat3 <- treedata(tree_out,LifehistData3, sort=T)

Lifehistphy3 <- Lifehistphydat3$phy
LifehistTraitDat3 <- as.data.frame(Lifehistphydat3$data)

str(LifehistTraitDat3)
length(LifehistTraitDat3$brood_value)
#102

### convert traits of interest to numeric

LifehistTraitDat3$Urban <- as.numeric(LifehistTraitDat3$Urban)
LifehistTraitDat3$Mass_log <- as.numeric(LifehistTraitDat3$Mass_log)
LifehistTraitDat3$brood_value <- as.numeric(LifehistTraitDat3$brood_value)




#lets run the model using Phyloglm!  

############################## Run Models

phyglm_UN_bv <- phyloglm( Urban ~ brood_value + Mass_log, 
                          data = LifehistTraitDat3, 
                          phy = Lifehistphy3, 
                          boot = 1000) 
# this model fails to converge
summary(phyglm_UN_bv)


# run the model with center and scaled variables
set.seed(seednums[8])
phyglm_UN_bv_scale <- phyloglm( Urban ~ scale(brood_value) + scale(Mass_log), 
                          data = LifehistTraitDat3, 
                          phy = Lifehistphy3, 
                          boot = 1000) 
summary(phyglm_UN_bv_scale) 
# this successfully converges 
# alpha is at upper bounds


# save model
saveRDS(phyglm_UN_bv_scale, here("Outputs", "phyglm_UN_bv_scale.rds"))
# load model
phyglm_UN_bv_scale <- readRDS(here("Outputs", "phyglm_UN_bv_scale.rds"))


# run a non-phylogenetic logistic model for comparison
glm_UN_bv <- logistf(Urban ~ scale(brood_value) + scale(Mass_log), 
  data = LifehistTraitDat3)            
summary(glm_UN_bv)
# we reach the same conclusions
# some differences in coefficients though


# get alpha, t, and half life for the model
(phyglm_UN_bv_scale$mean.tip.height) # t
(alpha_bv <- phyglm_UN_bv_scale$alpha) # alpha
(hl_bv <- log(2)/alpha_bv) # half-life
# small half life relative to t -> low phylogenetic signal


######################## UN and % clutch size ##########################


# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UNDataUT <- C_LifeHist_dat2 %>% filter(!is.na(Urban)) 
LifehistData6 <- UNDataUT %>% filter(!is.na(clutch_size)) 
length(LifehistData6$clutch_size)
#122 species with UN and clutch_size

###### add and pair tree

# add rownames to data
LifehistData6 <- as.data.frame(LifehistData6)
row.names(LifehistData6) <- LifehistData6$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Lifehistphydat6 <- treedata(tree_out,LifehistData6, sort=T)

Lifehistphy6 <- Lifehistphydat6$phy
LifehistTraitDat6 <- as.data.frame(Lifehistphydat6$data)

str(LifehistTraitDat6)
length(LifehistTraitDat6$clutch_size)
#122

### convert traits of interest to numeric

LifehistTraitDat6$Urban <- as.numeric(LifehistTraitDat6$Urban)
LifehistTraitDat6$Mass_log <- as.numeric(LifehistTraitDat6$Mass_log)
LifehistTraitDat6$clutch_size <- as.numeric(LifehistTraitDat6$clutch_size)



#lets run the model using Phyloglm!  

################################ Run Models

phyglm_UN_clutch <- phyloglm( Urban ~ clutch_size + Mass_log, 
                          data = LifehistTraitDat6, 
                          phy = Lifehistphy6, 
                          boot = 1000) 
# this model fails to converge


# try model with scaled and centered variables
phyglm_UN_clutch_scale <- phyloglm( Urban ~ scale(clutch_size) + scale(Mass_log), 
                              data = LifehistTraitDat6, 
                              phy = Lifehistphy6, 
                              boot = 1000) 

summary(phyglm_UN_clutch_scale) 
# still fails to converge


# print AIC values for models with different upper bounds
# intervals of 0.1 from 0 up to 4
for (i in seq(0, 4, by = 0.1)) {
  print(phyloglm(Urban ~ scale(clutch_size) + scale(Mass_log), 
                 data = LifehistTraitDat6, 
                 phy = Lifehistphy6,
                 log.alpha.bound = i)$aic)
}
# AIC values support models with larger values of alpha (low phylo signal)



# fix alpha at upper bounds
set.seed(seednums[9])
phyglm_UN_clutch_fix <- phyloglm( Urban ~ scale(clutch_size) + scale(Mass_log), 
                              data = LifehistTraitDat6, 
                              phy = Lifehistphy6, 
                              log.alpha.bound = 4,
                              start.alpha = 0.557,
                              boot = 1000) 

summary(phyglm_UN_clutch_fix)


# save model
saveRDS(phyglm_UN_clutch_fix, here("Outputs", "phyglm_UN_clutch_fix.rds"))
# load model
phyglm_UN_clutch_fix <- readRDS(here("Outputs", "phyglm_UN_clutch_fix.rds"))



# run non-phylogenetic glm for comparison
library(logistf)
glm_UN_clutch <- logistf(Urban ~ scale(clutch_size) + scale(Mass_log), 
                     data = LifehistTraitDat6)            
summary(glm_UN_clutch)
# we reach the same conclusions
# fairly similar coefficients to fixed model above



# get alpha, t, and half life for the model
(phyglm_UN_clutch_fix$mean.tip.height) # t
(alpha_clutch <- phyglm_UN_clutch_fix$alpha) # alpha
(hl_clutch <- log(2)/alpha_clutch) # half life
# small half life relative to t -> low phylogenetic signal

######################## UN and % longevity ##########################


# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UNDataUT <- C_LifeHist_dat2 %>% filter(!is.na(Urban)) 
LifehistData9 <- UNDataUT %>% filter(!is.na(longevity)) 
length(LifehistData9$longevity)
#129 species with UAI and longevity

###### add and pair tree

# add rownames to data
LifehistData9 <- as.data.frame(LifehistData9)
row.names(LifehistData9) <- LifehistData9$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Lifehistphydat9 <- treedata(tree_out,LifehistData9, sort=T)

Lifehistphy9 <- Lifehistphydat9$phy
LifehistTraitDat9 <- as.data.frame(Lifehistphydat9$data)

str(LifehistTraitDat9)
length(LifehistTraitDat9$longevity)
#129

### convert traits of interest to numeric

LifehistTraitDat9$Urban <- as.numeric(LifehistTraitDat9$Urban)
LifehistTraitDat9$Mass_log <- as.numeric(LifehistTraitDat9$Mass_log)
LifehistTraitDat9$longevity <- as.numeric(LifehistTraitDat9$longevity)



#lets run the model using Phyloglm!  

########################## Run Models

phyglm_UN_long <- phyloglm( Urban ~ longevity + Mass_log, 
                              data = LifehistTraitDat9, 
                              phy = Lifehistphy9, 
                              boot = 1000) 
summary(phyglm_UN_long)


# run model using scaled and centered variables
phyglm_UN_long_scale <- phyloglm( Urban ~ scale(longevity) + scale(Mass_log), 
                            data = LifehistTraitDat9, 
                            phy = Lifehistphy9, 
                            boot = 1000) 

# this model fails to converge
summary(phyglm_UN_long_scale)


# print AIC values for models with different upper bounds
# intervals of 0.1 from 0 up to 4
for (i in seq(0, 4, by = 0.1)) {
  print(phyloglm(Urban ~ scale(longevity) + scale(Mass_log), 
                 data = LifehistTraitDat9, 
                 phy = Lifehistphy9,
                 log.alpha.bound = i)$aic)
}
# AIC values support models with larger values of alpha (low phylo signal)


# try fixing alpha at upper bounds
set.seed(seednums[10])
phyglm_UN_long_fix <- phyloglm( Urban ~ scale(longevity) + scale(Mass_log), 
                                  data = LifehistTraitDat9, 
                                  phy = Lifehistphy9, 
                                  log.alpha.bound = 4,
                                  start.alpha = 0.55,
                                  boot = 1000)
summary(phyglm_UN_long_fix) # this converges


# save model
saveRDS(phyglm_UN_long_fix, here("Outputs", "phyglm_UN_long_fix.rds"))
# load model
phyglm_UN_long_fix <- readRDS(here("Outputs", "phyglm_UN_long_fix.rds"))


# because alpha at upper bounds, compare with non-phylogenetic logistic regression
glm_UN_long <- logistf(Urban ~ scale(longevity) + scale(Mass_log), 
                      data = LifehistTraitDat9)
summary(glm_UN_long)
# this is similar to model with log.alpha.bound fixed at 4


# get alpha, t, and half life for the model
(phyglm_UN_long_fix$mean.tip.height) # t
(alpha_long <- phyglm_UN_long_fix$alpha) # alpha
(hl_long <- log(2)/alpha_long) # half life
# compared to t, this is a small half-life


######################## UN and % developmental mode ##########################


# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UNDataUT <- C_LifeHist_dat2 %>% filter(!is.na(Urban)) 
LifehistData12 <- UNDataUT %>% filter(!is.na(developmental_mode)) 
length(LifehistData12$developmental_mode)
#129 species with MUTI and developmental_mode

###### add and pair tree

# add rownames to data
LifehistData12 <- as.data.frame(LifehistData12)
row.names(LifehistData12) <- LifehistData12$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Lifehistphydat12 <- treedata(tree_out,LifehistData12, sort=T)

Lifehistphy12 <- Lifehistphydat12$phy
LifehistTraitDat12 <- as.data.frame(Lifehistphydat12$data)

str(LifehistTraitDat12)
length(LifehistTraitDat12$developmental_mode)
#129

### convert traits of interest to numeric

LifehistTraitDat12$Urban <- as.numeric(LifehistTraitDat12$Urban)
LifehistTraitDat12$Mass_log <- as.numeric(LifehistTraitDat12$Mass_log)
LifehistTraitDat12$developmental_mode <- as.numeric(LifehistTraitDat12$developmental_mode)




#lets run the model using Phyloglm!  

################################## Run Models

phyglm_UN_develop <- phyloglm(Urban ~ developmental_mode + Mass_log, 
                            data = LifehistTraitDat12, 
                            phy = Lifehistphy12, 
                            boot = 1000) 
summary(phyglm_UN_develop) # this converges
# alpha at upper bound


# model with scaled and centered explanatory variables
# not scaling developmental mode as it is 0/1 binary variable
phyglm_UN_develop_scale <- phyloglm( Urban ~ developmental_mode + scale(Mass_log), 
                               data = LifehistTraitDat12, 
                               phy = Lifehistphy12,
                               boot = 1000) 
# this fails to converge
summary(phyglm_UN_develop_scale)


# print AIC values for models with different upper bounds
# intervals of 0.1 from 0 up to 4
for (i in seq(0, 4, by = 0.1)) {
  print(phyloglm(Urban ~ developmental_mode + scale(Mass_log), 
                 data = LifehistTraitDat12, 
                 phy = Lifehistphy12,
                 log.alpha.bound = i)$aic)
}
# AIC values support models with larger values of alpha (low phylo signal)


# try fixing model at upper bounds for alpha
set.seed(seednums[11])
phyglm_UN_develop_fix <- phyloglm( Urban ~ developmental_mode + scale(Mass_log), 
                                     data = LifehistTraitDat12, 
                                     phy = Lifehistphy12,
                                     log.alpha.bound = 4,
                                     start.alpha = 0.557,
                                    boot = 1000) 

# this converges
summary(phyglm_UN_develop_fix)
# 95% bootstrapped CI for developmental mode does not overlap 0, but p-value is 0.12
# I ran the model without the bootstrapping and it gives the exact same p-value
# I suspect p-values are derived from the original data (although I can't find anything in package documentation)
# while the 95% CI is derived from bootstrapping process and thus more reliable
# mark as a significant trend based on our criteria to use 95% CI


# save model
saveRDS(phyglm_UN_develop_fix, here("Outputs", "phyglm_UN_develop_fix.rds"))
# load model
phyglm_UN_develop_fix <- readRDS(here("Outputs", "phyglm_UN_develop_fix.rds"))

summary(phyglm_UN_develop_fix)
summary(phyglm_UN_develop_fix_OG )
summary(phyglm_UN_develop_fix_SARAH )
# compare results with a non-phylogenetic logistic model
glm_UN_develop <- logistf(Urban ~ developmental_mode + scale(Mass_log), 
                          data = LifehistTraitDat12)

summary(glm_UN_develop) # similar coefficients


# get alpha, t, and half life for the model
(phyglm_UN_develop_fix$mean.tip.height) # t
(alpha_dev <- phyglm_UN_develop_fix$alpha) # alpha
(hl_dev <- log(2)/alpha_dev) # half life
#compared to t, this is a small Half-Life


#######################################################
###################### NESTING TRAITS #################
#######################################################



###################### Prep


#load in "Coastal_Species_Nest.rds" - since this contains coastal species and all diet trait variables :) 

C_Nest_dat <- readRDS(here("Outputs", "Coastal_Species_Nest.rds"))
str(C_Nest_dat)
nrow(C_Nest_dat)
#807- correct 

C_Nest_dat2 <- C_Nest_dat %>%
  mutate(Species_Jetz  = str_replace(Species_Jetz, " ", "_"))
str(C_Nest_dat2)
nrow(C_Nest_dat)
#807 - correct

C_Nest_dat2$Urban <- ifelse(C_Nest_dat2$Urban == "U", 1, 0)
View(C_Nest_dat2)
colnames(C_Nest_dat2)
nrow(C_Nest_dat2)
#807 - good 



######################## UN and % Nest Strategy ##########################



# do not run --> only 8 species for UN with "enclosed" = 0 Nest Strategy 
# too small of an N for this category 




######################## UN and %  Nest Site LOW  ##########################
# 0 = not low
# 1 = LOW

# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UAIDataUT <- C_Nest_dat2 %>% filter(!is.na(Urban)) 
NestData6 <- UAIDataUT %>% filter(!is.na(NestSite_Low)) 
length(NestData6$NestSite_Low)
#129 species with Urban and NestSite_Low

#there is no Log_mass column in this .rds ... let's put one in! 

# add column for log transformed body mass
NestData6 <- NestData6 %>%
  mutate(Mass_log = log(Mass))

colnames(NestData6)

###### add and pair tree

# add rownames to data
row.names(NestData6) <- NestData6$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Nestphydat6 <- treedata(tree_out,NestData6, sort=T)

Nestphy6 <- Nestphydat6$phy
NestTraitDat6 <- as.data.frame(Nestphydat6$data)

str(NestTraitDat6)
length(NestTraitDat6$NestSite_Low)
#129


### convert traits of interest to numeric

NestTraitDat6$Urban <- as.numeric(NestTraitDat6$Urban)
NestTraitDat6$Mass_log <- as.numeric(NestTraitDat6$Mass_log)
NestTraitDat6$NestSite_Low <- as.numeric(NestTraitDat6$NestSite_Low)


#lets run a Phyloglm model for all UN species with Nest Site Low data 
# (aka, we will not exclude "ambiguous-nesting species" in this model)


############################## Run Models

phyglm_UN_nest_low <- phyloglm( Urban ~ NestSite_Low + Mass_log, 
                               data = NestTraitDat6, 
                               phy = Nestphy6, 
                               boot = 1000) 

summary(phyglm_UN_nest_low)
# converges. alpha at upper bounds



# model with scaled and centered mass
phyglm_UN_nest_low_scale <- phyloglm( Urban ~ NestSite_Low + scale(Mass_log), 
                                data = NestTraitDat6, 
                                phy = Nestphy6, 
                                boot = 1000) 

summary(phyglm_UN_nest_low_scale)
# fails to converge
# alpha at upper bounds


# print AIC values for models with different upper bounds
# intervals of 0.1 from 0 up to 4
for (i in seq(0, 4, by = 0.1)) {
  print(phyloglm(Urban ~  NestSite_Low + scale(Mass_log), 
                 data = NestTraitDat6, 
                 phy = Nestphy6,
                 log.alpha.bound = i)$aic)
}
# AIC values support models with larger values of alpha (low phylo signal)


# try fixing alpha at upper bounds
# I am giving the model a little more searching space for alpha because AIC is actually lowest slightly below log.alpha.bounds = 4
exp(3.7)/(phyglm_UN_nest_low_scale$mean.tip.height) # equals 0.41. Using this as start.alpha

set.seed(seednums[12])
#here
phyglm_UN_nest_low_fix <- phyloglm( Urban ~ NestSite_Low + scale(Mass_log), 
                                      data = NestTraitDat6, 
                                      phy = Nestphy6, 
                                      log.alpha.bound = 4,
                                      start.alpha = 0.41,
                                      boot = 1000)
summary(phyglm_UN_nest_low_fix)
# model converges
# alpha at upper bounds (= 0.559)


# save model
saveRDS(phyglm_UN_nest_low_fix, here("Outputs", "phyglm_UN_nest_low_fix.rds"))
# load model
phyglm_UN_nest_low_fix <- readRDS(here("Outputs", "phyglm_UN_nest_low_fix.rds"))


# compare with non-phylogenetic model
glm_UN_nest_low <- logistf(Urban ~ NestSite_Low + scale(Mass_log), 
                           data = NestTraitDat6)
summary(glm_UN_nest_low)
# some differences for coefficients but same conclusions reached


# get alpha, t, and half life for the model
(phyglm_UN_nest_low_fix$mean.tip.height) # t
(alpha_Nlow <- phyglm_UN_nest_low_fix$alpha) # alpha
(hl_NLow <- log(2)/alpha_Nlow) # half life
# compared to t, this is a small half life

 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#NOW, let's try REMOVING all species that have both HIGH and LOW nesting scores of 1 
#### idea - are the "non-ambiguous" nesting species the ones driving this significant relationship 
#### between nesting and UAI index? 

#filter 
UNDataUT <- C_Nest_dat2 %>% filter(!is.na(Urban)) 
NestData6_only <- UNDataUT %>% filter(!is.na(NestSite_Low)) %>% 
  filter(!(NestSite_Low == 1 & NestSite_High == 1))
length(NestData6_only$NestSite_Low)
#104 species with UAI and ONLY NestSite_Low (not also high nesters)
129-104 # = 25 --> there are 25 species that were both High and Low nesters and had UN scores 

#there is no Log_mass column in this .rds ... let's put one in! 

# add column for log transformed body mass
NestData6_only <- NestData6_only %>%
  mutate(Mass_log = log(Mass))

colnames(NestData6_only)

###### add and pair tree

# add rownames to data
row.names(NestData6_only) <- NestData6_only$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Nestphydat6_only <- treedata(tree_out,NestData6_only, sort=T)

Nestphy6_only <- Nestphydat6_only$phy
NestTraitDat6_only <- as.data.frame(Nestphydat6_only$data)

str(NestTraitDat6_only)
length(NestTraitDat6_only$NestSite_Low)
#104


### convert traits of interest to numeric

NestTraitDat6_only$Urban <- as.numeric(NestTraitDat6_only$Urban)
NestTraitDat6_only$Mass_log <- as.numeric(NestTraitDat6_only$Mass_log)
NestTraitDat6_only$NestSite_Low <- as.numeric(NestTraitDat6_only$NestSite_Low)



#lets run a Phyloglm model 

      #but this time, for species who ONLY nest "low" ! 


####################### Run Model
# model using scaled Mass_log
set.seed(seednums[13])
#here
phyglm_UN_nest_low_only_scale <- phyloglm( Urban ~ NestSite_Low + scale(Mass_log), 
                                data = NestTraitDat6_only, 
                                phy = Nestphy6_only,
                                boot = 1000)


summary(phyglm_UN_nest_low_only_scale)
# alpha at upper bounds
# puzzling... p-value for NestSite_low flagged as significant
# but the bootstrapped CI interval strongly overlaps zero
# this seems off
# I think the CI is more trustworthy
# I ran the model without the bootstrapping and it gives the exact same p-value
# I suspect p-values are calculated on original dataset only... 
# We should use the bootstrapped 95% CI here over the p-value


# save model
saveRDS(phyglm_UN_nest_low_only_scale, here("Outputs", "phyglm_UN_nest_low_only_scale.rds"))
# load model
phyglm_UN_nest_low_only_scale <- readRDS(here("Outputs", "phyglm_UN_nest_low_only_scale.rds"))


# get alpha, t, and half life for the model
(phyglm_UN_nest_low_only_scale$mean.tip.height) # t
(alpha_Nlowonly <- phyglm_UN_nest_low_only_scale$alpha) # alpha
(hl_Nlowonly <- log(2)/alpha_Nlowonly) # half life
# compared to t, this is a small half-life -> low phylogenetic signal



######################## UN and %  Nest Site HIGH ##########################
# 0 = not high
# 1 = HIGH

# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UNDataUT <- C_Nest_dat2 %>% filter(!is.na(Urban)) 
NestData9 <- UNDataUT %>% filter(!is.na(NestSite_High)) 
length(NestData9$NestSite_High)
#129 species with Urban and NestSite_High

#there is no Log_mass column in this .rds ... let's put one in! 

# add column for log transformed body mass
NestData9 <- NestData9 %>%
  mutate(Mass_log = log(Mass))

colnames(NestData9)

###### add and pair tree

# add rownames to data
row.names(NestData9) <- NestData9$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Nestphydat9 <- treedata(tree_out,NestData9, sort=T)

Nestphy9 <- Nestphydat9$phy
NestTraitDat9 <- as.data.frame(Nestphydat9$data)

str(NestTraitDat9)
length(NestTraitDat9$NestSite_High)
#129


### convert traits of interest to numeric

NestTraitDat9$Urban <- as.numeric(NestTraitDat9$Urban)
NestTraitDat9$Mass_log <- as.numeric(NestTraitDat9$Mass_log)
NestTraitDat9$NestSite_High <- as.numeric(NestTraitDat9$NestSite_High)


#lets run a Phyloglm model for all UN species with Nest Site Low data 
# (aka, we will not exclude "ambiguous-nesting species" in this model)


############################### Run Models

phyglm_UN_nest_high <- phyloglm( Urban ~ NestSite_High + Mass_log, 
                                data = NestTraitDat9, 
                                phy = Nestphy9, 
                                boot = 1000) 

summary(phyglm_UN_nest_high) # alpha at upper bounds


# scale and center mass to be consistent with other models
phyglm_UN_nest_high_scale <- phyloglm( Urban ~ NestSite_High + scale(Mass_log), 
                                 data = NestTraitDat9, 
                                 phy = Nestphy9, 
                                 boot = 1000)
summary(phyglm_UN_nest_high_scale) # this fails to converge



# print AIC values for models with different upper bounds
# intervals of 0.1 from 0 up to 4
for (i in seq(0, 4, by = 0.1)) {
  print(phyloglm(Urban ~  NestSite_High + scale(Mass_log), 
                 data = NestTraitDat9, 
                 phy = Nestphy9,
                 log.alpha.bound = i)$aic)
}
# AIC values support models with larger values of alpha (low phylo signal)



# try fixing alpha at upper bounds
# give the model a little more searching space for alpha because AIC is actually lowest slightly below log.alpha.bounds = 4
exp(3.8)/(phyglm_UN_nest_high_scale$mean.tip.height) # equals 0.45. Using this as start.alpha
set.seed(seednums[14])
#here
phyglm_UN_nest_high_fix <- phyloglm( Urban ~ NestSite_High + scale(Mass_log), 
                                    data = NestTraitDat9, 
                                    phy = Nestphy9, 
                                    log.alpha.bound = 4,
                                    start.alpha = 0.45,
                                    boot = 1000)
summary(phyglm_UN_nest_high_fix)
# model converges
# nest site high is important


# save model
saveRDS(phyglm_UN_nest_high_fix, here("Outputs", "phyglm_UN_nest_high_fix.rds"))
# load model
phyglm_UN_nest_high_fix <- readRDS(here("Outputs", "phyglm_UN_nest_high_fix.rds"))


# look at a non-phylogenetic logistic model
glm_UN_nest_high <- logistf(Urban ~ NestSite_High + scale(Mass_log), 
        data = NestTraitDat9)

summary(glm_UN_nest_high) # we reach same conclusions


# get alpha, t, and half life for the model
(phyglm_UN_nest_high_fix$mean.tip.height) # t
(alpha_Nhigh <- phyglm_UN_nest_high_fix$alpha) # alpha
(hl_Nhigh <- log(2)/alpha_Nhigh) # half life
#compared to t, this is a small half life


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#NOW, let's try REMOVING all species that have both HIGH and LOW nesting scores of 1 
#### idea - are the "non-ambiguous" nesting species the ones driving this significant relationship 
#### between nesting and UAI index? 



#filter 
UNDataUT <- C_Nest_dat2 %>% filter(!is.na(Urban)) 
NestData9_only <- UNDataUT %>% filter(!is.na(NestSite_High)) %>% 
  filter(!(NestSite_Low == 1 & NestSite_High == 1))
length(NestData9_only$NestSite_High)
#104 species with UAI and ONLY NestSite_High (not also high nesters)
129-104 # = 25 --> there are 25 species that were both High and Low nesters and had UN scores 

#there is no Log_mass column in this .rds ... let's put one in! 

# add column for log transformed body mass
NestData9_only <- NestData9_only %>%
  mutate(Mass_log = log(Mass))

colnames(NestData9_only)

###### add and pair tree

# add rownames to data
row.names(NestData9_only) <- NestData9_only$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Nestphydat9_only <- treedata(tree_out,NestData9_only, sort=T)

Nestphy9_only <- Nestphydat9_only$phy
NestTraitDat9_only <- as.data.frame(Nestphydat9_only$data)

str(NestTraitDat9_only)
length(NestTraitDat9_only$NestSite_High)
#104


### convert traits of interest to numeric

NestTraitDat9_only$Urban <- as.numeric(NestTraitDat9_only$Urban)
NestTraitDat9_only$Mass_log <- as.numeric(NestTraitDat9_only$Mass_log)
NestTraitDat9_only$NestSite_High <- as.numeric(NestTraitDat9_only$NestSite_High)



#lets run a Phyloglm model 

     #but this time, for species who ONLY nest "high" ! 


########################## Run Model

set.seed(seednums[15])
#here
phyglm_UN_nest_high_only_scale <- phyloglm( Urban ~ NestSite_High + scale(Mass_log), 
                                     data = NestTraitDat9_only, 
                                     phy = Nestphy9_only, 
                                     boot = 1000) 
# model converges
# alpha at upper bounds
summary(phyglm_UN_nest_high_only_scale)
# nest site high is no longer important


# save model
saveRDS(phyglm_UN_nest_high_only_scale, here("Outputs", "phyglm_UN_nest_high_only_scale.rds"))
# load model
phyglm_UN_nest_high_only_scale <- readRDS(here("Outputs", "phyglm_UN_nest_high_only_scale.rds"))


# get alpha, t, and half life for the model
(phyglm_UN_nest_high_only_scale$mean.tip.height) # t
(alpha_Nhighonly <- phyglm_UN_nest_high_only_scale$alpha) # alpha
(hl_Nhighonly <- log(2)/alpha_Nhighonly) # half life
# compared to t, this is a small half life -> low phylogenetic signal





######################## UN and Nest Safety ##########################


# lets first simplify a NEW database by removing records where we don't have an UN / brood_value
UNDataUT <- C_Nest_dat2 %>% filter(!is.na(Urban)) 
NestData12 <- UNDataUT %>% filter(!is.na(nest.safety)) 
length(NestData12$nest.safety)
#129 species with Urban and nest.safety

#there is no Log_mass column in this .rds ... let's put one in! 

# add column for log transformed body mass
NestData12 <- NestData12 %>%
  mutate(Mass_log = log(Mass))

colnames(NestData12)

###### add and pair tree

# add rownames to data
row.names(NestData12) <- NestData12$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Nestphydat12 <- treedata(tree_out,NestData12, sort=T)

Nestphy12 <- Nestphydat12$phy
NestTraitDat12 <- as.data.frame(Nestphydat12$data)

str(NestTraitDat12)
length(NestTraitDat12$nest.safety)
#129


### convert traits of interest to numeric

NestTraitDat12$Urban <- as.numeric(NestTraitDat12$Urban)
NestTraitDat12$Mass_log <- as.numeric(NestTraitDat12$Mass_log)
NestTraitDat12$nest.safety <- as.numeric(NestTraitDat12$nest.safety)



#lets run a Phyloglm model!


######################## Run Models

phyglm_UN_nest_safety <- phyloglm( Urban ~ nest.safety + Mass_log, 
                                 data = NestTraitDat12, 
                                 phy = Nestphy12, 
                                 boot = 1000)
# this model fails to converge


# try model with scaled and centered variables
set.seed(seednums[16])
#here
phyglm_UN_nest_safety_scale <- phyloglm( Urban ~ scale(nest.safety) + scale(Mass_log), 
                                   data = NestTraitDat12, 
                                   phy = Nestphy12, 
                                   boot = 1000) 
summary(phyglm_UN_nest_safety_scale)
# this model converges
# this is a positive trend for nest.safety based on bootstrapped 95% CI (don't use p-value which is less reliable)



# save model
saveRDS(phyglm_UN_nest_safety_scale, here("Outputs", "phyglm_UN_nest_safety_scale.rds"))
# load model
phyglm_UN_nest_safety_scale <- readRDS(here("Outputs", "phyglm_UN_nest_safety_scale.rds"))



# get alpha, t, and half life for the model
(phyglm_UN_nest_safety_scale$mean.tip.height) # t
(alpha_Nsafe <- phyglm_UN_nest_safety_scale$alpha) # alpha
(hl_Nsafe <- log(2)/alpha_Nsafe) # half life
# compared to t, this is small half life


############################################################
################## SEXUAL SELECTION TRAITS #################
###########################################################



###################### Prep

#load in "Coastal_Species_Nest.rds" - since this contains coastal species and all diet trait variables :) 

C_SSelect_dat <- readRDS(here("Outputs", "Coastal_Species_SSelect.rds"))
str(C_SSelect_dat)

C_SSelect_dat2 <- C_SSelect_dat %>%
  mutate(Species_Jetz  = str_replace(Species_Jetz, " ", "_"))
str(C_SSelect_dat2)

C_SSelect_dat2$Urban <- ifelse(C_SSelect_dat2$Urban == "U", 1, 0)
View(C_SSelect_dat2)
colnames(C_SSelect_dat2)




######################## UN and Dichromatism - BRIGHTNESS ##########################

# lets first simplify a NEW database by removing records where we don't have an UN / brood_value
UNDataUT <- C_SSelect_dat2 %>% filter(!is.na(Urban)) 
SSData3 <- UNDataUT %>% filter(!is.na(Dichrom_bright)) 
length(SSData3$Dichrom_bright)
#61 species with UN and Dichrom_bright

colnames(SSData3)

###### add and pair tree

# add rownames to data
row.names(SSData3) <- SSData3$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

SSphydat3 <- treedata(tree_out,SSData3, sort=T)

SSphy3 <- SSphydat3$phy
SSTraitDat3 <- as.data.frame(SSphydat3$data)

str(SSTraitDat3)
length(SSTraitDat3$Dichrom_bright)
#61


### convert traits of interest to numeric

SSTraitDat3$Urban <- as.numeric(SSTraitDat3$Urban)
SSTraitDat3$Mass_log <- as.numeric(SSTraitDat3$Mass_log)
SSTraitDat3$Dichrom_bright <- as.numeric(SSTraitDat3$Dichrom_bright)




#lets run a Phyloglm model!


################################## Run Models

phyglm_UN_brightness <- phyloglm( Urban ~ Dichrom_bright + Mass_log, 
                                   data = SSTraitDat3, 
                                   phy = SSphy3, 
                                   boot = 1000)
summary(phyglm_UN_brightness)


# run model with scaled and centered explanatory variables
set.seed(seednums[17])
#here
phyglm_UN_brightness_scale <- phyloglm( Urban ~ scale(Dichrom_bright) + scale(Mass_log), 
                                  data = SSTraitDat3, 
                                  phy = SSphy3, 
                                  boot = 1000) 
# this model converges
summary(phyglm_UN_brightness_scale)



# save model
saveRDS(phyglm_UN_brightness_scale, here("Outputs", "phyglm_UN_brightness_scale.rds"))
# load model
phyglm_UN_brightness_scale <- readRDS(here("Outputs", "phyglm_UN_brightness_scale.rds"))



# get alpha, t, and half life for the model
(phyglm_UN_brightness_scale$mean.tip.height) # t
(alpha_bright <- phyglm_UN_brightness_scale$alpha) # alpha
(hl_bright <- log(2)/alpha_bright) # half life
# compared to t, this is a small/moderate half life





######################## UN and Dichromatism - HUE ##########################

# lets first simplify a NEW database by removing records where we don't have an UN / brood_value
UNDataUT <- C_SSelect_dat2 %>% filter(!is.na(Urban)) 
SSData6 <- UNDataUT %>% filter(!is.na(Dichrom_hue)) 
length(SSData6$Dichrom_hue)
#61 species with UN and Dichrom_hue

colnames(SSData6)

###### add and pair tree

# add rownames to data
row.names(SSData6) <- SSData6$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

SSphydat6 <- treedata(tree_out,SSData6, sort=T)

SSphy6 <- SSphydat6$phy
SSTraitDat6 <- as.data.frame(SSphydat6$data)

str(SSTraitDat6)
length(SSTraitDat6$Dichrom_hue)
#61


### convert traits of interest to numeric

SSTraitDat6$Urban <- as.numeric(SSTraitDat6$Urban)
SSTraitDat6$Mass_log <- as.numeric(SSTraitDat6$Mass_log)
SSTraitDat6$Dichrom_hue <- as.numeric(SSTraitDat6$Dichrom_hue)



#lets run a Phyloglm model!


################################## Run Models 

phyglm_UN_hue <- phyloglm( Urban ~ Dichrom_hue + Mass_log, 
                                  data = SSTraitDat6, 
                                  phy = SSphy6, 
                                  boot = 1000) 
summary(phyglm_UN_hue)



# model with scaled and centered variables
set.seed(seednums[18])
#here
phyglm_UN_hue_scale <- phyloglm( Urban ~ scale(Dichrom_hue) + scale(Mass_log), 
                           data = SSTraitDat6, 
                           phy = SSphy6, 
                            boot = 1000)

summary(phyglm_UN_hue_scale) # this model converges successfully



# save model
saveRDS(phyglm_UN_hue_scale, here("Outputs", "phyglm_UN_hue_scale.rds"))
# load model
phyglm_UN_hue_scale <- readRDS(here("Outputs", "phyglm_UN_hue_scale.rds"))


# get alpha, t, and half life for the model
(phyglm_UN_hue_scale$mean.tip.height) # t
(alpha_hue <- phyglm_UN_hue_scale$alpha) # alpha
(hl_hue <- log(2)/alpha_hue) # half life for the model 
# small half life compared with t -> low phylogenetic signal



######################## UAI and sexual-selection Male ##########################

# lets first simplify a NEW database by removing records where we don't have an UN / brood_value
UNDataUT <- C_SSelect_dat2 %>% filter(!is.na(Urban)) 
SSData9 <- UNDataUT %>% filter(!is.na(sex.sel.m)) 
length(SSData9$sex.sel.m)
#129 species with UN and sex.sel.m

colnames(SSData9)

###### add and pair tree

# add rownames to data
row.names(SSData9) <- SSData9$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

SSphydat9 <- treedata(tree_out,SSData9, sort=T)

SSphy9 <- SSphydat9$phy
SSTraitDat9 <- as.data.frame(SSphydat9$data)

str(SSTraitDat9)
length(SSTraitDat9$sex.sel.m)
#129


### convert traits of interest to numeric

SSTraitDat9$Urban <- as.numeric(SSTraitDat9$Urban)
SSTraitDat9$Mass_log <- as.numeric(SSTraitDat9$Mass_log)
SSTraitDat9$sex.sel.m <- as.numeric(SSTraitDat9$sex.sel.m)



#lets run a Phyloglm model!


################################ Run Models

phyglm_UN_ssm <- phyloglm( Urban ~ sex.sel.m + Mass_log, 
                           data = SSTraitDat9, 
                           phy = SSphy9, 
                           boot = 1000)

summary(phyglm_UN_ssm)
# this model fails to converge



# run model with scaled and centered explanatory variables
phyglm_UN_ssm_scale <- phyloglm( Urban ~ scale(sex.sel.m) + scale(Mass_log), 
                           data = SSTraitDat9, 
                           phy = SSphy9, 
                           boot = 1000) 
# this also fails to converge
summary(phyglm_UN_ssm_scale)



# print AIC values for models with different upper bounds
# intervals of 0.1 from 0 up to 4
for (i in seq(0, 4, by = 0.1)) {
  print(phyloglm(Urban ~ scale(sex.sel.m) + scale(Mass_log), 
                 data = SSTraitDat9, 
                 phy = SSphy9,
                 log.alpha.bound = i)$aic)
}
# AIC values support models with larger values of alpha (low phylo signal)


# try fixing alpha at upper bounds
# give the model a little more searching space for alpha because AIC is actually lowest slightly below log.alpha.bounds = 4
exp(3.8)/(phyglm_UN_ssm_scale$mean.tip.height) # equals 0.45. Using this as start.alpha
set.seed(seednums[19])
#here
phyglm_UN_ssm_fix <- phyloglm( Urban ~ scale(sex.sel.m) + scale(Mass_log), 
                                    data = SSTraitDat9, 
                                    phy = SSphy9, 
                                     log.alpha.bound = 4,
                                     start.alpha = 0.45,
                                      boot = 1000)
summary(phyglm_UN_ssm_fix)
# this model converges
# alpha at upper boundary


# save model
saveRDS(phyglm_UN_ssm_fix, here("Outputs", "phyglm_UN_ssm_fix.rds"))
# load model
phyglm_UN_ssm_fix <- readRDS(here("Outputs", "phyglm_UN_ssm_fix.rds"))


# compare results with a non-phylogenetic logistic model
glm_UN_ssm <- logistf(Urban ~ scale(sex.sel.m) + scale(Mass_log), 
                      data = SSTraitDat9)
summary(glm_UN_ssm)



# get alpha, t, and half life for the model
(phyglm_UN_ssm_fix$mean.tip.height) # t
(alpha_ssm <- phyglm_UN_ssm_fix$alpha) # alpha
(hl_ssm <- log(2)/alpha_ssm) # half life  
# small half life relative to t



######################## UN and Sexual Selection - Female ##########################

# lets first simplify a NEW database by removing records where we don't have an UN / brood_value
UNDataUT <- C_SSelect_dat2 %>% filter(!is.na(Urban)) 
SSData12 <- UNDataUT %>% filter(!is.na(sex.sel.f)) 
length(SSData12$sex.sel.f)
#129 species with UN and sex.sel.f

colnames(SSData12)

###### add and pair tree

# add rownames to data
row.names(SSData12) <- SSData12$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

SSphydat12 <- treedata(tree_out,SSData12, sort=T)

SSphy12 <- SSphydat12$phy
SSTraitDat12 <- as.data.frame(SSphydat12$data)

str(SSTraitDat12)
length(SSTraitDat12$sex.sel.f)
#129


### convert traits of interest to numeric

SSTraitDat12$Urban <- as.numeric(SSTraitDat12$Urban)
SSTraitDat12$Mass_log <- as.numeric(SSTraitDat12$Mass_log)
SSTraitDat12$sex.sel.f <- as.numeric(SSTraitDat12$sex.sel.f)



#lets run a Phyloglm model!


######################################## Run Models

phyglm_UN_ssf <- phyloglm( Urban ~ sex.sel.f + Mass_log, 
                           data = SSTraitDat12, 
                           phy = SSphy12, 
                           boot = 1000) 
summary(phyglm_UN_ssf)
# model converges
# alpha at upper bounds


# run model with scaled and centered variables
phyglm_UN_ssf_scale <- phyloglm( Urban ~ scale(sex.sel.f) + Mass_log, 
                           data = SSTraitDat12, 
                           phy = SSphy12, 
                           boot = 1000) 

summary(phyglm_UN_ssf_scale) 
# this version fails to converge
# alpha at upper bounds


# print AIC values for models with different upper bounds
# intervals of 0.1 from 0 up to 4
for (i in seq(0, 4, by = 0.1)) {
  print(phyloglm(Urban ~ scale(sex.sel.f) + scale(Mass_log), 
                 data = SSTraitDat12, 
                 phy = SSphy12,
                 log.alpha.bound = i)$aic)
}
# AIC values support models with larger values of alpha (low phylo signal)


# try fixing alpha at upper bounds
# had to go slightly above log.alpha.bound of 4 to get models that would converge
set.seed(seednums[20])
#here
phyglm_UN_ssf_fix_4.05 <- phyloglm( Urban ~ scale(sex.sel.f) + scale(Mass_log), 
                               data = SSTraitDat12, 
                               phy = SSphy12, 
                               log.alpha.bound = 4.05,
                               start.alpha = 0.56,
                               boot = 1000
                          )

phyglm_UN_ssf_fix_4.1 <- phyloglm( Urban ~ scale(sex.sel.f) + scale(Mass_log), 
                                    data = SSTraitDat12, 
                                    phy = SSphy12, 
                                    log.alpha.bound = 4.1,
                                    start.alpha = 0.6,
                                    boot = 1000
)
summary(phyglm_UN_ssf_fix_4.05)
summary(phyglm_UN_ssf_fix_4.1)
# both model converge and produce similar estimates


# save model
saveRDS(phyglm_UN_ssf_fix_4.05, here("Outputs", "phyglm_UN_ssf_fix.rds"))
# load model
phyglm_UN_ssf_fix_4.05 <- readRDS(here("Outputs", "phyglm_UN_ssf_fix.rds"))


# look at non-phylogenetic model
glm_UN_ssf <- logistf(Urban ~ scale(sex.sel.f) + scale(Mass_log), 
                      data = SSTraitDat12)
summary(glm_UN_ssf)
# agrees with two models above


# get alpha, t, and half life for the model
(phyglm_UN_ssf_fix_4.05$mean.tip.height) # t
(alpha_ssf <- phyglm_UN_ssf_fix_4.05$alpha) # alpha
(hl_ssf <- log(2)/alpha_ssf) # half life
# small half life compared to t


#######################################################
###################### SOCIAL TRAITS ##################
#######################################################



###################### Prep

#load in "Coastal_Species_Social.rds" - since this contains coastal species and all diet trait variables :) 

C_Social_dat <- readRDS(here("Outputs", "Coastal_Species_Social.rds"))
str(C_Social_dat)

C_Social_dat2 <- C_Social_dat %>%
  mutate(Species_Jetz  = str_replace(Species_Jetz, " ", "_"))
str(C_Social_dat2)

C_Social_dat2$Urban <- ifelse(C_Social_dat2$Urban == "U", 1, 0)
View(C_Social_dat2)
colnames(C_Social_dat2)



######################## UN and territoriality ##########################

# lets first simplify a NEW database by removing records where we don't have an UN / brood_value
UNDataUT <- C_Social_dat2 %>% filter(!is.na(Urban)) 
SocialData3 <- UNDataUT %>% 
  filter(!is.na(territoriality)) %>% 
  filter(territoriality != 2)
length(SocialData3$territoriality)
#122 species with UN and territoriality

colnames(SocialData3)

###### add and pair tree

# add rownames to data
row.names(SocialData3) <- SocialData3$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Socialphydat3 <- treedata(tree_out,SocialData3, sort=T)

Socialphy3 <- Socialphydat3$phy
SocialTraitDat3 <- as.data.frame(Socialphydat3$data)

str(SocialTraitDat3)
length(SocialTraitDat3$territoriality)
#122


### convert traits of interest to numeric

SocialTraitDat3$Urban <- as.numeric(SocialTraitDat3$Urban)
SocialTraitDat3$Mass_log <- as.numeric(SocialTraitDat3$Mass_log)
SocialTraitDat3$territoriality <- as.numeric(SocialTraitDat3$territoriality)



#lets run a Phyloglm model!


################################# Run Models

phyglm_UN_territorial <- phyloglm( Urban ~ territoriality + Mass_log, 
                           data = SocialTraitDat3, 
                           phy = Socialphy3, 
                           boot = 1000) 

summary(phyglm_UN_territorial)


# run model with scaled and centered explanatory variables
set.seed(seednums[21])
#here
phyglm_UN_territorial_scale <- phyloglm( Urban ~ scale(territoriality) + scale(Mass_log), 
                                   data = SocialTraitDat3, 
                                   phy = Socialphy3, 
                                    boot = 1000) 

summary(phyglm_UN_territorial_scale)



# save model
saveRDS(phyglm_UN_territorial_scale, here("Outputs", "phyglm_UN_territorial_scale.rds"))
# load model
phyglm_UN_territorial_scale <- readRDS(here("Outputs", "phyglm_UN_territorial_scale.rds"))



# get alpha, t, and half life for the model
(phyglm_UN_territorial_scale$mean.tip.height) # t
(alpha_terr <- phyglm_UN_territorial_scale$alpha) # alpha
(hl_terr <- log(2)/alpha_terr) # half life
# small half life compared to t -> low phylogenetic signal




######################## UN and cooperative ##########################

# do not do --> too small of N for Cooperative = 1 (only N = 8)


