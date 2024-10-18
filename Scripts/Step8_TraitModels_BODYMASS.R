##########################################################
######### Coastal Birds Urban Tolerance Project ##########
##########################################################
# Step 8 Phylogenetic Trait Models - BODY MASS
# Authors: Sarah L. Jennings, Emma M. Garrison
##########################################################
# The objective of this script is to run phylogenetic models using JUST body mass


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

###################### Prep

#load in "Coastal_Species_w_Mass.rds" - this contains coastal species and body mass

C_mass_dat <- readRDS(here("Outputs", "Coastal_Species_w_Mass.rds"))
str(C_mass_dat)

C_mass_dat2 <- C_mass_dat %>%
  mutate(Species_Jetz  = str_replace(Species_Jetz, " ", "_"))
str(C_mass_dat2)

C_mass_dat2$Urban <- ifelse(C_mass_dat2$Urban == "U", 1, 0)
View(C_mass_dat2)
colnames(C_mass_dat2)


######################## UAI and body mass ##########################

# create a new data frame by removing species with no UAI value or that are missing body mass
UAI_Mass_Data <- C_mass_dat2 %>% filter(!is.na(aveUAI)) %>% filter(!is.na(Mass_log)) 
length(UAI_Mass_Data$Mass_log)
#798 species with UAI and Mass_log

colnames(MassData1)

###### add and pair tree

# add rownames to data
row.names(UAI_Mass_Data) <- UAI_Mass_Data$Species_Jetz

# import tree
tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_Mass_phydat <- treedata(tree_out, UAI_Mass_Data, sort=T)

UAI_Mass_phy <- UAI_Mass_phydat$phy
UAI_Mass_dat <- as.data.frame(UAI_Mass_phydat$data)

str(UAI_Mass_dat)
length(UAI_Mass_dat$Mass_log)
#798


### convert traits of interest to numeric
UAI_Mass_dat$aveUAI <- as.numeric(UAI_Mass_dat$aveUAI)
UAI_Mass_dat$Mass_log <- as.numeric(UAI_Mass_dat$Mass_log)


# Run phylogenetic linear model for UAI
UAI_GLS_mass <- gls(aveUAI~ Mass_log, data = UAI_Mass_dat, 
                         correlation = corPagel(0.5, phy= UAI_Mass_phy,fixed=F, form = ~Species_Jetz), 
                         method = "ML") 

# model summary
summary(UAI_GLS_mass) 
confint(UAI_GLS_mass) # 95% CI
confint(UAI_GLS_mass, level = 0.85) # 85% CI

# model diagnostics
check_model(UAI_GLS_mass) 
qqnorm(resid(UAI_GLS_mass)) 
qqline(resid(UAI_GLS_mass)) 
hist(resid(UAI_GLS_mass)) 

# save model for easy retrieval 
saveRDS(UAI_GLS_mass, here("Models/UAI", "UAI_GLS_mass.rds"))

######################## MUTI and body mass ##########################

# create a new data frame by removing species with no MUTI value or that are missing body mass
MUTI_Mass_Data <- C_mass_dat2 %>% filter(!is.na(MUTIscore)) %>% filter(!is.na(Mass_log)) 
length(MUTI_Mass_Data$Mass_log)
# 130 species with MUTI and Mass_log

###### add and pair tree

# add rownames to data
row.names(MUTI_Mass_Data) <- MUTI_Mass_Data$Species_Jetz

# import tree
tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_Mass_phydat <- treedata(tree_out, MUTI_Mass_Data, sort=T)

MUTI_Mass_phy <- MUTI_Mass_phydat$phy
MUTI_Mass_dat <- as.data.frame(MUTI_Mass_phydat$data)

str(MUTI_Mass_dat)
length(MUTI_Mass_dat$Mass_log)
#130

### convert traits of interest to numeric
MUTI_Mass_dat$MUTIscore <- as.numeric(MUTI_Mass_dat$MUTIscore)
MUTI_Mass_dat$Mass_log <- as.numeric(MUTI_Mass_dat$Mass_log)

# Run phylogenetic linear model for MUTI
MUTI_GLS_mass <- gls(MUTIscore~ Mass_log, data = MUTI_Mass_dat, 
                    correlation = corPagel(0.5, phy = MUTI_Mass_phy, fixed=F, form = ~Species_Jetz), 
                    method = "ML") 

# model summary and results
summary(MUTI_GLS_mass) 
confint(MUTI_GLS_mass)
confint(MUTI_GLS_mass, level = 0.85)

# model diagnostics
check_model(MUTI_GLS_mass) 
qqnorm(resid(MUTI_GLS_mass)) 
qqline(resid(MUTI_GLS_mass)) 
hist(resid(MUTI_GLS_mass)) 

# save model for easy retrieval 
saveRDS(MUTI_GLS_mass, here("Models/MUTI", "MUTI_GLS_mass.rds"))

######################## UN and body mass ##########################


# create a new data frame by removing species with no UN value or that are missing body mass
UN_Mass_Data <- C_mass_dat2 %>% filter(!is.na(Urban)) %>%
  filter(!is.na(Mass_log)) %>%
  column_to_rownames(., var="Species_Jetz")
length(UN_Mass_Data$Mass_log)
#129 species with UN and Mass_log

###### add and pair tree

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UN_Mass_phydat <- treedata(tree_out, UN_Mass_Data, sort=T)

UN_Mass_phy <- UN_Mass_phydat$phy
UN_Mass_dat <- as.data.frame(UN_Mass_phydat$data)

str(UN_Mass_dat)
length(UN_Mass_dat$Mass_log)
#129


### convert traits of interest to numeric
UN_Mass_dat$Urban <- as.numeric(UN_Mass_dat$Urban)
UN_Mass_dat$Mass_log <- as.numeric(UN_Mass_dat$Mass_log)

# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
set.seed(309)
phyglm_UN_Mass_scale <- phyloglm( Urban ~ scale(Mass_log), 
                                  data = UN_Mass_dat, 
                                  phy = UN_Mass_phy, boot=1000) 

summary(phyglm_UN_Mass_scale)
# some of the bootstrapped models fail to converge

# we will try to determine which values of alpha produce the best fit models using AIC
# to do this, we can print AIC values for models with different upper bounds
# using intervals of 0.1 from 0 up to 4 for log.alpha.bound
# 4 is the default setting for the upper bounds of log.alpha.bound
for (i in seq(0, 4, by = 0.1)) {
  print(phyloglm(Urban ~ scale(Mass_log), 
                 data = UN_Mass_dat, 
                 phy = UN_Mass_phy,
                 log.alpha.bound = i)$aic)
}

# there is higher support for larger alpha values (later in the printed list) that correspond to low phylogenetic signal

# Next, we try to fit a bootstrapped model that converges by fixing alpha at upper bounds
# the upper bounds for alpha in our models is:
t <- phyglm_UN_Mass_scale$mean.tip.height # set the mean tip height of the tree as a variable called t
t # t for our tree is 97.5606

# use t to find the upper bound of alpha for our model
exp(4)/t # equals 0.5596322
exp(-4)/t # this is the lower bound = 0.0001877
# Note: alpha can take on values from 0 to infinity, but around log.alpha.bound = 4 the model should produce results similar to a non-phylogenetic logistic model

# using log.alpha.bound = 4 in phyloglm() sets the upper bound at this value (0.5596332)
# we can also set start.alpha, the alpha value where the model begins its search
# we can't set start.alpha to exactly 0.5596332 and have the run model (so it won't let us actually fix the value of alpha)
# instead, we can give a start.alpha value very close to the upper bound (e.g., 0.557)
# this will constrains the model's search area for the optimal alpha value within a very small range of values (similar to fixing it)
set.seed(848)
phyglm_UN_Mass_fix <- phyloglm(Urban ~ scale(Mass_log), 
                               data = UN_Mass_dat, 
                               phy = UN_Mass_phy, 
                               start.alpha = 0.557,
                               log.alpha.bound = 4, 
                               boot=1000)
summary(phyglm_UN_Mass_fix)

# save model for easy retrieval 
saveRDS(phyglm_UN_Mass_fix, here("Models/UN", "phyglm_UN_Mass_fix.rds"))
# load model
phyglm_UN_Mass_fix <- readRDS(here("Models/UN", "phyglm_UN_Mass_fix.rds"))

# as alpha is at upper bounds, we can also look at a regular logistic model
# we will use logistf package for this because this runs a logistic regression with Firth's correction (a bias reduction method)
# Ives and Garland 2010 recommend log.alpha.bound = 4 as the limit (and this is default setting in phyloglm)
# they specify 4 because when the model reaches the upper bounds of this limit,
# the model estimates should be very similar to a logistic model using Firth's correction
# therefore, this next step should produce results similar to phyglm_UN_Mass_fix
glm_UN_Mass <- logistf(Urban ~ scale(Mass_log), data = UN_Mass_Data)
summary(glm_UN_Mass)
summary(phyglm_UN_Mass_fix) # compare estimates
# very similar coefficients

# get values for results table for the fixed alpha model with scale(Mass_log)
summary(phyglm_UN_Mass_fix)
# use the bootstrapped CIs in the summary

# get alpha, t, and the half life for the model
(phyglm_UN_Mass_fix$mean.tip.height) # this is t (mean tip height) of the tree
(alpha_mass <- phyglm_UN_Mass_fix$alpha) # this is alpha
(hl_mass <- log(2)/alpha_mass) # this is the half-life for the model
# compare the value for half life with the mean tip height of the tree
# compared to t, the half life is small -> therefore we conclude there is low phylogenetic signal
