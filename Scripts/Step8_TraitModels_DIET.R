##########################################################
######### Coastal Birds Urban Tolerance Project ##########
##########################################################
# Step 8 Phylogenetic Trait Models - DIET
# Authors: Sarah L. Jennings, Emma M. Garrison
##########################################################
# The objective of this script is to run phylogenetic models using diet traits
# 4 diet traits: % invert, % vert, % plant/seeds, % fruit/nectar


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


#load in "Coastal_Species_Diet.rds" - this contains coastal species and all diet trait variables

C_Diet_dat <- readRDS(here("Outputs", "Coastal_Species_Diet.rds"))
str(C_Diet_dat)

C_Diet_dat2 <- C_Diet_dat %>%
  mutate(Species_Jetz  = str_replace(Species_Jetz, " ", "_"))
str(C_Diet_dat2)

C_Diet_dat2$Urban <- ifelse(C_Diet_dat2$Urban == "U", 1, 0)
View(C_Diet_dat2)
colnames(C_Diet_dat2)


######################## UAI and % Diet Invertebrates ##########################

# create a new data frame by removing species with no UAI value or that are missing Diet % Inv
UAI_DietInvert <- C_Diet_dat2 %>% filter(!is.na(aveUAI)) %>% 
  filter(!is.na(Diet.Inv)) %>% as.data.frame()
length(UAI_DietInvert$Diet.Inv)
#798 species with UAI and Diet Invert

###### add and pair tree

# add rownames to data
row.names(UAI_DietInvert) <- UAI_DietInvert$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_DietInvert_phydat <- geiger::treedata(tree_out, UAI_DietInvert, sort=T)

UAI_DietInvert_phy <- UAI_DietInvert_phydat$phy
UAI_DietInvert_dat <- as.data.frame(UAI_DietInvert_phydat$data)

str(UAI_DietInvert_dat)
length(UAI_DietInvert_dat$Diet.Inv)
#798

### convert traits of interest to numeric
UAI_DietInvert_dat$aveUAI <- as.numeric(UAI_DietInvert_dat$aveUAI)
UAI_DietInvert_dat$Mass_log <- as.numeric(UAI_DietInvert_dat$Mass_log)
UAI_DietInvert_dat$Diet.Inv <- as.numeric(UAI_DietInvert_dat$Diet.Inv)

# Run phylogenetic linear model for UAI
UAI_GLS_invert <- gls(aveUAI~ Diet.Inv + Mass_log, data = UAI_DietInvert_dat, 
                   correlation = corPagel(0.5, phy = UAI_DietInvert_phy, fixed=F, form = ~Species_Jetz), 
                   method = "ML") 


# model summary and results
summary(UAI_GLS_invert) 
confint(UAI_GLS_invert)

# model diagnostics
check_model(UAI_GLS_invert) 
qqnorm(resid(UAI_GLS_invert)) 
qqline(resid(UAI_GLS_invert)) 
hist(resid(UAI_GLS_invert)) 

# save model for easy retrieval 
saveRDS(UAI_GLS_invert, here("Models/UAI", "UAI_GLS_invert.rds"))

######################## MUTI and % Diet Invertebrates ##########################

# create a new data frame by removing species with no MUTI value or that are missing Diet % Inv
MUTI_DietInvert <- C_Diet_dat2 %>% filter(!is.na(MUTIscore)) %>% 
  filter(!is.na(Diet.Inv)) %>% as.data.frame()
length(MUTI_DietInvert$Diet.Inv)
# 130 species with MUTI and Diet Invert

###### add and pair tree

# add rownames to data
row.names(MUTI_DietInvert) <- MUTI_DietInvert$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_DietInvert_phydat <- geiger::treedata(tree_out, MUTI_DietInvert, sort=T)

MUTI_DietInvert_phy <- MUTI_DietInvert_phydat$phy
MUTI_DietInvert_dat <- as.data.frame(MUTI_DietInvert_phydat$data)

str(MUTI_DietInvert_dat)
length(MUTI_DietInvert_dat$Diet.Inv)
# 130

### convert traits of interest to numeric
MUTI_DietInvert_dat$MUTIscore <- as.numeric(MUTI_DietInvert_dat$MUTIscore)
MUTI_DietInvert_dat$Mass_log <- as.numeric(MUTI_DietInvert_dat$Mass_log)
MUTI_DietInvert_dat$Diet.Inv <- as.numeric(MUTI_DietInvert_dat$Diet.Inv)

# Run phylogenetic linear model for MUTI
MUTI_GLS_invert <- gls(MUTIscore ~ Diet.Inv + Mass_log, data = MUTI_DietInvert_dat, 
                      correlation = corPagel(0.5, phy = MUTI_DietInvert_phy, fixed=F, form = ~Species_Jetz), 
                      method = "ML") 


# model summary and results
summary(MUTI_GLS_invert) 
confint(MUTI_GLS_invert)
confint(MUTI_GLS_invert, level = 0.85)

# model diagnostics
check_model(MUTI_GLS_invert) 
qqnorm(resid(MUTI_GLS_invert)) 
qqline(resid(MUTI_GLS_invert)) 
hist(resid(MUTI_GLS_invert)) 

# save model for easy retrieval 
saveRDS(MUTI_GLS_invert, here("Models/MUTI", "MUTI_GLS_invert.rds"))

######################## UN and % Diet Invertebrates ##########################

# create a new data frame by removing species with no UN value or that are missing % diet invert
UN_DietInvert <- C_Diet_dat2 %>% filter(!is.na(Urban)) %>% 
  filter(!is.na(Diet.Inv)) %>%
  column_to_rownames(., var="Species_Jetz")
length(UN_DietInvert$Diet.Inv)
#129 species with UN and Diet Invert

###### add and pair tree
tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UN_DietInvert_phydat <- treedata(tree_out, UN_DietInvert, sort=T)

UN_DietInvert_phy <- UN_DietInvert_phydat$phy
UN_DietInvert_dat <- as.data.frame(UN_DietInvert_phydat$data)

str(UN_DietInvert_dat)
length(UN_DietInvert_dat$Diet.Inv)
#129

### convert traits of interest to numeric
UN_DietInvert_dat$Urban <- as.numeric(UN_DietInvert_dat$Urban)
UN_DietInvert_dat$Mass_log <- as.numeric(UN_DietInvert_dat$Mass_log)
UN_DietInvert_dat$Diet.Inv <- as.numeric(UN_DietInvert_dat$Diet.Inv)


# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
phyglm_UN_Invert_scale <- phyloglm( Urban ~ scale(Diet.Inv) + scale(Mass_log), 
                                    data = UN_DietInvert_dat, 
                                    phy = UN_DietInvert_phy, 
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
                 data = UN_DietInvert_dat, 
                 phy = UN_DietInvert_phy,
                 log.alpha.bound = i)$aic)
}

# higher support for larger alpha values (later in the printed list) that correspond to low phylogenetic signal


# try to fix alpha at exp(4)/t (at the default upper bounds)
# this fails to converge
phyglm_UN_Invert_fix_4 <- phyloglm( Urban ~ scale(Diet.Inv) + scale(Mass_log), 
                                    data = UN_DietInvert_dat, 
                                    phy = UN_DietInvert_phy,
                                    start.alpha = 0.55,
                                    log.alpha.bound = 4)

# try increasing upper bounds a small amount over 4
# both these models converge
set.seed(568)
phyglm_UN_Invert_fix_4.05 <- phyloglm( Urban ~ scale(Diet.Inv) + scale(Mass_log), 
                                       data = UN_DietInvert_dat, 
                                       phy = UN_DietInvert_phy, 
                                       start.alpha = 0.55,
                                       log.alpha.bound = 4.05, boot=1000)

summary(phyglm_UN_Invert_fix_4.05) # this converges

phyglm_UN_Invert_fix_4.1 <- phyloglm( Urban ~ scale(Diet.Inv) + scale(Mass_log), 
                                      data = UN_DietInvert_dat, 
                                      phy = UN_DietInvert_phy,
                                      start.alpha = 0.55,
                                      log.alpha.bound = 4.1, boot=1000)
summary(phyglm_UN_Invert_fix_4.1) # also converges

# coefficients fairly stable across 4.05 and 4.1 models
# save model with log.alpha.bound = 4.05
saveRDS(phyglm_UN_Invert_fix_4.05, here("Models/UN", "phyglm_UN_Invert_fix.rds"))
# load model
phyglm_UN_Invert_fix_4.05 <- readRDS(here("Models/UN", "phyglm_UN_Invert_fix.rds"))

# as alpha is at upper bounds, we can also compare results to a non-phylogenetic logistic model
glm_UN_Invert <- logistf(Urban ~ scale(Diet.Inv) + scale(Mass_log), 
                         data = UN_DietInvert)
summary(glm_UN_Invert)
# we reach similar conclusions as two models above

# get alpha, t, and half life for the model
# using phyglm_UN_Invert_fix_4.05 as final model
(phyglm_UN_Invert_fix_4.05$mean.tip.height) # t
(alpha_invert <- phyglm_UN_Invert_fix_4.05$alpha) # alpha
(hl_invert <- log(2)/alpha_invert) # half life
# small half life -> low phylogenetic signal


#############################################################################
#############################################################################

######################## UAI and % Diet Vertebrates #########################

# create a new data frame by removing species with no UAI value or that are missing % diet vert
UAI_DietVert <- C_Diet_dat2 %>% filter(!is.na(aveUAI))  %>% 
  filter(!is.na(Diet.Vert)) %>% as.data.frame()
length(UAI_DietVert$Diet.Vert)
#798 species with UAI and Diet Vert

###### add and pair tree

# add rownames to data
row.names(UAI_DietVert) <- UAI_DietVert$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_DietVert_phydat <- treedata(tree_out, UAI_DietVert, sort=T)

UAI_DietVert_phy <- UAI_DietVert_phydat$phy
UAI_DietVert_dat <- as.data.frame(UAI_DietVert_phydat$data)

str(UAI_DietVert_dat)
length(UAI_DietVert_dat$Diet.Vert)
#798

### convert traits of interest to numeric
UAI_DietVert_dat$aveUAI <- as.numeric(UAI_DietVert_dat$aveUAI)
UAI_DietVert_dat$Mass_log <- as.numeric(UAI_DietVert_dat$Mass_log)
UAI_DietVert_dat$Diet.Vert <- as.numeric(UAI_DietVert_dat$Diet.Vert)

# Run phylogenetic linear model for UAI
UAI_GLS_vert <- gls(aveUAI~ Diet.Vert + Mass_log, data = UAI_DietVert_dat, 
                      correlation = corPagel(0.5, phy = UAI_DietVert_phy, fixed=F, form = ~Species_Jetz), 
                      method = "ML") 

# model summary and results 
summary(UAI_GLS_vert) 
confint(UAI_GLS_vert)

# model diagnostics
check_model(UAI_GLS_vert) 
qqnorm(resid(UAI_GLS_vert)) 
qqline(resid(UAI_GLS_vert)) 
hist(resid(UAI_GLS_vert)) 

# save model for easy retrieval 
saveRDS(UAI_GLS_vert, here("Models/UAI", "UAI_GLS_vert.rds"))

######################## MUTI and % Diet Vertebrates ##########################


# create a new data frame by removing species with no MUTI value or that are missing % diet vert
MUTI_DietVert <- C_Diet_dat2 %>% filter(!is.na(MUTIscore))  %>% 
  filter(!is.na(Diet.Vert)) %>% as.data.frame()
length(MUTI_DietVert$Diet.Vert)
# 130 species with MUTI and Diet Vert

###### add and pair tree

# add rownames to data
row.names(MUTI_DietVert) <- MUTI_DietVert$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_DietVert_phydat <- treedata(tree_out, MUTI_DietVert, sort=T)

MUTI_DietVert_phy <- MUTI_DietVert_phydat$phy
MUTI_DietVert_dat <- as.data.frame(MUTI_DietVert_phydat$data)

str(MUTI_DietVert_dat)
length(MUTI_DietVert_dat$Diet.Vert)
#798

### convert traits of interest to numeric
MUTI_DietVert_dat$MUTIscore <- as.numeric(MUTI_DietVert_dat$MUTIscore)
MUTI_DietVert_dat$Mass_log <- as.numeric(MUTI_DietVert_dat$Mass_log)
MUTI_DietVert_dat$Diet.Vert <- as.numeric(MUTI_DietVert_dat$Diet.Vert)

# Run phylogenetic linear model for MUTI
MUTI_GLS_vert <- gls(MUTIscore ~ Diet.Vert + Mass_log, data = MUTI_DietVert_dat, 
                    correlation = corPagel(0.5, phy = MUTI_DietVert_phy, fixed=F, form = ~Species_Jetz), 
                    method = "ML") 

# model summary and results 
summary(MUTI_GLS_vert) 
confint(MUTI_GLS_vert)

# model diagnostics
check_model(MUTI_GLS_vert) 
qqnorm(resid(MUTI_GLS_vert)) 
qqline(resid(MUTI_GLS_vert)) 
hist(resid(MUTI_GLS_vert)) 

# save model for easy retrieval 
saveRDS(MUTI_GLS_vert, here("Models/MUTI", "MUTI_GLS_vert.rds"))

######################## UN and % Diet Vertebrates ##########################

# create a new data frame by removing species with no UN value or that are missing % diet vert
UN_DietVert <- C_Diet_dat2 %>% filter(!is.na(Urban)) %>% 
  filter(!is.na(Diet.Vert)) %>%
  column_to_rownames(., var="Species_Jetz")
length(UN_DietVert$Diet.Vert)
#129 species with UN and Diet Vert

###### add and pair tree
tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UN_DietVert_phydat <- treedata(tree_out, UN_DietVert, sort=T)

UN_DietVert_phy <- UN_DietVert_phydat$phy
UN_DietVert_dat <- as.data.frame(UN_DietVert_phydat$data)

str(UN_DietVert_dat)
length(UN_DietVert_dat$Diet.Vert)
#129

### convert traits of interest to numeric
UN_DietVert_dat$Urban <- as.numeric(UN_DietVert_dat$Urban)
UN_DietVert_dat$Mass_log <- as.numeric(UN_DietVert_dat$Mass_log)
UN_DietVert_dat$Diet.Vert <- as.numeric(UN_DietVert_dat$Diet.Vert)


# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
phyglm_UN_Vert_scale <- phyloglm( Urban ~ scale(Diet.Vert) + scale(Mass_log), 
                                  data = UN_DietVert_dat, 
                                  phy = UN_DietVert_phy, 
                                  boot = 1000)
summary(phyglm_UN_Vert_scale)
# this fails to converge

# print AIC values for models with different upper bounds
# intervals of 0.1 from 0 up to 4
for (i in seq(0, 4, by = 0.1)) {
  print(phyloglm(Urban ~ scale(Diet.Vert) + scale(Mass_log), 
                 data = UN_DietVert_dat, 
                 phy = UN_DietVert_phy, 
                 log.alpha.bound = i)$aic)
}
# AIC values support models with larger values of alpha (low phylo signal)

# try to fix alpha at upper bounds
phyglm_UN_Vert_fix <- phyloglm( Urban ~ scale(Diet.Vert) + scale(Mass_log), 
                                data = UN_DietVert_dat, 
                                phy = UN_DietVert_phy, 
                                start.alpha = 0.55,
                                log.alpha.bound = 4, boot = 1000)
summary(phyglm_UN_Vert_fix) # this fails to converge



# try increasing upper bounds a small amount over 4
# both these models converge
set.seed(351)
phyglm_UN_Vert_fix_4.05 <- phyloglm( Urban ~ scale(Diet.Vert) + scale(Mass_log), 
                                       data = UN_DietVert_dat, 
                                       phy = UN_DietVert_phy,  
                                       start.alpha = 0.55,
                                       log.alpha.bound = 4.05, boot=1000)

summary(phyglm_UN_Vert_fix_4.05) # this converges

phyglm_UN_Vert_fix_4.1 <- phyloglm( Urban ~ scale(Diet.Vert) + scale(Mass_log), 
                                      data = UN_DietVert_dat, 
                                      phy = UN_DietVert_phy, 
                                      start.alpha = 0.55,
                                      log.alpha.bound = 4.1, boot=1000)
summary(phyglm_UN_Vert_fix_4.1) # also converges


# coefficients fairly stable across 4.05 and 4.1 models
# save model with log.alpha.bound = 4.05
saveRDS(phyglm_UN_Vert_fix_4.05, here("Models/UN", "phyglm_UN_Vert_fix.rds"))
# load model
phyglm_UN_Vert_fix_4.05 <- readRDS(here("Models/UN", "phyglm_UN_Vert_fix.rds"))

# as alpha is at upper bound, look at a non-phylogenetic model
glm_UN_Vert <- logistf(Urban ~ scale(Diet.Vert) + scale(Mass_log), 
                       data = UN_DietVert)
summary(glm_UN_Vert)
# very similar

# get alpha, t, and half life for the model
(phyglm_UN_Vert_fix$mean.tip.height) # t
(alpha_vert <- phyglm_UN_Vert_fix$alpha) # alpha
(hl_vert <- log(2)/alpha_vert) # half life
# compared to t, this is a small half life


#############################################################################
#############################################################################

######################## UAI and % Diet Plant/Seed ##########################

# create a new data frame by removing species with no UAI value or that are missing Diet % Plant/Seed
UAI_DietPS <- C_Diet_dat2 %>% filter(!is.na(aveUAI))  %>% 
  filter(!is.na(Diet.PS)) %>% as.data.frame()
length(UAI_DietPS$Diet.PS)
#798 species with UAI and Diet Plant/Seed

###### add and pair tree

# add rownames to data
row.names(UAI_DietPS) <- UAI_DietPS$Species_Jetz

tree_out <- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_DietPS_phydat <- treedata(tree_out, UAI_DietPS, sort=T)

UAI_DietPS_phy <- UAI_DietPS_phydat$phy
UAI_DietPS_dat <- as.data.frame(UAI_DietPS_phydat$data)

str(UAI_DietPS_dat)
length(UAI_DietPS_dat$Diet.PS)
#798

### convert traits of interest to numeric
UAI_DietPS_dat$aveUAI <- as.numeric(UAI_DietPS_dat$aveUAI)
UAI_DietPS_dat$Mass_log <- as.numeric(UAI_DietPS_dat$Mass_log)
UAI_DietPS_dat$Diet.PS <- as.numeric(UAI_DietPS_dat$Diet.PS)

# run phylogenetic linear model for UAI using gls()
UAI_GLS_PS <- gls(aveUAI ~ Diet.PS + Mass_log, data = UAI_DietPS_dat, 
                     correlation = corPagel(0.5, phy = UAI_DietPS_phy, fixed=F, form = ~Species_Jetz), 
                     method = "ML") 

# model summary and results
summary(UAI_GLS_PS) 
confint(UAI_GLS_PS)

# model diagnostics
check_model(UAI_GLS_PS) 
qqnorm(resid(UAI_GLS_PS)) 
qqline(resid(UAI_GLS_PS)) 
hist(resid(UAI_GLS_PS)) 

# save model for easy retrieval 
saveRDS(UAI_GLS_PS, here("Models/UAI", "UAI_GLS_PS.rds"))

######################## MUTI and % Diet Plant/Seed ##########################

# create a new data frame by removing species with no MUTI value or that are missing Diet % Plant/Seed
MUTI_DietPS <- C_Diet_dat2 %>% filter(!is.na(MUTIscore))  %>% 
  filter(!is.na(Diet.PS)) %>% as.data.frame()
length(MUTI_DietPS$Diet.PS)
# 130 species with MUTI and Diet Plant/Seed

###### add and pair tree

# add rownames to data
row.names(MUTI_DietPS) <- MUTI_DietPS$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_DietPS_phydat <- treedata(tree_out, MUTI_DietPS, sort=T)

MUTI_DietPS_phy <- MUTI_DietPS_phydat$phy
MUTI_DietPS_dat <- as.data.frame(MUTI_DietPS_phydat$data)

str(MUTI_DietPS_dat)
length(MUTI_DietPS_dat$Diet.PS)
#130

### convert traits of interest to numeric
MUTI_DietPS_dat$MUTIscore <- as.numeric(MUTI_DietPS_dat$MUTIscore)
MUTI_DietPS_dat$Mass_log <- as.numeric(MUTI_DietPS_dat$Mass_log)
MUTI_DietPS_dat$Diet.PS <- as.numeric(MUTI_DietPS_dat$Diet.PS)


# run phylogenetic linear model for MUTI using gls()
MUTI_GLS_PS <- gls(MUTIscore ~ Diet.PS + Mass_log, data = MUTI_DietPS_dat, 
                  correlation = corPagel(0.5, phy = MUTI_DietPS_phy, fixed=F, form = ~Species_Jetz), 
                  method = "ML") 

# model summary and results
summary(MUTI_GLS_PS) 
confint(MUTI_GLS_PS)

#check out the model
check_model(MUTI_GLS_PS) 
qqnorm(resid(MUTI_GLS_PS)) 
qqline(resid(MUTI_GLS_PS)) 
hist(resid(MUTI_GLS_PS)) 

# save model for easy retrieval 
saveRDS(MUTI_GLS_PS, here("Models/MUTI", "MUTI_GLS_PS.rds"))


######################## UN and % Diet Plant/Seed ##########################

# create a new data frame by removing species with no UN value or that are missing Diet % Plant/Seed
UN_DietPS <- C_Diet_dat2 %>% filter(!is.na(Urban)) %>%
  filter(!is.na(Diet.PS)) %>% column_to_rownames(., var = "Species_Jetz")
length(UN_DietPS$Diet.PS)
#129 species with UN and Diet Plant/Seed

###### add and pair tree
tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UN_DietPS_phydat <- treedata(tree_out, UN_DietPS, sort=T)

UN_DietPS_phy <- UN_DietPS_phydat$phy
UN_DietPS_dat <- as.data.frame(UN_DietPS_phydat$data)

str(UN_DietPS_dat)
length(UN_DietPS_dat$Diet.PS)
#129

### convert traits of interest to numeric
UN_DietPS_dat$Urban <- as.numeric(UN_DietPS_dat$Urban)
UN_DietPS_dat$Mass_log <- as.numeric(UN_DietPS_dat$Mass_log)
UN_DietPS_dat$Diet.PS <- as.numeric(UN_DietPS_dat$Diet.PS)

# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
set.seed(382)
phyglm_UN_PS_scale <- phyloglm( Urban ~ scale(Diet.PS) + scale(Mass_log), 
                                data = UN_DietPS_dat, 
                                phy = UN_DietPS_phy, 
                                boot = 1000) 

# this converges
summary(phyglm_UN_PS_scale)

# save model
saveRDS(phyglm_UN_PS_scale, here("Models/UN", "phyglm_UN_PS_scale.rds"))
# load model
phyglm_UN_PS_scale <- readRDS(here("Models/UN", "phyglm_UN_PS_scale.rds"))

# get alpha, t, and half life for the model
(phyglm_UN_PS_scale$mean.tip.height) # t
(alpha_PS <- phyglm_UN_PS_scale$alpha) # alpha
(hl_PS <- log(2)/alpha_PS) # half life
# small compared to t -> low phylogenetic signal


#############################################################################
#############################################################################

######################## UAI and % Diet Fruit/Nectar ##########################

# create a new data frame by removing species with no UAI value or that are missing Diet % Fruit/Nectar
UAI_DietFN <- C_Diet_dat2 %>% filter(!is.na(aveUAI))  %>% 
  filter(!is.na(Diet.FN)) %>% as.data.frame()
length(UAI_DietFN$Diet.FN)
#798 species with UAI and and Diet Fruit Nectar

###### add and pair tree

# add rownames to data
row.names(UAI_DietFN) <- UAI_DietFN$Species_Jetz

# import tree
tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_DietFN_phydat <- treedata(tree_out, UAI_DietFN, sort=T)

UAI_DietFN_phy <- UAI_DietFN_phydat$phy
UAI_DietFN_dat <- as.data.frame(UAI_DietFN_phydat$data)

str(UAI_DietFN_dat)
length(UAI_DietFN_dat$Diet.FN)
#798

### convert traits of interest to numeric
UAI_DietFN_dat$aveUAI <- as.numeric(UAI_DietFN_dat$aveUAI)
UAI_DietFN_dat$Mass_log <- as.numeric(UAI_DietFN_dat$Mass_log)
UAI_DietFN_dat$Diet.FN <- as.numeric(UAI_DietFN_dat$Diet.FN)


# run phylogenetic linear model for UAI using gls()
UAI_GLS_FN <- gls(aveUAI ~ Diet.FN + Mass_log, data = UAI_DietFN_dat, 
                   correlation = corPagel(0.5, phy = UAI_DietFN_phy, fixed=F, form = ~Species_Jetz), 
                   method = "ML") 


# model summary and results
summary(UAI_GLS_FN) 
confint(UAI_GLS_FN)

# model diagnostics
check_model(UAI_GLS_FN) 
qqnorm(resid(UAI_GLS_FN)) 
qqline(resid(UAI_GLS_FN)) 
hist(resid(UAI_GLS_FN))

# save model for easy retrieval 
saveRDS(UAI_GLS_FN, here("Models/UAI", "UAI_GLS_FN.rds"))


######################## MUTI and % Diet Fruit/Nectar ##########################


# create a new data frame by removing species with no MUTI value or that are missing Diet % Fruit/Nectar
MUTI_DietFN <- C_Diet_dat2 %>% filter(!is.na(MUTIscore)) %>% 
  filter(!is.na(Diet.FN)) %>% as.data.frame()
length(MUTI_DietFN$Diet.FN)
#130 species with MUTI and Diet Fruit Nectar

###### add and pair tree

# add rownames to data
row.names(MUTI_DietFN) <- MUTI_DietFN$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_DietFN_phydat <- treedata(tree_out, MUTI_DietFN, sort=T)

MUTI_DietFN_phy <- MUTI_DietFN_phydat$phy
MUTI_DietFN_dat <- as.data.frame(MUTI_DietFN_phydat$data)

str(MUTI_DietFN_dat)
length(MUTI_DietFN_dat$Diet.FN)
#798

### convert traits of interest to numeric
MUTI_DietFN_dat$MUTIscore <- as.numeric(MUTI_DietFN_dat$MUTIscore)
MUTI_DietFN_dat$Mass_log <- as.numeric(MUTI_DietFN_dat$Mass_log)
MUTI_DietFN_dat$Diet.FN <- as.numeric(MUTI_DietFN_dat$Diet.FN)


# run phylogenetic linear model for MUTI using gls()
MUTI_GLS_FN <- gls(MUTIscore ~ Diet.FN + Mass_log, data = MUTI_DietFN_dat, 
                  correlation = corPagel(0.5, phy = MUTI_DietFN_phy,fixed=F, form = ~Species_Jetz), 
                  method = "ML") 

# model summary and results
summary(MUTI_GLS_FN) 
confint(MUTI_GLS_FN)

# model diagnostics
check_model(MUTI_GLS_FN) 
qqnorm(resid(MUTI_GLS_FN)) 
qqline(resid(MUTI_GLS_FN)) 
hist(resid(MUTI_GLS_FN)) 

# save model for easy retrieval 
saveRDS(MUTI_GLS_FN, here("Models/MUTI", "MUTI_GLS_FN.rds"))

######################## UN and % Diet Fruit/Nectar ##########################

# create a new data frame by removing species with no UN value or that are missing % Diet Fruit/Nectar
UN_DietFN <- C_Diet_dat2 %>% filter(!is.na(Urban))  %>% filter(!is.na(Diet.FN)) %>%
  column_to_rownames(., var="Species_Jetz")
length(UN_DietFN$Diet.FN)
#129 species with UN and Diet Fruit/Nectar

###### add and pair tree

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UN_DietFN_phydat <- treedata(tree_out, UN_DietFN, sort=T)

UN_DietFN_phy <- UN_DietFN_phydat$phy
UN_DietFN_dat <- as.data.frame(UN_DietFN_phydat$data)

str(UN_DietFN_dat)
length(UN_DietFN_dat$Diet.FN)
#129

### convert traits of interest to numeric
UN_DietFN_dat$Urban <- as.numeric(UN_DietFN_dat$Urban)
UN_DietFN_dat$Mass_log <- as.numeric(UN_DietFN_dat$Mass_log)
UN_DietFN_dat$Diet.FN <- as.numeric(UN_DietFN_dat$Diet.FN)

# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
set.seed(380)
phyglm_UN_FN_scale <- phyloglm( Urban ~ scale(Diet.FN) + scale(Mass_log), 
                                data = UN_DietFN_dat, 
                                phy = UN_DietFN_phy, 
                                boot = 1000)

# this model converges
summary(phyglm_UN_FN_scale)

# save model
saveRDS(phyglm_UN_FN_scale, here("Models/UN", "phyglm_UN_FN_scale.rds"))
# load model
phyglm_UN_FN_scale <- readRDS(here("Models/UN", "phyglm_UN_FN_scale.rds"))


# get alpha, t, and half life for the model
(phyglm_UN_FN_scale$mean.tip.height) # t
(alpha_FN <- phyglm_UN_FN_scale$alpha) # alpha
(hl_FN <- log(2)/alpha_FN) # half life
# compared to t, this is a small half life
