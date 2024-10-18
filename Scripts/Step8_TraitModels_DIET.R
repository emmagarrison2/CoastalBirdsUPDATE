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
UAIDataUT <- C_Diet_dat2 %>% filter(!is.na(aveUAI)) 
DietData1 <- UAIDataUT %>% filter(!is.na(Diet.Inv)) 
length(DietData1$Diet.Inv)
#798 species with UAI and CT

###### add and pair tree

# add rownames to data
row.names(DietData1) <- DietData1$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Dietphydat1 <- geiger::treedata(tree_out,DietData1, sort=T)

Dietphy1 <- Dietphydat1$phy
DietTraitDat1 <- as.data.frame(Dietphydat1$data)

str(DietTraitDat1)
length(DietTraitDat1$Diet.Inv)
#798

### convert traits of interest to numeric
DietTraitDat1$aveUAI <- as.numeric(DietTraitDat1$aveUAI)
DietTraitDat1$Mass_log <- as.numeric(DietTraitDat1$Mass_log)
DietTraitDat1$Diet.Inv <- as.numeric(DietTraitDat1$Diet.Inv)

# Run phylogenetic linear model for UAI
UAI_GLS_invert <- gls(aveUAI~ Diet.Inv + Mass_log, data = DietTraitDat1, 
                   correlation = corPagel(0.5, phy=Dietphy1,fixed=F, form = ~Species_Jetz), 
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

# create a new data frame by removing species with no MUTI value or that are missing % diet invert
MUTIDataUT <- C_Diet_dat2 %>% filter(!is.na(MUTIscore)) 
DietData2 <- MUTIDataUT %>% filter(!is.na(Diet.Inv)) 
length(DietData2$Diet.Inv)
#798 species with UAI and CT

###### add and pair tree

# add rownames to data
row.names(DietData2) <- DietData2$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Dietphydat2 <- geiger::treedata(tree_out,DietData2, sort=T)

Dietphy2 <- Dietphydat2$phy
DietTraitDat2 <- as.data.frame(Dietphydat2$data)

str(DietTraitDat2)
length(DietTraitDat2$Diet.Inv)
#130

### convert traits of interest to numeric
DietTraitDat2$MUTIscore <- as.numeric(DietTraitDat2$MUTIscore)
DietTraitDat2$Mass_log <- as.numeric(DietTraitDat2$Mass_log)
DietTraitDat2$Diet.Inv <- as.numeric(DietTraitDat2$Diet.Inv)

# Run phylogenetic linear model for MUTI
MUTI_GLS_invert <- gls(MUTIscore~ Diet.Inv + Mass_log, data = DietTraitDat2, 
                      correlation = corPagel(0.5, phy=Dietphy2,fixed=F, form = ~Species_Jetz), 
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
UNDataUT <- C_Diet_dat2 %>% filter(!is.na(Urban)) 
DietData3 <- UNDataUT %>% filter(!is.na(Diet.Inv)) 
length(DietData3$Diet.Inv)
#129 species with UAI and CT

###### add and pair tree

DietData3 <- as.data.frame(DietData3)
# add rownames to data
row.names(DietData3) <- DietData3$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Dietphydat3 <- geiger::treedata(tree_out,DietData3, sort=T)

Dietphy3 <- Dietphydat3$phy
DietTraitDat3 <- as.data.frame(Dietphydat3$data)

str(DietTraitDat3)
length(DietTraitDat3$Diet.Inv)
#129

### convert traits of interest to numeric
DietTraitDat3$Urban <- as.numeric(DietTraitDat3$Urban)
DietTraitDat3$Mass_log <- as.numeric(DietTraitDat3$Mass_log)
DietTraitDat3$Diet.Inv <- as.numeric(DietTraitDat3$Diet.Inv)


# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
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


# try to fix alpha at exp(4)/t (at the default upper bounds)
# this fails to converge
phyglm_UN_Invert_fix_4 <- phyloglm( Urban ~ scale(Diet.Inv) + scale(Mass_log), 
                                    data = DietTraitDat3, 
                                    phy = Dietphy3, 
                                    start.alpha = 0.55,
                                    log.alpha.bound = 4)

# try increasing upper bounds a small amount over 4
# both these models converge
set.seed(568)
phyglm_UN_Invert_fix_4.05 <- phyloglm( Urban ~ scale(Diet.Inv) + scale(Mass_log), 
                                       data = DietTraitDat3, 
                                       phy = Dietphy3, 
                                       start.alpha = 0.55,
                                       log.alpha.bound = 4.05, boot=1000)

summary(phyglm_UN_Invert_fix_4.05) # this converges

phyglm_UN_Invert_fix_4.1 <- phyloglm( Urban ~ scale(Diet.Inv) + scale(Mass_log), 
                                      data = DietTraitDat3, 
                                      phy = Dietphy3, 
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
                         data = DietTraitDat3)
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
UAIDataUT <- C_Diet_dat2 %>% filter(!is.na(aveUAI)) 
DietData4 <- UAIDataUT %>% filter(!is.na(Diet.Vert)) 
length(DietData4$Diet.Vert)
#798 species with UAI and CT

###### add and pair tree

DietData4 <- as.data.frame(DietData4)
# add rownames to data
row.names(DietData4) <- DietData4$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Dietphydat4 <- geiger::treedata(tree_out,DietData4, sort=T)

Dietphy4 <- Dietphydat4$phy
DietTraitDat4 <- as.data.frame(Dietphydat4$data)

str(DietTraitDat4)
length(DietTraitDat4$Diet.Vert)
#798

### convert traits of interest to numeric
DietTraitDat4$aveUAI <- as.numeric(DietTraitDat4$aveUAI)
DietTraitDat4$Mass_log <- as.numeric(DietTraitDat4$Mass_log)
DietTraitDat4$Diet.Vert <- as.numeric(DietTraitDat4$Diet.Vert)

# Run phylogenetic linear model for UAI
UAI_GLS_vert <- gls(aveUAI~ Diet.Vert + Mass_log, data = DietTraitDat4, 
                      correlation = corPagel(0.5, phy=Dietphy4,fixed=F, form = ~Species_Jetz), 
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
MUTIDataUT <- C_Diet_dat2 %>% filter(!is.na(MUTIscore)) 
DietData5 <- MUTIDataUT %>% filter(!is.na(Diet.Vert)) 
length(DietData5$Diet.Vert)
#130 species with UAI and CT

###### add and pair tree

DietData5 <- as.data.frame(DietData5)
# add rownames to data
row.names(DietData5) <- DietData5$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Dietphydat5 <- geiger::treedata(tree_out,DietData5, sort=T)

Dietphy5 <- Dietphydat5$phy
DietTraitDat5 <- as.data.frame(Dietphydat5$data)

str(DietTraitDat5)
length(DietTraitDat5$Diet.Vert)
#130

### convert traits of interest to numeric
DietTraitDat5$MUTIscore <- as.numeric(DietTraitDat5$MUTIscore)
DietTraitDat5$Mass_log <- as.numeric(DietTraitDat5$Mass_log)
DietTraitDat5$Diet.Vert <- as.numeric(DietTraitDat5$Diet.Vert)


# Run phylogenetic linear model for MUTI
MUTI_GLS_vert <- gls(MUTIscore~ Diet.Vert + Mass_log, data = DietTraitDat5, 
                    correlation = corPagel(0.5, phy=Dietphy5,fixed=F, form = ~Species_Jetz), 
                    method = "ML") 

# model summary and results
summary(MUTI_GLS_vert)
confint(MUTI_GLS_vert)
confint(MUTI_GLS_vert, level = 0.85)

# model diagnostics
check_model(MUTI_GLS_vert)  
qqnorm(resid(MUTI_GLS_vert)) 
qqline(resid(MUTI_GLS_vert)) 
hist(resid(MUTI_GLS_vert)) 

# save model for easy retrieval 
saveRDS(MUTI_GLS_vert, here("Models/MUTI", "MUTI_GLS_vert.rds"))


######################## UN and % Diet Vertebrates ##########################


# create a new data frame by removing species with no UN value or that are missing % diet vert
UNDataUT <- C_Diet_dat2 %>% filter(!is.na(Urban)) 
DietData6 <- UNDataUT %>% filter(!is.na(Diet.Vert)) 
length(DietData6$Diet.Vert)
#129 species with UAI and CT

###### add and pair tree

DietData6 <- as.data.frame(DietData6)
# add rownames to data
row.names(DietData6) <- DietData6$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Dietphydat6 <- geiger::treedata(tree_out,DietData6, sort=T)

Dietphy6 <- Dietphydat6$phy
DietTraitDat6 <- as.data.frame(Dietphydat6$data)

str(DietTraitDat6)
length(DietTraitDat6$Diet.Vert)
#129

### convert traits of interest to numeric
DietTraitDat6$Urban <- as.numeric(DietTraitDat6$Urban)
DietTraitDat6$Mass_log <- as.numeric(DietTraitDat6$Mass_log)
DietTraitDat6$Diet.Vert <- as.numeric(DietTraitDat6$Diet.Vert)


# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
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

# try to fix alpha at upper bounds
phyglm_UN_Vert_fix <- phyloglm( Urban ~ scale(Diet.Vert) + scale(Mass_log), 
                                data = DietTraitDat6, 
                                phy = Dietphy6, 
                                start.alpha = 0.55,
                                log.alpha.bound = 4, boot = 1000)
summary(phyglm_UN_Vert_fix) # this fails to converge



# try increasing upper bounds a small amount over 4
# both these models converge
set.seed(351)
phyglm_UN_Vert_fix_4.05 <- phyloglm( Urban ~ scale(Diet.Vert) + scale(Mass_log), 
                                       data = DietTraitDat6, 
                                       phy = Dietphy6, 
                                       start.alpha = 0.55,
                                       log.alpha.bound = 4.05, boot=1000)

summary(phyglm_UN_Vert_fix_4.05) # this converges

phyglm_UN_Vert_fix_4.1 <- phyloglm( Urban ~ scale(Diet.Vert) + scale(Mass_log), 
                                      data = DietTraitDat6, 
                                      phy = Dietphy6, 
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
                       data = DietTraitDat6)
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
UAIDataUT <- C_Diet_dat2 %>% filter(!is.na(aveUAI)) 
DietData7 <- UAIDataUT %>% filter(!is.na(Diet.PS)) 
length(DietData7$Diet.PS)
#798 species with UAI and CT

###### add and pair tree

DietData7 <- as.data.frame(DietData7)
# add rownames to data
row.names(DietData7) <- DietData7$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Dietphydat7 <- geiger::treedata(tree_out,DietData7, sort=T)

Dietphy7 <- Dietphydat7$phy
DietTraitDat7 <- as.data.frame(Dietphydat7$data)

str(DietTraitDat7)
length(DietTraitDat7$Diet.PS)
#798

### convert traits of interest to numeric
DietTraitDat7$aveUAI <- as.numeric(DietTraitDat7$aveUAI)
DietTraitDat7$Mass_log <- as.numeric(DietTraitDat7$Mass_log)
DietTraitDat7$Diet.PS <- as.numeric(DietTraitDat7$Diet.PS)

# run phylogenetic linear model for UAI using gls()
UAI_GLS_PS <- gls(aveUAI ~ Diet.PS + Mass_log, data = DietTraitDat7, 
                     correlation = corPagel(0.5, phy=Dietphy7,fixed=F, form = ~Species_Jetz), 
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
MUTIDataUT <- C_Diet_dat2 %>% filter(!is.na(MUTIscore)) 
DietData8 <- MUTIDataUT %>% filter(!is.na(Diet.PS)) 
length(DietData8$Diet.PS)
#798 species with UAI and CT

###### add and pair tree

DietData8 <- as.data.frame(DietData8)
# add rownames to data
row.names(DietData8) <- DietData8$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Dietphydat8 <- geiger::treedata(tree_out,DietData8, sort=T)

Dietphy8 <- Dietphydat8$phy
DietTraitDat8 <- as.data.frame(Dietphydat8$data)

str(DietTraitDat8)
length(DietTraitDat8$Diet.PS)
#130

### convert traits of interest to numeric
DietTraitDat8$MUTIscore <- as.numeric(DietTraitDat8$MUTIscore)
DietTraitDat8$Mass_log <- as.numeric(DietTraitDat8$Mass_log)
DietTraitDat8$Diet.PS <- as.numeric(DietTraitDat8$Diet.PS)


# run phylogenetic linear model for MUTI using gls()
MUTI_GLS_PS <- gls(MUTIscore ~ Diet.PS + Mass_log, data = DietTraitDat8, 
                  correlation = corPagel(0.5, phy=Dietphy8,fixed=F, form = ~Species_Jetz), 
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
UNDataUT <- C_Diet_dat2 %>% filter(!is.na(Urban)) 
DietData9 <- UNDataUT %>% filter(!is.na(Diet.PS)) 
length(DietData9$Diet.PS)
#129 species with UAI and CT

###### add and pair tree

DietData9 <- as.data.frame(DietData9)
# add rownames to data
row.names(DietData9) <- DietData9$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Dietphydat9 <- geiger::treedata(tree_out,DietData9, sort=T)

Dietphy9 <- Dietphydat9$phy
DietTraitDat9 <- as.data.frame(Dietphydat9$data)

str(DietTraitDat9)
length(DietTraitDat9$Diet.PS)
#129

### convert traits of interest to numeric
DietTraitDat9$Urban <- as.numeric(DietTraitDat9$Urban)
DietTraitDat9$Mass_log <- as.numeric(DietTraitDat9$Mass_log)
DietTraitDat9$Diet.PS <- as.numeric(DietTraitDat9$Diet.PS)

# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
set.seed(382)
phyglm_UN_PS_scale <- phyloglm( Urban ~ scale(Diet.PS) + scale(Mass_log), 
                                data = DietTraitDat9, 
                                phy = Dietphy9, 
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
UAIDataUT <- C_Diet_dat2 %>% filter(!is.na(aveUAI)) 
DietData10 <- UAIDataUT %>% filter(!is.na(Diet.FN)) 
length(DietData10$Diet.FN)
#798 species with UAI and CT

###### add and pair tree
DietData10 <- as.data.frame(DietData10)
# add rownames to data
row.names(DietData10) <- DietData10$Species_Jetz

# import tree
tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Dietphydat10 <- geiger::treedata(tree_out,DietData10, sort=T)

Dietphy10 <- Dietphydat10$phy
DietTraitDat10 <- as.data.frame(Dietphydat10$data)

str(DietTraitDat10)
length(DietTraitDat10$Diet.FN)
#798

### convert traits of interest to numeric
DietTraitDat10$aveUAI <- as.numeric(DietTraitDat10$aveUAI)
DietTraitDat10$Mass_log <- as.numeric(DietTraitDat10$Mass_log)
DietTraitDat10$Diet.FN <- as.numeric(DietTraitDat10$Diet.FN)


# run phylogenetic linear model for UAI using gls()
UAI_GLS_FN <- gls(aveUAI ~ Diet.FN + Mass_log, data = DietTraitDat10, 
                   correlation = corPagel(0.5, phy=Dietphy10,fixed=F, form = ~Species_Jetz), 
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
MUTIDataUT <- C_Diet_dat2 %>% filter(!is.na(MUTIscore)) 
DietData11 <- MUTIDataUT %>% filter(!is.na(Diet.FN)) 
length(DietData11$Diet.FN)
#130 species with UAI and CT

###### add and pair tree
DietData11 <- as.data.frame(DietData11)
# add rownames to data
row.names(DietData11) <- DietData11$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Dietphydat11 <- geiger::treedata(tree_out,DietData11, sort=T)

Dietphy11 <- Dietphydat11$phy
DietTraitDat11 <- as.data.frame(Dietphydat11$data)

str(DietTraitDat11)
length(DietTraitDat11$Diet.FN)
#798

### convert traits of interest to numeric
DietTraitDat11$MUTIscore <- as.numeric(DietTraitDat11$MUTIscore)
DietTraitDat11$Mass_log <- as.numeric(DietTraitDat11$Mass_log)
DietTraitDat11$Diet.FN <- as.numeric(DietTraitDat11$Diet.FN)


# run phylogenetic linear model for MUTI using gls()
MUTI_GLS_FN <- gls(MUTIscore ~ Diet.FN + Mass_log, data = DietTraitDat11, 
                  correlation = corPagel(0.5, phy=Dietphy11,fixed=F, form = ~Species_Jetz), 
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
UNDataUT <- C_Diet_dat2 %>% filter(!is.na(Urban)) 
DietData12 <- UNDataUT %>% filter(!is.na(Diet.FN)) 
length(DietData12$Diet.FN)
#129 species with UAI and CT

###### add and pair tree
DietData12 <- as.data.frame(DietData12)
# add rownames to data
row.names(DietData12) <- DietData12$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Dietphydat12 <- geiger::treedata(tree_out,DietData12, sort=T)

Dietphy12 <- Dietphydat12$phy
DietTraitDat12 <- as.data.frame(Dietphydat12$data)

str(DietTraitDat12)
length(DietTraitDat12$Diet.FN)
#129

### convert traits of interest to numeric
DietTraitDat12$Urban <- as.numeric(DietTraitDat12$Urban)
DietTraitDat12$Mass_log <- as.numeric(DietTraitDat12$Mass_log)
DietTraitDat12$Diet.FN <- as.numeric(DietTraitDat12$Diet.FN)

# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
set.seed(380)
phyglm_UN_FN_scale <- phyloglm( Urban ~ scale(Diet.FN) + scale(Mass_log), 
                                data = DietTraitDat12, 
                                phy = Dietphy12, 
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
