##########################################################
######### Coastal Birds Urban Tolerance Project ##########
##########################################################
# Step 8 Phylogenetic Trait Models - LIFE HISTORY TRAITS
# Authors: Sarah L. Jennings, Emma M. Garrison
##########################################################
# The objective of this script is to run phylogenetic models using life history traits
# 4 life history traits: brood value, clutch size, longevity and developmental mode


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

# load "Coastal_Species_LifeHistory.rds"

C_LifeHist_dat <- readRDS(here("Outputs", "Coastal_Species_LifeHistory.rds"))
str(C_LifeHist_dat)

C_LifeHist_dat2 <- C_LifeHist_dat %>%
  mutate(Species_Jetz  = str_replace(Species_Jetz, " ", "_"))
str(C_LifeHist_dat2)

C_LifeHist_dat2$Urban <- ifelse(C_LifeHist_dat2$Urban == "U", 1, 0)
View(C_LifeHist_dat2)
colnames(C_LifeHist_dat2)


######################## UAI and Brood Value ##########################

# create a new data frame that contains only species with both UAI and brood values
UAI_Brood <- C_LifeHist_dat2 %>% filter(!is.na(aveUAI)) %>% 
  filter(!is.na(brood_value)) %>% as.data.frame()
length(UAI_Brood$brood_value)
#480 species with UAI and brood_value

###### add and pair tree

# add rownames to data
row.names(UAI_Brood) <- UAI_Brood$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_Brood_phydat <- treedata(tree_out, UAI_Brood, sort=T)

UAI_Brood_phy <- UAI_Brood_phydat$phy
UAI_Brood_dat <- as.data.frame(UAI_Brood_phydat$data)

str(UAI_Brood_dat)
length(UAI_Brood_dat$brood_value)
#480

### convert traits of interest to numeric
UAI_Brood_dat$aveUAI <- as.numeric(UAI_Brood_dat$aveUAI)
UAI_Brood_dat$Mass_log <- as.numeric(UAI_Brood_dat$Mass_log)
UAI_Brood_dat$brood_value <- as.numeric(UAI_Brood_dat$brood_value)

# Run phylogenetic linear model
UAI_GLS_bv <- gls(aveUAI ~ brood_value + Mass_log, data = UAI_Brood_dat, 
                      correlation = corPagel(0.5, phy = UAI_Brood_phy, fixed=F, form = ~Species_Jetz), 
                      method = "ML") 

# model summary and results
summary(UAI_GLS_bv) 
confint(UAI_GLS_bv)

# model diagnostics
check_model(UAI_GLS_bv) 
qqnorm(resid(UAI_GLS_bv)) 
qqline(resid(UAI_GLS_bv)) 
hist(resid(UAI_GLS_bv)) 

# save model for easy retrieval 
saveRDS(UAI_GLS_bv, here("Models/UAI", "UAI_GLS_bv.rds"))

# There is one species with an extreme brood value
# Re-run model with UAI and Brood Value without this species to see if this relationship is still significantly negative 

# Filter out brood_value that are less than -5
UAI_Brood_Filtered <- UAI_Brood %>% filter(brood_value >= -5)

# Add rownames to filtered data
row.names(UAI_Brood_Filtered) <- UAI_Brood_Filtered$Species_Jetz

# Pair with the tree
UAI_Brood_Fil_phydat <- treedata(tree_out, UAI_Brood_Filtered, sort=T)

UAI_Brood_Fil_phy <- UAI_Brood_Fil_phydat$phy
UAI_Brood_Fil_dat <- as.data.frame(UAI_Brood_Fil_phydat$data)

# Convert traits of interest to numeric
UAI_Brood_Fil_dat$aveUAI <- as.numeric(UAI_Brood_Fil_dat$aveUAI)
UAI_Brood_Fil_dat$Mass_log <- as.numeric(UAI_Brood_Fil_dat$Mass_log)
UAI_Brood_Fil_dat$brood_value <- as.numeric(UAI_Brood_Fil_dat$brood_value)

# Run phylogenetic linear model with extreme brood value removed
UAI_GLS_bv_filtered <- gls(aveUAI ~ brood_value + Mass_log, data = UAI_Brood_Fil_dat, 
                           correlation = corPagel(0.5, phy = UAI_Brood_Fil_phy, fixed = F, form = ~Species_Jetz), 
                           method = "ML")

# model summary and results
summary(UAI_GLS_bv_filtered) # still significant (p = 0.0288)
confint(UAI_GLS_bv_filtered) # still significant (C.I. = (-0.349, -0.0195)

# model diagnostics
check_model(UAI_GLS_bv_filtered)
qqnorm(resid(UAI_GLS_bv_filtered)) 
qqline(resid(UAI_GLS_bv_filtered))
hist(resid(UAI_GLS_bv_filtered))

# save model for easy retrieval 
saveRDS(UAI_GLS_bv_filtered, here("Models/UAI", "UAI_GLS_bv_filtered.rds"))


######################## MUTI and Brood Value ##########################

# create a new data frame that contains only species with both MUTI and brood values
MUTI_Brood <- C_LifeHist_dat2 %>% filter(!is.na(MUTIscore)) %>% 
  filter(!is.na(brood_value)) %>% as.data.frame()
length(MUTI_Brood$brood_value)
#121 species with MUTI and brood_value

###### add and pair tree

# add rownames to data
row.names(MUTI_Brood) <- MUTI_Brood$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_Brood_phydat <- treedata(tree_out, MUTI_Brood, sort=T)

MUTI_Brood_phy <- MUTI_Brood_phydat$phy
MUTI_Brood_dat <- as.data.frame(MUTI_Brood_phydat$data)

str(MUTI_Brood_dat)
length(MUTI_Brood_dat$brood_value)
#121

### convert traits of interest to numeric
MUTI_Brood_dat$MUTIscore <- as.numeric(MUTI_Brood_dat$MUTIscore)
MUTI_Brood_dat$Mass_log <- as.numeric(MUTI_Brood_dat$Mass_log)
MUTI_Brood_dat$brood_value <- as.numeric(MUTI_Brood_dat$brood_value)


# Run phylogenetic linear model
MUTI_GLS_bv <- gls(MUTIscore~ brood_value + Mass_log, data = MUTI_Brood_dat, 
                  correlation = corPagel(0.5, phy = MUTI_Brood_phy, fixed=F, form = ~Species_Jetz), 
                  method = "ML") 

# model summary and results
summary(MUTI_GLS_bv) 
confint(MUTI_GLS_bv)

# model diagnostics
check_model(MUTI_GLS_bv) 
qqnorm(resid(MUTI_GLS_bv)) 
qqline(resid(MUTI_GLS_bv)) 
hist(resid(MUTI_GLS_bv)) 

# save model for easy retrieval 
saveRDS(MUTI_GLS_bv, here("Models/MUTI", "MUTI_GLS_bv.rds"))

######################## UN and Brood Value ##########################

# create a new data frame that contains only species with both UN and brood values
UN_Brood <- C_LifeHist_dat2 %>% filter(!is.na(Urban)) %>% 
  filter(!is.na(brood_value)) %>%
  column_to_rownames(., var="Species_Jetz")
length(UN_Brood$brood_value)
#102 species with UN and brood_value

###### add and pair tree

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UN_Brood_phydat <- treedata(tree_out, UN_Brood, sort=T)

UN_Brood_phy <- UN_Brood_phydat$phy
UN_Brood_dat <- as.data.frame(UN_Brood_phydat$data)

str(UN_Brood_dat)
length(UN_Brood_dat$brood_value)
#102

### convert traits of interest to numeric
UN_Brood_dat$Urban <- as.numeric(UN_Brood_dat$Urban)
UN_Brood_dat$Mass_log <- as.numeric(UN_Brood_dat$Mass_log)
UN_Brood_dat$brood_value <- as.numeric(UN_Brood_dat$brood_value)

# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
set.seed(305)
phyglm_UN_bv_scale <- phyloglm( Urban ~ scale(brood_value) + scale(Mass_log), 
                                data = UN_Brood_dat, 
                                phy = UN_Brood_phy, 
                                boot = 1000) 
summary(phyglm_UN_bv_scale) 
# this successfully converges 
# alpha is at upper bounds

# save model
saveRDS(phyglm_UN_bv_scale, here("Models/UN", "phyglm_UN_bv_scale.rds"))
# load model
phyglm_UN_bv_scale <- readRDS(here("Models/UN", "phyglm_UN_bv_scale.rds"))


# run a non-phylogenetic logistic model for comparison
glm_UN_bv <- logistf(Urban ~ scale(brood_value) + scale(Mass_log), 
                     data = UN_Brood)            
summary(glm_UN_bv)
# we reach the same conclusions
# some differences in coefficients though


# get alpha, t, and half life for the model
(phyglm_UN_bv_scale$mean.tip.height) # t
(alpha_bv <- phyglm_UN_bv_scale$alpha) # alpha
(hl_bv <- log(2)/alpha_bv) # half-life
# small half life relative to t -> low phylogenetic signal



#############################################################################
#############################################################################

########################## UAI and Clutch Size ##############################


# create a new data frame that contains only species with both UAI and clutch size
UAI_Clutch <- C_LifeHist_dat2 %>% filter(!is.na(aveUAI)) %>% 
  filter(!is.na(clutch_size)) %>% as.data.frame()
length(UAI_Clutch$clutch_size)
#738 species with UAI and clutch_size

###### add and pair tree

# add rownames to data
row.names(UAI_Clutch) <- UAI_Clutch$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_Clutch_phydat <- treedata(tree_out, UAI_Clutch, sort=T)

UAI_Clutch_phy <- UAI_Clutch_phydat$phy
UAI_Clutch_dat <- as.data.frame(UAI_Clutch_phydat$data)

str(UAI_Clutch_dat)
length(UAI_Clutch_dat$clutch_size)
#738

### convert traits of interest to numeric
UAI_Clutch_dat$aveUAI <- as.numeric(UAI_Clutch_dat$aveUAI)
UAI_Clutch_dat$Mass_log <- as.numeric(UAI_Clutch_dat$Mass_log)
UAI_Clutch_dat$clutch_size <- as.numeric(UAI_Clutch_dat$clutch_size)

# Run phylogenetic linear model
UAI_GLS_clutch <- gls(aveUAI~ clutch_size + Mass_log, data = UAI_Clutch_dat, 
                   correlation = corPagel(0.5, phy = UAI_Clutch_phy, fixed=F, form = ~Species_Jetz), 
                   method = "ML") 

# model summary and results
summary(UAI_GLS_clutch) 
confint(UAI_GLS_clutch)

# model diagnostics
check_model(UAI_GLS_clutch) 
qqnorm(resid(UAI_GLS_clutch)) 
qqline(resid(UAI_GLS_clutch))
hist(resid(UAI_GLS_clutch)) 

# save model for easy retrieval 
saveRDS(UAI_GLS_clutch, here("Models/UAI", "UAI_GLS_clutch.rds"))


######################## MUTI and Clutch Size ##########################


# create a new data frame that contains only species with both MUTI and clutch size
MUTI_Clutch <- C_LifeHist_dat2 %>% filter(!is.na(MUTIscore)) %>% 
  filter(!is.na(clutch_size)) %>% as.data.frame()
length(MUTI_Clutch$clutch_size)
# 126 species with MUTI and clutch_size

###### add and pair tree

# add rownames to data
row.names(MUTI_Clutch) <- MUTI_Clutch$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_Clutch_phydat <- treedata(tree_out, MUTI_Clutch, sort=T)

MUTI_Clutch_phy <- MUTI_Clutch_phydat$phy
MUTI_Clutch_dat <- as.data.frame(MUTI_Clutch_phydat$data)

str(MUTI_Clutch_dat)
length(MUTI_Clutch_dat$clutch_size)
# 126

### convert traits of interest to numeric
MUTI_Clutch_dat$MUTIscore <- as.numeric(MUTI_Clutch_dat$MUTIscore)
MUTI_Clutch_dat$Mass_log <- as.numeric(MUTI_Clutch_dat$Mass_log)
MUTI_Clutch_dat$clutch_size <- as.numeric(MUTI_Clutch_dat$clutch_size)

# Run phylogenetic linear model
MUTI_GLS_clutch <- gls(MUTIscore~ clutch_size + Mass_log, data = MUTI_Clutch_dat, 
                      correlation = corPagel(0.5, phy = MUTI_Clutch_phy, fixed=F, form = ~Species_Jetz), 
                      method = "ML") 

# model summary and results
summary(MUTI_GLS_clutch) 
confint(MUTI_GLS_clutch)

# model diagnostics
check_model(MUTI_GLS_clutch) 
qqnorm(resid(MUTI_GLS_clutch)) 
qqline(resid(MUTI_GLS_clutch))
hist(resid(MUTI_GLS_clutch)) 

# save model for easy retrieval 
saveRDS(MUTI_GLS_clutch, here("Models/MUTI", "MUTI_GLS_clutch.rds"))



######################## UN and Clutch Size ##########################

# create a new data frame that contains only species with both UN and clutch size
UN_Clutch <- C_LifeHist_dat2 %>% filter(!is.na(Urban)) %>% 
  filter(!is.na(clutch_size)) %>% column_to_rownames(., var = "Species_Jetz")
length(UN_Clutch$clutch_size)
#122 species with UN and clutch_size

###### add and pair tree

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UN_Clutch_phydat <- treedata(tree_out, UN_Clutch, sort=T)

UN_Clutch_phy <- UN_Clutch_phydat$phy
UN_Clutch_dat <- as.data.frame(UN_Clutch_phydat$data)

str(UN_Clutch_dat)
length(UN_Clutch_dat$clutch_size)
#122

### convert traits of interest to numeric
UN_Clutch_dat$Urban <- as.numeric(UN_Clutch_dat$Urban)
UN_Clutch_dat$Mass_log <- as.numeric(UN_Clutch_dat$Mass_log)
UN_Clutch_dat$clutch_size <- as.numeric(UN_Clutch_dat$clutch_size)

# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
phyglm_UN_clutch_scale <- phyloglm( Urban ~ scale(clutch_size) + scale(Mass_log), 
                                    data = UN_Clutch_dat, 
                                    phy = UN_Clutch_phy, 
                                    boot = 1000) 

summary(phyglm_UN_clutch_scale) 
# fails to converge


# print AIC values for models with different upper bounds
# intervals of 0.1 from 0 up to 4
for (i in seq(0, 4, by = 0.1)) {
  print(phyloglm(Urban ~ scale(clutch_size) + scale(Mass_log), 
                 data = UN_Clutch_dat, 
                 phy = UN_Clutch_phy,
                 log.alpha.bound = i)$aic)
}
# AIC values support models with larger values of alpha (low phylo signal)


# fix alpha at upper bounds
set.seed(770)
phyglm_UN_clutch_fix <- phyloglm(Urban ~ scale(clutch_size) + scale(Mass_log), 
                                  data = UN_Clutch_dat, 
                                  phy = UN_Clutch_phy,
                                  log.alpha.bound = 4,
                                  start.alpha = 0.53,
                                  boot = 1000) 

summary(phyglm_UN_clutch_fix)


# save model
saveRDS(phyglm_UN_clutch_fix, here("Models/UN", "phyglm_UN_clutch_fix.rds"))
# load model
phyglm_UN_clutch_fix <- readRDS(here("Models/UN", "phyglm_UN_clutch_fix.rds"))


# run non-phylogenetic glm for comparison
library(logistf)
glm_UN_clutch <- logistf(Urban ~ scale(clutch_size) + scale(Mass_log), 
                         data = UN_Clutch)            
summary(glm_UN_clutch)
# we reach the same conclusions
# fairly similar coefficients to fixed model above


# get alpha, t, and half life for the model
(phyglm_UN_clutch_fix$mean.tip.height) # t
(alpha_clutch <- phyglm_UN_clutch_fix$alpha) # alpha
(hl_clutch <- log(2)/alpha_clutch) # half life
# small half life relative to t -> low phylogenetic signal

#############################################################################
#############################################################################

######################## UAI and Longevity ##########################


# create a new data frame that contains only species with both UAI and longevity
UAI_Long <- C_LifeHist_dat2 %>% filter(!is.na(aveUAI)) %>% 
  filter(!is.na(longevity)) %>% as.data.frame()
length(UAI_Long$longevity)
# 796 species with UAI and longevity

###### add and pair tree

# add rownames to data
row.names(UAI_Long) <- UAI_Long$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_Long_phydat <- treedata(tree_out, UAI_Long, sort=T)

UAI_Long_phy <- UAI_Long_phydat$phy
UAI_Long_dat <- as.data.frame(UAI_Long_phydat$data)

str(UAI_Long_dat)
length(UAI_Long_dat$longevity)
# 796

### convert traits of interest to numeric
UAI_Long_dat$aveUAI <- as.numeric(UAI_Long_dat$aveUAI)
UAI_Long_dat$Mass_log <- as.numeric(UAI_Long_dat$Mass_log)
UAI_Long_dat$longevity <- as.numeric(UAI_Long_dat$longevity)

# Run phylogenetic linear model
UAI_GLS_long <- gls(aveUAI~ longevity + Mass_log, data = UAI_Long_dat, 
                     correlation = corPagel(0.5, phy = UAI_Long_phy, fixed=F, form = ~Species_Jetz), 
                     method = "ML") 

# model summary and results
summary(UAI_GLS_long) 
confint(UAI_GLS_long)

# model diagnostics
check_model(UAI_GLS_long) 
qqnorm(resid(UAI_GLS_long)) 
qqline(resid(UAI_GLS_long))
hist(resid(UAI_GLS_long)) 

# save model for easy retrieval 
saveRDS(UAI_GLS_long, here("Models/UAI", "UAI_GLS_long.rds"))

######################## MUTI and Longevity ##########################


# create a new data frame that contains only species with both MUTI and longevity
MUTI_Long <- C_LifeHist_dat2 %>% filter(!is.na(MUTIscore)) %>% 
  filter(!is.na(longevity)) %>% as.data.frame()
length(MUTI_Long$longevity)
# 130 species with MUTI and longevity

###### add and pair tree

# add rownames to data
row.names(MUTI_Long) <- MUTI_Long$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_Long_phydat <- treedata(tree_out, MUTI_Long, sort=T)

MUTI_Long_phy <- MUTI_Long_phydat$phy
MUTI_Long_dat <- as.data.frame(MUTI_Long_phydat$data)

str(MUTI_Long_dat)
length(MUTI_Long_dat$longevity)
# 130

### convert traits of interest to numeric
MUTI_Long_dat$MUTIscore <- as.numeric(MUTI_Long_dat$MUTIscore)
MUTI_Long_dat$Mass_log <- as.numeric(MUTI_Long_dat$Mass_log)
MUTI_Long_dat$longevity <- as.numeric(MUTI_Long_dat$longevity)

# Run phylogenetic linear model
MUTI_GLS_long <- gls(MUTIscore~ longevity + Mass_log, data = MUTI_Long_dat, 
                       correlation = corPagel(0.5, phy = MUTI_Long_phy, fixed=F, form = ~Species_Jetz), 
                       method = "ML") 

# model summary and results
summary(MUTI_GLS_long) 
confint(MUTI_GLS_long)

# model diagnostics
check_model(MUTI_GLS_long) 
qqnorm(resid(MUTI_GLS_long)) 
qqline(resid(MUTI_GLS_long))
hist(resid(MUTI_GLS_long)) 

# save model for easy retrieval 
saveRDS(MUTI_GLS_long, here("Models/MUTI", "MUTI_GLS_long.rds"))

######################## UN and Longevity ##########################


# create a new data frame that contains only species with both UN and longevity
UN_Long <- C_LifeHist_dat2 %>% filter(!is.na(Urban)) %>% 
  filter(!is.na(longevity)) %>% column_to_rownames(., var="Species_Jetz")
length(UN_Long$longevity)
#129 species with UAI and longevity

###### add and pair tree

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UN_Long_phydat <- treedata(tree_out, UN_Long, sort=T)

UN_Long_phy <- UN_Long_phydat$phy
UN_Long_dat <- as.data.frame(UN_Long_phydat$data)

str(UN_Long_dat)
length(UN_Long_dat$longevity)
#129

### convert traits of interest to numeric
UN_Long_dat$Urban <- as.numeric(UN_Long_dat$Urban)
UN_Long_dat$Mass_log <- as.numeric(UN_Long_dat$Mass_log)
UN_Long_dat$longevity <- as.numeric(UN_Long_dat$longevity)


# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
phyglm_UN_long_scale <- phyloglm( Urban ~ scale(longevity) + scale(Mass_log), 
                                  data = UN_Long_dat, 
                                  phy = UN_Long_phy, 
                                  boot = 1000) 

# this model fails to converge
summary(phyglm_UN_long_scale)


# print AIC values for models with different upper bounds
# intervals of 0.1 from 0 up to 4
for (i in seq(0, 4, by = 0.1)) {
  print(phyloglm(Urban ~ scale(longevity) + scale(Mass_log), 
                 data = UN_Long_dat, 
                 phy = UN_Long_phy,
                 log.alpha.bound = i)$aic)
}
# AIC values support models with larger values of alpha (low phylo signal)


# try fixing alpha at upper bounds
set.seed(989)
phyglm_UN_long_fix <- phyloglm( Urban ~ scale(longevity) + scale(Mass_log), 
                                data = UN_Long_dat, 
                                phy = UN_Long_phy,
                                log.alpha.bound = 4,
                                start.alpha = 0.55,
                                boot = 1000)
summary(phyglm_UN_long_fix) # this converges


# save model
saveRDS(phyglm_UN_long_fix, here("Models/UN", "phyglm_UN_long_fix.rds"))

# load model
phyglm_UN_long_fix <- readRDS(here("Models/UN", "phyglm_UN_long_fix.rds"))


# Compare with non-phylogenetic logistic regression
glm_UN_long <- logistf(Urban ~ scale(longevity) + scale(Mass_log), 
                       data = UN_Long)
summary(glm_UN_long)
# this is similar to model with log.alpha.bound fixed at 4


# get alpha, t, and half life for the model
(phyglm_UN_long_fix$mean.tip.height) # t
(alpha_long <- phyglm_UN_long_fix$alpha) # alpha
(hl_long <- log(2)/alpha_long) # half life
# compared to t, this is a small half-life


#############################################################################
#############################################################################


######################## UAI and Developmental Mode ##########################


# create a new data frame that contains only species with both UAI and developmental mode
UAI_Develop <- C_LifeHist_dat2 %>% filter(!is.na(aveUAI)) %>% 
  filter(!is.na(developmental_mode)) %>% as.data.frame()
length(UAI_Develop$developmental_mode)
#766 species with UAI and developmental_mode

###### add and pair tree

# add rownames to data
row.names(UAI_Develop) <- UAI_Develop$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_Develop_phydat <- treedata(tree_out, UAI_Develop, sort=T)

UAI_Develop_phy <- UAI_Develop_phydat$phy
UAI_Develop_dat <- as.data.frame(UAI_Develop_phydat$data)

str(UAI_Develop_dat)
length(UAI_Develop_dat$developmental_mode)
#766

### convert traits of interest to numeric
UAI_Develop_dat$aveUAI <- as.numeric(UAI_Develop_dat$aveUAI)
UAI_Develop_dat$Mass_log <- as.numeric(UAI_Develop_dat$Mass_log)
UAI_Develop_dat$developmental_mode <- as.numeric(UAI_Develop_dat$developmental_mode)


# Run phylogenetic linear model
UAI_GLS_develop <- gls(aveUAI~ developmental_mode + Mass_log, data = UAI_Develop_dat, 
                     correlation = corPagel(0.5, phy=UAI_Develop_phy, fixed=F, form = ~Species_Jetz), 
                     method = "ML") 

# model summary and results
summary(UAI_GLS_develop) 
confint(UAI_GLS_develop)

# model diagnostics
check_model(UAI_GLS_develop) 
qqnorm(resid(UAI_GLS_develop)) 
qqline(resid(UAI_GLS_develop)) 
hist(resid(UAI_GLS_develop)) 


# save model for easy retrieval 
saveRDS(UAI_GLS_develop, here("Models/UAI", "UAI_GLS_develop.rds"))

######################## MUTI and Developmental Mode ##########################

MUTI_Develop <- C_LifeHist_dat2 %>% filter(!is.na(MUTIscore)) %>% 
  filter(!is.na(developmental_mode)) %>% as.data.frame()
length(MUTI_Develop$developmental_mode)
# 127 species with MUTI and developmental_mode

###### add and pair tree

# add rownames to data
row.names(MUTI_Develop) <- MUTI_Develop$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_Develop_phydat <- treedata(tree_out, MUTI_Develop, sort=T)

MUTI_Develop_phy <- MUTI_Develop_phydat$phy
MUTI_Develop_dat <- as.data.frame(MUTI_Develop_phydat$data)

str(MUTI_Develop_dat)
length(MUTI_Develop_dat$developmental_mode)
# 127

### convert traits of interest to numeric
MUTI_Develop_dat$MUTIscore <- as.numeric(MUTI_Develop_dat$MUTIscore)
MUTI_Develop_dat$Mass_log <- as.numeric(MUTI_Develop_dat$Mass_log)
MUTI_Develop_dat$developmental_mode <- as.numeric(MUTI_Develop_dat$developmental_mode)


# Run phylogenetic linear model
MUTI_GLS_develop <- gls(MUTIscore~ developmental_mode + Mass_log, data = MUTI_Develop_dat, 
                       correlation = corPagel(0.5, phy=MUTI_Develop_phy, fixed=F, form = ~Species_Jetz), 
                       method = "ML") 

# model summary and results
summary(MUTI_GLS_develop) 
confint(MUTI_GLS_develop)

# model diagnostics
check_model(MUTI_GLS_develop) 
qqnorm(resid(MUTI_GLS_develop)) 
qqline(resid(MUTI_GLS_develop)) 
hist(resid(MUTI_GLS_develop)) 

# save model for easy retrieval 
saveRDS(MUTI_GLS_develop, here("Models/MUTI", "MUTI_GLS_develop.rds"))

######################## UN and Developmental Mode ##########################


# create a new data frame that contains only species with both UN and developmental mode
UN_Develop <- C_LifeHist_dat2 %>% filter(!is.na(Urban)) %>% 
  filter(!is.na(developmental_mode)) %>% column_to_rownames(., var = "Species_Jetz")
length(UN_Develop$developmental_mode)
#129 species with MUTI and developmental_mode

###### add and pair tree
tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UN_Develop_phydat <- treedata(tree_out, UN_Develop, sort=T)

UN_Develop_phy <- UN_Develop_phydat$phy
UN_Develop_dat <- as.data.frame(UN_Develop_phydat$data)

str(UN_Develop_dat)
length(UN_Develop_dat$developmental_mode)
#129

### convert traits of interest to numeric
UN_Develop_dat$Urban <- as.numeric(UN_Develop_dat$Urban)
UN_Develop_dat$Mass_log <- as.numeric(UN_Develop_dat$Mass_log)
UN_Develop_dat$developmental_mode <- as.numeric(UN_Develop_dat$developmental_mode)

# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
# not scaling developmental mode as it is 0/1 binary variable
phyglm_UN_develop_scale <- phyloglm( Urban ~ developmental_mode + scale(Mass_log), 
                                     data = UN_Develop_dat, 
                                     phy = UN_Develop_phy,
                                     boot = 1000) 
# this fails to converge
summary(phyglm_UN_develop_scale)


# print AIC values for models with different upper bounds
# intervals of 0.1 from 0 up to 4
for (i in seq(0, 4, by = 0.1)) {
  print(phyloglm(Urban ~ developmental_mode + scale(Mass_log), 
                 data = UN_Develop_dat, 
                 phy = UN_Develop_phy,
                 log.alpha.bound = i)$aic)
}
# AIC values support models with larger values of alpha (low phylo signal)


# try fixing model at upper bounds for alpha
set.seed(434)
phyglm_UN_develop_fix <- phyloglm( Urban ~ developmental_mode + scale(Mass_log), 
                                   data = UN_Develop_dat, 
                                   phy = UN_Develop_phy,
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
saveRDS(phyglm_UN_develop_fix, here("Models/UN", "phyglm_UN_develop_fix.rds"))
# load model
phyglm_UN_develop_fix <- readRDS(here("Models/UN", "phyglm_UN_develop_fix.rds"))

# compare results with a non-phylogenetic logistic model
glm_UN_develop <- logistf(Urban ~ developmental_mode + scale(Mass_log), 
                          data = UN_Develop)

summary(glm_UN_develop) # similar coefficients


# get alpha, t, and half life for the model
(phyglm_UN_develop_fix$mean.tip.height) # t
(alpha_dev <- phyglm_UN_develop_fix$alpha) # alpha
(hl_dev <- log(2)/alpha_dev) # half life
#compared to t, this is a small Half-Life

