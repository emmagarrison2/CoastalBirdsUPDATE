##########################################################
######### Coastal Birds Urban Tolerance Project ##########
##########################################################
# Step 8 Phylogenetic Trait Models - SEXUAL SELECTION TRAITS
# Authors: Sarah L. Jennings, Emma M. Garrison
##########################################################
# The objective of this script is to run phylogenetic models using sexual selection traits
# 4 sexual selection traits: brightness dichromatism, hue dichromatism, intensity on males, intensity on females


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

#load in "Coastal_Species_SSelect.rds"

C_SSelect_dat <- readRDS(here("Outputs", "Coastal_Species_SSelect.rds"))
str(C_SSelect_dat)

C_SSelect_dat2 <- C_SSelect_dat %>%
  mutate(Species_Jetz  = str_replace(Species_Jetz, " ", "_"))
str(C_SSelect_dat2)

C_SSelect_dat2$Urban <- ifelse(C_SSelect_dat2$Urban == "U", 1, 0)
View(C_SSelect_dat2)
colnames(C_SSelect_dat2)



######################## UAI and Dichromatism - BRIGHTNESS ##########################

# create a new data frame that contains only species with both UAI and brightness values
UAI_Bright <- C_SSelect_dat2 %>% filter(!is.na(aveUAI)) %>% 
  filter(!is.na(Dichrom_bright)) %>% as.data.frame()
length(UAI_Bright$Dichrom_bright)
#201 species with UAI and Dichrom_bright

###### add and pair tree

# add rownames to data
row.names(UAI_Bright) <- UAI_Bright$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_Bright_phydat <- treedata(tree_out, UAI_Bright, sort=T)

UAI_Bright_phy <- UAI_Bright_phydat$phy
UAI_Bright_dat <- as.data.frame(UAI_Bright_phydat$data)

str(UAI_Bright_dat)
length(UAI_Bright_dat$Dichrom_bright)
#201

### convert traits of interest to numeric
UAI_Bright_dat$aveUAI <- as.numeric(UAI_Bright_dat$aveUAI)
UAI_Bright_dat$Mass_log <- as.numeric(UAI_Bright_dat$Mass_log)
UAI_Bright_dat$Dichrom_bright <- as.numeric(UAI_Bright_dat$Dichrom_bright)

# Run phylogenetic linear model
UAI_GLS_bright <- gls(aveUAI~ Dichrom_bright + Mass_log, data = UAI_Bright_dat, 
                       correlation = corPagel(0.5, phy = UAI_Bright_phy, fixed=F, form = ~Species_Jetz), 
                       method = "ML") 

# model summary and results
summary(UAI_GLS_bright) 
confint(UAI_GLS_bright)

# model diagnostics
check_model(UAI_GLS_bright)
qqnorm(resid(UAI_GLS_bright)) 
qqline(resid(UAI_GLS_bright)) 
hist(resid(UAI_GLS_bright)) 

# save model for easy retrieval 
saveRDS(UAI_GLS_bright, here("Models/UAI", "UAI_GLS_bright.rds"))

######################## MUTI and Dichromatism - BRIGHTNESS ##########################

# create a new data frame that contains only species with both MUTI and brightness values
MUTI_Bright <- C_SSelect_dat2 %>% filter(!is.na(MUTIscore)) %>% 
  filter(!is.na(Dichrom_bright)) %>% as.data.frame()
length(MUTI_Bright$Dichrom_bright)
#67 species with MUTI and Dichrom_bright


###### add and pair tree

# add rownames to data
row.names(MUTI_Bright) <- MUTI_Bright$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_Bright_phydat <- treedata(tree_out, MUTI_Bright, sort=T)

MUTI_Bright_phy <- MUTI_Bright_phydat$phy
MUTI_Bright_dat <- as.data.frame(MUTI_Bright_phydat$data)

str(MUTI_Bright_dat)
length(MUTI_Bright_dat$Dichrom_bright)
#67

### convert traits of interest to numeric
MUTI_Bright_dat$MUTIscore <- as.numeric(MUTI_Bright_dat$MUTIscore)
MUTI_Bright_dat$Mass_log <- as.numeric(MUTI_Bright_dat$Mass_log)
MUTI_Bright_dat$Dichrom_bright <- as.numeric(MUTI_Bright_dat$Dichrom_bright)

# Run phylogenetic linear model
# model not converging when corPagel has starting point = 0.5 and fixed = F 
# Find a lambda value to fix in the model 

# Create an empty vector to store AIC values
AIC_values <- numeric()

# Loop through different values of the parameter for corPagel
for (i in seq(0, 1, by = 0.1)) {
  # Fit the gls model with the current value of i
  model <- gls(MUTIscore ~ Dichrom_bright + Mass_log, 
               data = MUTI_Bright_dat, 
               correlation = corPagel(i, phy = MUTI_Bright_phy, fixed = TRUE, form = ~Species_Jetz), 
               method = "ML")
  # Extract AIC value and store it in the vector
  AIC_values <- c(AIC_values, AIC(model))
}

# Print AIC values
print(AIC_values)
#0.4 = best AIC score 


# Run phylogenetic linear model with fixed lambda value
MUTI_GLS_bright <- gls(MUTIscore~ Dichrom_bright + Mass_log, data = MUTI_Bright_dat, 
                      correlation = corPagel(0.4, phy = MUTI_Bright_phy, fixed=T, form = ~Species_Jetz), 
                      method = "ML") 

# model summary and results
summary(MUTI_GLS_bright) 
confint(MUTI_GLS_bright)


# model diagnostics
check_model(MUTI_GLS_bright)
qqnorm(resid(MUTI_GLS_bright)) 
qqline(resid(MUTI_GLS_bright)) 
hist(resid(MUTI_GLS_bright)) 


# save model for easy retrieval 
saveRDS(MUTI_GLS_bright, here("Models/MUTI", "MUTI_GLS_bright.rds"))


######################## UN and Dichromatism - BRIGHTNESS ##########################

# create a new data frame that contains only species with both UN and brightness values
UN_Bright <- C_SSelect_dat2 %>% filter(!is.na(Urban)) %>% 
  filter(!is.na(Dichrom_bright)) %>% column_to_rownames(., var = "Species_Jetz")
length(UN_Bright$Dichrom_bright)
#61 species with UN and Dichrom_bright

###### add and pair tree

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UN_Bright_phydat <- treedata(tree_out, UN_Bright, sort=T)

UN_Bright_phy <- UN_Bright_phydat$phy
UN_Bright_dat <- as.data.frame(UN_Bright_phydat$data)

str(UN_Bright_dat)
length(UN_Bright_dat$Dichrom_bright)
#61


### convert traits of interest to numeric
UN_Bright_dat$Urban <- as.numeric(UN_Bright_dat$Urban)
UN_Bright_dat$Mass_log <- as.numeric(UN_Bright_dat$Mass_log)
UN_Bright_dat$Dichrom_bright <- as.numeric(UN_Bright_dat$Dichrom_bright)


# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
set.seed(624)
phyglm_UN_brightness_scale <- phyloglm( Urban ~ scale(Dichrom_bright) + scale(Mass_log), 
                                        data = UN_Bright_dat, 
                                        phy = UN_Bright_phy, 
                                        boot = 1000) 
# this model converges
summary(phyglm_UN_brightness_scale)

# save model
saveRDS(phyglm_UN_brightness_scale, here("Models/UN", "phyglm_UN_brightness_scale.rds"))
# load model
phyglm_UN_brightness_scale <- readRDS(here("Models/UN", "phyglm_UN_brightness_scale.rds"))


# get alpha, t, and half life for the model
(phyglm_UN_brightness_scale$mean.tip.height) # t
(alpha_bright <- phyglm_UN_brightness_scale$alpha) # alpha
(hl_bright <- log(2)/alpha_bright) # half life
# compared to t, this is a small/moderate half life



##############################################################################
##############################################################################

######################## UAI and Dichromatism - HUE ##########################

# create a new data frame that contains only species with both UAI and hue values
UAI_Hue <- C_SSelect_dat2 %>% filter(!is.na(aveUAI)) %>% 
  filter(!is.na(Dichrom_hue)) %>% as.data.frame()
length(UAI_Hue$Dichrom_hue)
#201 species with UAI and Dichrom_hue

###### add and pair tree

# add rownames to data
row.names(UAI_Hue) <- UAI_Hue$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_Hue_phydat <- treedata(tree_out, UAI_Hue, sort=T)

UAI_Hue_phy <- UAI_Hue_phydat$phy
UAI_Hue_dat <- as.data.frame(UAI_Hue_phydat$data)

str(UAI_Hue_dat)
length(UAI_Hue_dat$Dichrom_hue)
#201


### convert traits of interest to numeric
UAI_Hue_dat$aveUAI <- as.numeric(UAI_Hue_dat$aveUAI)
UAI_Hue_dat$Mass_log <- as.numeric(UAI_Hue_dat$Mass_log)
UAI_Hue_dat$Dichrom_hue <- as.numeric(UAI_Hue_dat$Dichrom_hue)

# Run phylogenetic linear model
UAI_GLS_hue <- gls(aveUAI~ Dichrom_hue + Mass_log, data = UAI_Hue_dat, 
                      correlation = corPagel(0.5, phy= UAI_Hue_phy, fixed=F, form = ~Species_Jetz), 
                      method = "ML") 

# model summary and results
summary(UAI_GLS_hue) 
confint(UAI_GLS_hue)

# model diagnostics
check_model(UAI_GLS_hue) 
qqnorm(resid(UAI_GLS_hue)) 
qqline(resid(UAI_GLS_hue)) 
hist(resid(UAI_GLS_hue)) 


# save model for easy retrieval 
saveRDS(UAI_GLS_hue, here("Models/UAI", "UAI_GLS_hue.rds"))


######################## MUTI and Dichromatism - HUE ##########################

# create a new data frame that contains only species with both MUTI and hue values
MUTI_Hue <- C_SSelect_dat2 %>% filter(!is.na(MUTIscore)) %>% 
  filter(!is.na(Dichrom_hue)) %>% as.data.frame()
length(MUTI_Hue$Dichrom_hue)
#67 species with MUTI and Dichrom_hue

###### add and pair tree

# add rownames to data
row.names(MUTI_Hue) <- MUTI_Hue$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_Hue_phydat <- treedata(tree_out, MUTI_Hue, sort=T)

MUTI_Hue_phy <- MUTI_Hue_phydat$phy
MUTI_Hue_dat <- as.data.frame(MUTI_Hue_phydat$data)

str(MUTI_Hue_dat)
length(MUTI_Hue_dat$Dichrom_hue)
#67

### convert traits of interest to numeric
MUTI_Hue_dat$MUTIscore <- as.numeric(MUTI_Hue_dat$MUTIscore)
MUTI_Hue_dat$Mass_log <- as.numeric(MUTI_Hue_dat$Mass_log)
MUTI_Hue_dat$Dichrom_hue <- as.numeric(MUTI_Hue_dat$Dichrom_hue)

# Run phylogenetic linear model
MUTI_GLS_hue <- gls(MUTIscore~ Dichrom_hue + Mass_log, data = MUTI_Hue_dat, 
                   correlation = corPagel(0.5, phy = MUTI_Hue_phy, fixed = F, form = ~Species_Jetz), 
                   method = "ML") 

# model summary and results
summary(MUTI_GLS_hue)
confint(MUTI_GLS_hue)

# model diagnostics
check_model(MUTI_GLS_hue) 
qqnorm(resid(MUTI_GLS_hue)) 
qqline(resid(MUTI_GLS_hue)) 
hist(resid(MUTI_GLS_hue)) 


# save model for easy retrieval 
saveRDS(MUTI_GLS_hue, here("Models/MUTI", "MUTI_GLS_hue.rds"))


######################## UN and Dichromatism - HUE ##########################

# create a new data frame that contains only species with both UN and hue values
UN_Hue <- C_SSelect_dat2 %>% filter(!is.na(Urban)) %>% 
  filter(!is.na(Dichrom_hue)) %>% column_to_rownames(., var = "Species_Jetz")
length(UN_Hue$Dichrom_hue)
#61 species with UN and Dichrom_hue

###### add and pair tree

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UN_Hue_phydat <- treedata(tree_out, UN_Hue, sort=T)

UN_Hue_phy <- UN_Hue_phydat$phy
UN_Hue_dat <- as.data.frame(UN_Hue_phydat$data)

str(UN_Hue_dat)
length(UN_Hue_dat$Dichrom_hue)
#61

### convert traits of interest to numeric
UN_Hue_dat$Urban <- as.numeric(UN_Hue_dat$Urban)
UN_Hue_dat$Mass_log <- as.numeric(UN_Hue_dat$Mass_log)
UN_Hue_dat$Dichrom_hue <- as.numeric(UN_Hue_dat$Dichrom_hue)


# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
set.seed(296)
phyglm_UN_hue_scale <- phyloglm( Urban ~ scale(Dichrom_hue) + scale(Mass_log), 
                                 data = UN_Hue_dat, 
                                 phy = UN_Hue_phy, 
                                 boot = 1000)

summary(phyglm_UN_hue_scale) # this model converges successfully

# save model
saveRDS(phyglm_UN_hue_scale, here("Models/UN", "phyglm_UN_hue_scale.rds"))
# load model
phyglm_UN_hue_scale <- readRDS(here("Models/UN", "phyglm_UN_hue_scale.rds"))


# get alpha, t, and half life for the model
(phyglm_UN_hue_scale$mean.tip.height) # t
(alpha_hue <- phyglm_UN_hue_scale$alpha) # alpha
(hl_hue <- log(2)/alpha_hue) # half life for the model 
# small half life compared with t -> low phylogenetic signal



##############################################################################
##############################################################################

############## UAI and and Sexual Selection Intensity on Males ###############

# create a new data frame that contains only species with both UAI and sexual selection intensity for males
UAI_SSM <- C_SSelect_dat2 %>% filter(!is.na(aveUAI))  %>% 
  filter(!is.na(sex.sel.m)) %>% as.data.frame()
length(UAI_SSM$sex.sel.m)
#766 species with UAI and sex.sel.m

###### add and pair tree

# add rownames to data
row.names(UAI_SSM) <- UAI_SSM$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_SSM_phydat <- treedata(tree_out, UAI_SSM, sort=T)

UAI_SSM_phy <- UAI_SSM_phydat$phy
UAI_SSM_dat <- as.data.frame(UAI_SSM_phydat$data)

str(UAI_SSM_dat)
length(UAI_SSM_dat$sex.sel.m)
#766


### convert traits of interest to numeric
UAI_SSM_dat$aveUAI <- as.numeric(UAI_SSM_dat$aveUAI)
UAI_SSM_dat$Mass_log <- as.numeric(UAI_SSM_dat$Mass_log)
UAI_SSM_dat$sex.sel.m <- as.numeric(UAI_SSM_dat$sex.sel.m)


# Run phylogenetic linear model
UAI_GLS_ssm <- gls(aveUAI~ sex.sel.m + Mass_log, data = UAI_SSM_dat, 
                   correlation = corPagel(0.5, phy = UAI_SSM_phy, fixed = F, form = ~Species_Jetz), 
                   method = "ML") 

# model summary and results
summary(UAI_GLS_ssm) 
confint(UAI_GLS_ssm)

# model diagnostics
check_model(UAI_GLS_ssm)
qqnorm(resid(UAI_GLS_ssm)) 
qqline(resid(UAI_GLS_ssm))
hist(resid(UAI_GLS_ssm)) 

# save model for easy retrieval 
saveRDS(UAI_GLS_ssm, here("Models/UAI", "UAI_GLS_ssm.rds"))

######################## MUTI and Sexual Selection Intensity on Males ##########################

# create a new data frame that contains only species with both MUTI and sexual selection intensity for males
MUTI_SSM <- C_SSelect_dat2 %>% filter(!is.na(MUTIscore)) %>% 
  filter(!is.na(sex.sel.m)) %>% as.data.frame()
length(MUTI_SSM$sex.sel.m)
#127 species with MUTI and sex.sel.m

colnames(MUTI_SSM)

###### add and pair tree

# add rownames to data
row.names(MUTI_SSM) <- MUTI_SSM$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_SSM_phydat <- treedata(tree_out, MUTI_SSM, sort=T)

MUTI_SSM_phy <- MUTI_SSM_phydat$phy
MUTI_SSM_dat <- as.data.frame(MUTI_SSM_phydat$data)

str(MUTI_SSM_dat)
length(MUTI_SSM_dat$sex.sel.m)
#127


### convert traits of interest to numeric
MUTI_SSM_dat$MUTIscore <- as.numeric(MUTI_SSM_dat$MUTIscore)
MUTI_SSM_dat$Mass_log <- as.numeric(MUTI_SSM_dat$Mass_log)
MUTI_SSM_dat$sex.sel.m <- as.numeric(MUTI_SSM_dat$sex.sel.m)


# Run phylogenetic linear model
MUTI_GLS_ssm <- gls(MUTIscore~ sex.sel.m + Mass_log, data = MUTI_SSM_dat, 
                   correlation = corPagel(0.5, phy = MUTI_SSM_phy, fixed=F, form = ~Species_Jetz), 
                   method = "ML") 


# model summary and results
summary(MUTI_GLS_ssm) 
confint(MUTI_GLS_ssm)

# model diagnostics
check_model(MUTI_GLS_ssm)
qqnorm(resid(MUTI_GLS_ssm)) 
qqline(resid(MUTI_GLS_ssm))
hist(resid(MUTI_GLS_ssm)) 

# save model for easy retrieval 
saveRDS(MUTI_GLS_ssm, here("Models/MUTI", "MUTI_GLS_ssm.rds"))


######################## UN and Sexual Selection Intensity on Males ##########################

# create a new data frame that contains only species with both UN and sexual selection intensity for males
UN_SSM <- C_SSelect_dat2 %>% filter(!is.na(Urban)) %>% 
  filter(!is.na(sex.sel.m)) %>% column_to_rownames(., var = "Species_Jetz")
length(UN_SSM$sex.sel.m)
# 129 species with UN and sex.sel.m


###### add and pair tree
tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UN_SSM_phydat <- treedata(tree_out, UN_SSM, sort=T)

UN_SSM_phy <- UN_SSM_phydat$phy
UN_SSM_dat <- as.data.frame(UN_SSM_phydat$data)

str(UN_SSM_dat)
length(UN_SSM_dat$sex.sel.m)
#129


### convert traits of interest to numeric
UN_SSM_dat$Urban <- as.numeric(UN_SSM_dat$Urban)
UN_SSM_dat$Mass_log <- as.numeric(UN_SSM_dat$Mass_log)
UN_SSM_dat$sex.sel.m <- as.numeric(UN_SSM_dat$sex.sel.m)


# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
phyglm_UN_ssm_scale <- phyloglm( Urban ~ scale(sex.sel.m) + scale(Mass_log), 
                                 data = UN_SSM_dat, 
                                 phy = UN_SSM_phy, 
                                 boot = 1000) 
# this fails to converge
summary(phyglm_UN_ssm_scale)


# print AIC values for models with different upper bounds
# intervals of 0.1 from 0 up to 4
for (i in seq(0, 4, by = 0.1)) {
  print(phyloglm(Urban ~ scale(sex.sel.m) + scale(Mass_log), 
                 data = UN_SSM_dat, 
                 phy = UN_SSM_phy, 
                 log.alpha.bound = i)$aic)
}
# AIC values support models with larger values of alpha (low phylo signal)


# try model with restricted alpha 
# give the model a little more searching space for alpha because AIC is actually lowest slightly below log.alpha.bounds = 4
exp(3.8)/(phyglm_UN_ssm_scale$mean.tip.height) # equals 0.45. Using this as start.alpha

set.seed(641)
phyglm_UN_ssm_fix <- phyloglm( Urban ~ scale(sex.sel.m) + scale(Mass_log), 
                               data = UN_SSM_dat, 
                               phy = UN_SSM_phy, 
                               log.alpha.bound = 4,
                               start.alpha = 0.45,
                               boot = 1000)
summary(phyglm_UN_ssm_fix)


# save model
saveRDS(phyglm_UN_ssm_fix, here("Models/UN", "phyglm_UN_ssm_fix.rds"))
# load model
phyglm_UN_ssm_fix <- readRDS(here("Models/UN", "phyglm_UN_ssm_fix.rds"))


# compare results with a non-phylogenetic logistic model
glm_UN_ssm <- logistf(Urban ~ scale(sex.sel.m) + scale(Mass_log), 
                      data = UN_SSM)
summary(glm_UN_ssm)


# get alpha, t, and half life for the model
(phyglm_UN_ssm_fix$mean.tip.height) # t
(alpha_ssm <- phyglm_UN_ssm_fix$alpha) # alpha
(hl_ssm <- log(2)/alpha_ssm) # half life  
# small half life relative to t



##############################################################################
##############################################################################


################## UAI and Sexual Selection Intensity on Females #############

# create a new data frame that contains only species with both UAI and sexual selection intensity for females
UAI_SSF <- C_SSelect_dat2 %>% filter(!is.na(aveUAI)) %>% 
  filter(!is.na(sex.sel.f)) %>% as.data.frame()
length(UAI_SSF$sex.sel.f)
#766 species with UAI and sex.sel.f

###### add and pair tree

# add rownames to data
row.names(UAI_SSF) <- UAI_SSF$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_SSF_phydat <- treedata(tree_out, UAI_SSF, sort=T)

UAI_SSF_phy <- UAI_SSF_phydat$phy
UAI_SSF_dat <- as.data.frame(UAI_SSF_phydat$data)

str(UAI_SSF_dat)
length(UAI_SSF_dat$sex.sel.f)
#766


### convert traits of interest to numeric
UAI_SSF_dat$aveUAI <- as.numeric(UAI_SSF_dat$aveUAI)
UAI_SSF_dat$Mass_log <- as.numeric(UAI_SSF_dat$Mass_log)
UAI_SSF_dat$sex.sel.f <- as.numeric(UAI_SSF_dat$sex.sel.f)


# Run phylogenetic linear model
UAI_GLS_ssf <- gls(aveUAI~ sex.sel.f + Mass_log, data = UAI_SSF_dat, 
                    correlation = corPagel(0.5, phy = UAI_SSF_phy, fixed=F, form = ~Species_Jetz), 
                    method = "ML") 

# model summary and results
summary(UAI_GLS_ssf) 
confint(UAI_GLS_ssf)

# model diagnostics
check_model(UAI_GLS_ssf)
qqnorm(resid(UAI_GLS_ssf)) 
qqline(resid(UAI_GLS_ssf))
hist(resid(UAI_GLS_ssf)) 

# save model for easy retrieval 
saveRDS(UAI_GLS_ssf, here("Models/UAI", "UAI_GLS_ssf.rds"))


######################## MUTI and Sexual Selection Intensity on Females ###################

# create a new data frame that contains only species with both MUTI and sexual selection intensity for females
MUTI_SSF <- C_SSelect_dat2 %>% filter(!is.na(MUTIscore)) %>% 
  filter(!is.na(sex.sel.f)) %>% as.data.frame()
length(MUTI_SSF$sex.sel.f)
#127 species with MUTI and sex.sel.f


###### add and pair tree

# add rownames to data
row.names(MUTI_SSF) <- MUTI_SSF$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_SSF_phydat <- treedata(tree_out, MUTI_SSF, sort=T)

MUTI_SSF_phy <- MUTI_SSF_phydat$phy
MUTI_SSF_dat <- as.data.frame(MUTI_SSF_phydat$data)

str(MUTI_SSF_dat)
length(MUTI_SSF_dat$sex.sel.f)
#127


### convert traits of interest to numeric
MUTI_SSF_dat$MUTIscore <- as.numeric(MUTI_SSF_dat$MUTIscore)
MUTI_SSF_dat$Mass_log <- as.numeric(MUTI_SSF_dat$Mass_log)
MUTI_SSF_dat$sex.sel.f <- as.numeric(MUTI_SSF_dat$sex.sel.f)


# Run phylogenetic linear model
MUTI_GLS_ssf <- gls(MUTIscore~ sex.sel.f + Mass_log, data = MUTI_SSF_dat, 
                   correlation = corPagel(0.5, phy = MUTI_SSF_phy, fixed=F, form = ~Species_Jetz), 
                   method = "ML") 

# model summary and results
summary(MUTI_GLS_ssf) 
confint(MUTI_GLS_ssf)

# model diagnostics
check_model(MUTI_GLS_ssf)
qqnorm(resid(MUTI_GLS_ssf)) 
qqline(resid(MUTI_GLS_ssf))
hist(resid(MUTI_GLS_ssf)) 

# save model for easy retrieval 
saveRDS(MUTI_GLS_ssf, here("Models/MUTI", "MUTI_GLS_ssf.rds"))


######################## UN and Sexual Selection Intensity on Females ##########################

# create a new data frame that contains only species with both UN and sexual selection intensity for females
UN_SSF <- C_SSelect_dat2 %>% filter(!is.na(Urban)) %>% 
  filter(!is.na(sex.sel.f)) %>% column_to_rownames(., var="Species_Jetz")
length(UN_SSF$sex.sel.f)
#129 species with UN and sex.sel.f


###### add and pair tree

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UN_SSF_phydat <- treedata(tree_out, UN_SSF, sort=T)

UN_SSF_phy <- UN_SSF_phydat$phy
UN_SSF_dat <- as.data.frame(UN_SSF_phydat$data)

str(UN_SSF_dat)
length(UN_SSF_dat$sex.sel.f)
#129


### convert traits of interest to numeric
UN_SSF_dat$Urban <- as.numeric(UN_SSF_dat$Urban)
UN_SSF_dat$Mass_log <- as.numeric(UN_SSF_dat$Mass_log)
UN_SSF_dat$sex.sel.f <- as.numeric(UN_SSF_dat$sex.sel.f)


# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
phyglm_UN_ssf_scale <- phyloglm( Urban ~ scale(sex.sel.f) + Mass_log, 
                                 data = UN_SSF_dat, 
                                 phy = UN_SSF_phy, 
                                 boot = 1000) 

summary(phyglm_UN_ssf_scale) 
# this version fails to converge
# alpha at upper bounds


# print AIC values for models with different upper bounds
# intervals of 0.1 from 0 up to 4
for (i in seq(0, 4, by = 0.1)) {
  print(phyloglm(Urban ~ scale(sex.sel.f) + scale(Mass_log), 
                 data = UN_SSF_dat, 
                 phy = UN_SSF_phy, 
                 log.alpha.bound = i)$aic)
}
# AIC values support models with larger values of alpha (low phylo signal)


# try fixing alpha at upper bounds
# had to go slightly above log.alpha.bound of 4 to get models that would converge
set.seed(786)
phyglm_UN_ssf_fix_4.05 <- phyloglm( Urban ~ scale(sex.sel.f) + scale(Mass_log), 
                                    data = UN_SSF_dat, 
                                    phy = UN_SSF_phy, 
                                    log.alpha.bound = 4.05,
                                    start.alpha = 0.56,
                                    boot = 1000
)

phyglm_UN_ssf_fix_4.1 <- phyloglm( Urban ~ scale(sex.sel.f) + scale(Mass_log), 
                                   data = UN_SSF_dat, 
                                   phy = UN_SSF_phy, 
                                   log.alpha.bound = 4.1,
                                   start.alpha = 0.6,
                                   boot = 1000
)
summary(phyglm_UN_ssf_fix_4.05)
summary(phyglm_UN_ssf_fix_4.1)
# both model converge and produce similar estimates


# save model
saveRDS(phyglm_UN_ssf_fix_4.05, here("Models/UN", "phyglm_UN_ssf_fix.rds"))
# load model
phyglm_UN_ssf_fix_4.05 <- readRDS(here("Models/UN", "phyglm_UN_ssf_fix.rds"))


# look at non-phylogenetic model
glm_UN_ssf <- logistf(Urban ~ scale(sex.sel.f) + scale(Mass_log), 
                      data = UN_SSF)
summary(glm_UN_ssf)
# agrees with two models above


# get alpha, t, and half life for the model
(phyglm_UN_ssf_fix_4.05$mean.tip.height) # t
(alpha_ssf <- phyglm_UN_ssf_fix_4.05$alpha) # alpha
(hl_ssf <- log(2)/alpha_ssf) # half life
# small half life compared to t

