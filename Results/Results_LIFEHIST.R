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

Lifehistphydat5 <- treedata(tree_out,LifehistData5, sort=T)

Lifehistphy5 <- Lifehistphydat5$phy
LifehistTraitDat5 <- as.data.frame(Lifehistphydat5$data)

str(LifehistTraitDat5)
length(LifehistTraitDat5$clutch_size)
#126

### convert traits of interest to numeric
LifehistTraitDat5$MUTIscore <- as.numeric(LifehistTraitDat5$MUTIscore)
LifehistTraitDat5$Mass_log <- as.numeric(LifehistTraitDat5$Mass_log)
LifehistTraitDat5$clutch_size <- as.numeric(LifehistTraitDat5$clutch_size)


# Run phylogenetic linear model
MUTI_GLS_clutch <- gls(MUTIscore~ clutch_size + Mass_log, data = LifehistTraitDat5, 
                      correlation = corPagel(0.5, phy=Lifehistphy5,fixed=F, form = ~Species_Jetz), 
                      method = "ML") 


# model diagnostics
check_model(MUTI_GLS_clutch) ## low collinearity - which is good! - normality of residuals line does not fall on line, but is in a straight line 
qqnorm(resid(MUTI_GLS_clutch)) 
qqline(resid(MUTI_GLS_clutch)) #most points fall on the line 
hist(resid(MUTI_GLS_clutch)) #roughly normal dist of residuals

# model summary and results
summary(MUTI_GLS_clutch) # strong phylogenetic relationship between traits and response, but no remaining influence of any predictor trait and response
confint(MUTI_GLS_clutch)

# save model for easy retrieval 
saveRDS(MUTI_GLS_clutch, here("Models/MUTI", "MUTI_GLS_clutch.rds"))


######################## UN and Clutch Size ##########################


# create a new data frame that contains only species with both UN and clutch size
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



#############################################################################
#############################################################################

######################## UAI and Longevity ##########################


# create a new data frame that contains only species with both UAI and longevity
UAIDataUT <- C_LifeHist_dat2 %>% filter(!is.na(aveUAI)) 
LifehistData7 <- UAIDataUT %>% filter(!is.na(longevity)) %>% as.data.frame()
length(LifehistData7$longevity)
#796 species with UAI and longevity

###### add and pair tree

# add rownames to data
LifehistData7 <- as.data.frame(LifehistData7)
row.names(LifehistData7) <- LifehistData7$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Lifehistphydat7 <- treedata(tree_out,LifehistData7, sort=T)

Lifehistphy7 <- Lifehistphydat7$phy
LifehistTraitDat7 <- as.data.frame(Lifehistphydat7$data)

str(LifehistTraitDat7)
length(LifehistTraitDat7$longevity)
#796

### convert traits of interest to numeric
LifehistTraitDat7$aveUAI <- as.numeric(LifehistTraitDat7$aveUAI)
LifehistTraitDat7$Mass_log <- as.numeric(LifehistTraitDat7$Mass_log)
LifehistTraitDat7$longevity <- as.numeric(LifehistTraitDat7$longevity)


# Run phylogenetic linear model
UAI_GLS_long <- gls(aveUAI~ longevity + Mass_log, data = LifehistTraitDat7, 
                      correlation = corPagel(0.5, phy=Lifehistphy7,fixed=F, form = ~Species_Jetz), 
                      method = "ML") 

# model diagnostics
check_model(UAI_GLS_long) 
qqnorm(resid(UAI_GLS_long)) 
qqline(resid(UAI_GLS_long)) 
hist(resid(UAI_GLS_long)) 

# model summary and results
summary(UAI_GLS_long) 
confint(UAI_GLS_long)



#because 0.05 < p value < 0.20 , let's check the 85% CI to see if this relationship is still notable 
Confidence_Interval_85 <- confint(UAI_GLS_long, level = 0.85)
Confidence_Interval_85

#crosses 0... let's not consider this model as "notable" in our results 

# save model for easy retrieval 
saveRDS(UAI_GLS_long, here("Models/UAI", "UAI_GLS_long.rds"))

######################## MUTI and Longevity ##########################


# create a new data frame that contains only species with both MUTI and longevity
MUTIDataUT <- C_LifeHist_dat2 %>% filter(!is.na(MUTIscore)) 
LifehistData8 <- MUTIDataUT %>% filter(!is.na(longevity)) %>% as.data.frame()
length(LifehistData8$longevity)
#130 species with UAI and longevity

###### add and pair tree

# add rownames to data
LifehistData8 <- as.data.frame(LifehistData8)
row.names(LifehistData8) <- LifehistData8$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Lifehistphydat8 <- treedata(tree_out,LifehistData8, sort=T)

Lifehistphy8 <- Lifehistphydat8$phy
LifehistTraitDat8 <- as.data.frame(Lifehistphydat8$data)

str(LifehistTraitDat8)
length(LifehistTraitDat8$longevity)
#130

### convert traits of interest to numeric
LifehistTraitDat8$MUTIscore <- as.numeric(LifehistTraitDat8$MUTIscore)
LifehistTraitDat8$Mass_log <- as.numeric(LifehistTraitDat8$Mass_log)
LifehistTraitDat8$longevity <- as.numeric(LifehistTraitDat8$longevity)


# Run phylogenetic linear model
MUTI_GLS_long <- gls(MUTIscore~ longevity + Mass_log, data = LifehistTraitDat8, 
                    correlation = corPagel(0.5, phy=Lifehistphy8,fixed=F, form = ~Species_Jetz), 
                    method = "ML") 

# model diagnostics
check_model(MUTI_GLS_long) 
qqnorm(resid(MUTI_GLS_long)) 
qqline(resid(MUTI_GLS_long)) 
hist(resid(MUTI_GLS_long)) 

# model summary and results
summary(MUTI_GLS_long) 
confint(MUTI_GLS_long)


# save model for easy retrieval 
saveRDS(MUTI_GLS_long, here("Models/MUTI", "MUTI_GLS_long.rds"))

######################## UN and Longevity ##########################


# create a new data frame that contains only species with both UN and longevity
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







#############################################################################
#############################################################################


######################## UAI and Developmental Mode ##########################


# create a new data frame that contains only species with both UAI and developmental mode
UAIDataUT <- C_LifeHist_dat2 %>% filter(!is.na(aveUAI)) 
LifehistData10 <- UAIDataUT %>% filter(!is.na(developmental_mode)) %>% as.data.frame()
length(LifehistData10$developmental_mode)
#766 species with UAI and developmental_mode

###### add and pair tree

# add rownames to data
LifehistData10 <- as.data.frame(LifehistData10)
row.names(LifehistData10) <- LifehistData10$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Lifehistphydat10 <- treedata(tree_out,LifehistData10, sort=T)

Lifehistphy10 <- Lifehistphydat10$phy
LifehistTraitDat10 <- as.data.frame(Lifehistphydat10$data)

str(LifehistTraitDat10)
length(LifehistTraitDat10$developmental_mode)
#766

### convert traits of interest to numeric
LifehistTraitDat10$aveUAI <- as.numeric(LifehistTraitDat10$aveUAI)
LifehistTraitDat10$Mass_log <- as.numeric(LifehistTraitDat10$Mass_log)
LifehistTraitDat10$developmental_mode <- as.numeric(LifehistTraitDat10$developmental_mode)



# Run phylogenetic linear model
UAI_GLS_develop <- gls(aveUAI~ developmental_mode + Mass_log, data = LifehistTraitDat10, 
                     correlation = corPagel(0.5, phy=Lifehistphy10,fixed=F, form = ~Species_Jetz), 
                     method = "ML") 

# model diagnostics
check_model(UAI_GLS_develop) 
qqnorm(resid(UAI_GLS_develop)) 
qqline(resid(UAI_GLS_develop)) 
hist(resid(UAI_GLS_develop)) 

# model summary and results
summary(UAI_GLS_develop) 
confint(UAI_GLS_develop)


# save model for easy retrieval 
saveRDS(UAI_GLS_develop, here("Models/UAI", "UAI_GLS_develop.rds"))

######################## MUTI and Developmental Mode ##########################


# create a new data frame that contains only species with both MUTI and developmental mode
MUTIDataUT <- C_LifeHist_dat2 %>% filter(!is.na(MUTIscore)) 
LifehistData11 <- MUTIDataUT %>% filter(!is.na(developmental_mode)) %>% as.data.frame()
length(LifehistData11$developmental_mode)
#127 species with MUTI and developmental_mode

###### add and pair tree

# add rownames to data
LifehistData11 <- as.data.frame(LifehistData11)
row.names(LifehistData11) <- LifehistData11$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Lifehistphydat11 <- treedata(tree_out,LifehistData11, sort=T)

Lifehistphy11 <- Lifehistphydat11$phy
LifehistTraitDat11 <- as.data.frame(Lifehistphydat11$data)

str(LifehistTraitDat11)
length(LifehistTraitDat11$developmental_mode)
#127

### convert traits of interest to numeric

LifehistTraitDat11$MUTIscore <- as.numeric(LifehistTraitDat11$MUTIscore)
LifehistTraitDat11$Mass_log <- as.numeric(LifehistTraitDat11$Mass_log)
LifehistTraitDat11$developmental_mode <- as.numeric(LifehistTraitDat11$developmental_mode)



# Run phylogenetic linear model
MUTI_GLS_develop <- gls(MUTIscore~ developmental_mode + Mass_log, data = LifehistTraitDat11, 
                       correlation = corPagel(0.5, phy=Lifehistphy11,fixed=F, form = ~Species_Jetz), 
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




