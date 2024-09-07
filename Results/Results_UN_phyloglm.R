# The objective of this script is to run UN models with phyloglm(), 
# which is a logistic regression, NOT linear. This is more fitting 
# for binary response variables. 

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
library(phylolm)

###################### BODY MASS ######################
#######################################################
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

######## method = logistic_IG10 ######## 

UN_M_mass_IG10 <- phyloglm(Urban ~ Mass_log, 
                               data = MassTraitDat3, 
                               phy = Massphy3, 
                               method = "logistic_IG10")
#Warning message:
  #In phyloglm(Urban ~ Mass_log, data = MassTraitDat3, phy = Massphy3,  :
               # the boundary of the linear predictor has been reached during the optimization procedure.
              #You can increase this bound by increasing 'btol'.

# time to check out the model 
qqnorm(UN_M_mass_IG10$residuals)
qqline(UN_M_mass_IG10$residuals)
hist(UN_M_mass_IG10$residuals, breaks = 12) 

#lets get those values for our results table 
summary(UN_M_mass_IG10)
confint(UN_M_mass_IG10)

# z value = -0.116
# SE = 0.114
# p value = 0.907
# 95% CI = (-0.236, 0.209)
# alpha = 0.545


######## method = logistic_MPLE ######## 

UN_M_mass_MLPE <- phyloglm(Urban ~ Mass_log, 
                           data = MassTraitDat3, 
                           phy = Massphy3, 
                           method = "logistic_MPLE")
#Warning message:
  #In phyloglm(Urban ~ Mass_log, data = MassTraitDat3, phy = Massphy3,  :
                #the estimate of 'alpha' (0.555065022314764) reached the upper bound (0.559633190138608).
             # This may simply reflect a flat likelihood at large alpha values,
             # meaning that the phylogenetic correlation is estimated to be negligible.


# time to check out the model 
qqnorm(UN_M_mass_MLPE$residuals)
qqline(UN_M_mass_MLPE$residuals)
hist(UN_M_mass_MLPE$residuals, breaks = 12) 

#lets get those values for our results table 
summary(UN_M_mass_MLPE)
confint(UN_M_mass_MLPE)

# z value = 1.177
# SE = 0.114
# p value = 0.239
# 95% CI = (-0.0890, 0.357)
# alpha = 0.555



###################### DIET TRAITS ######################
#######################################################
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

# let's run a model using phyloglm 

######## method = logistic_IG10 ######## 

UN_M_invert_IG10 <- phyloglm(Urban ~ Diet.Inv + Mass_log, 
                               data = DietTraitDat3, 
                               phy = Dietphy3, 
                               method = "logistic_IG10")
#Warning message:
  #In phyloglm(Urban ~ Diet.Inv + Mass_log, data = DietTraitDat3, phy = Dietphy3,  :
               #phyloglm failed to converge.

# time to check out the model 
qqnorm(UN_M_invert_IG10$residuals)
qqline(UN_M_invert_IG10$residuals)
hist(UN_M_invert_IG10$residuals, breaks = 20) 

#lets get those values for our results table 
summary(UN_M_invert_IG10)
confint(UN_M_invert_IG10)

# z value = 0.0218
# SE = 0.00826
# p value = 0.983
# 95% CI = (-0.0160, 0.0164)
# alpha = 0.528

######## method = logistic_MPLE ######## 

UN_M_invert_MPLE <- phyloglm(Urban ~ Diet.Inv + Mass_log, 
                             data = DietTraitDat3, 
                             phy = Dietphy3, 
                             method = "logistic_MPLE")
#Warning message:
  #In phyloglm(Urban ~ Diet.Inv + Mass_log, data = DietTraitDat3, phy = Dietphy3,  :
                #the boundary of the linear predictor has been reached during the optimization procedure.
             # You can increase this bound by increasing 'btol'.

# time to check out the model 
qqnorm(UN_M_invert_MPLE$residuals)
qqline(UN_M_invert_MPLE$residuals)
hist(UN_M_invert_MPLE$residuals, breaks = 20) 

#lets get those values for our results table 
summary(UN_M_invert_MPLE)
confint(UN_M_invert_MPLE)

# z value = 1.080
# SE = 0.00880
# p value = 0.280
# 95% CI = (-0.00774, 0.0267)
# alpha = 0.0107

