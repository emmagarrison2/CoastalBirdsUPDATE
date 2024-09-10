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

######## method = default should be "logistic_MPLE" ######## 

phyglm_UN_Mass <- phyloglm( Urban ~ Mass_log, 
                               data = MassTraitDat3, 
                               phy = Massphy3, 
                               boot = 1000) # just add the boot = argument to the function

#Warning message:
  #In phyloglm(Urban ~ Mass_log, data = MassTraitDat3, phy = Massphy3,  :
                #the estimate of 'alpha' (0.555065022314764) reached the upper bound (0.559633190138608).
             # This may simply reflect a flat likelihood at large alpha values,
             # meaning that the phylogenetic correlation is estimated to be negligible.


# time to check out the model 
qqnorm(phyglm_UN_Mass$residuals)
qqline(phyglm_UN_Mass$residuals)
hist(phyglm_UN_Mass$residuals, breaks = 12) 

#lets get those values for our results table 
summary(phyglm_UN_Mass)
confint(phyglm_UN_Mass)

# z value = 1.177
# SE = 0.114
# p value = 0.239
# 95% CI = (-0.0890, 0.357)

#let's calculate alpha and half-life as measurements of phylogenetic signal 

alpha <- phyglm_UN_Mass$alpha # where phyglm_UN_Mass is the name of your model
alpha #555065
log(2)/alpha # Half-Life for the model = 1.249
#compared to T of 97.561, this is a small Half-Life (small phylogenetic signal) 


###################### SENSORY TRAITS ######################
#######################################################
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


######## method = logistic_MPLE ######## 


phyglm_UN_CT <- phyloglm( Urban ~ C.T + Mass_log, 
                              data = SensoryTraitDat2, 
                              phy = Sensoryphy2, 
                              boot = 1000) # just add the boot = argument to the function

#Warning message:
 # In phyloglm(Urban ~ C.T + Mass_log, data = SensoryTraitDat2, phy = Sensoryphy2,  :
  #              the boundary of the linear predictor has been reached during the optimization procedure.
   #           You can increase this bound by increasing 'btol'.

# time to check out the model 
qqnorm(phyglm_UN_CT$residuals)
qqline(phyglm_UN_CT$residuals)
hist(phyglm_UN_CT$residuals, breaks = 20) 

#lets get those values for our results table 
summary(phyglm_UN_CT)
confint(phyglm_UN_CT)

alpha <- phyglm_UN_CT$alpha 
alpha # 0.5172797 
log(2)/alpha # Half-Life for the model = 1.339985 
#compared to T of 97.561, this is a small Half-Life 





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


######## method = logistic_MPLE ######## 


phyglm_UN_pf <- phyloglm( Urban ~ peak_freq + Mass_log, 
                          data = SensoryTraitDat5, 
                          phy = Sensoryphy5, 
                          boot = 1000) # just add the boot = argument to the function

#Warning message:
# In phyloglm(Urban ~ C.T + Mass_log, data = SensoryTraitDat2, phy = Sensoryphy2,  :
#              the boundary of the linear predictor has been reached during the optimization procedure.
#           You can increase this bound by increasing 'btol'.

# time to check out the model 
qqnorm(phyglm_UN_CT$residuals)
qqline(phyglm_UN_CT$residuals)
hist(phyglm_UN_CT$residuals, breaks = 20) 

#lets get those values for our results table 
summary(phyglm_UN_CT)
confint(phyglm_UN_CT)

alpha <- phyglm_UN_CT$alpha 
alpha # 0.5172797 
log(2)/alpha # Half-Life for the model = 1.339985 
#compared to T of 97.561, this is a small Half-Life 





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

# let's run a model using Phyloglm 


######## method = logistic_MPLE ######## 


phyglm_UN_Invert <- phyloglm( Urban ~ Diet.Inv + Mass_log, 
                            data = DietTraitDat3, 
                            phy = Dietphy3, 
                            boot = 1000) # just add the boot = argument to the function

#Warning message:
  #In phyloglm(Urban ~ Diet.Inv + Mass_log, data = DietTraitDat3, phy = Dietphy3,  :
                #the boundary of the linear predictor has been reached during the optimization procedure.
             # You can increase this bound by increasing 'btol'.

# time to check out the model 
qqnorm(phyglm_UN_Invert$residuals)
qqline(phyglm_UN_Invert$residuals)
hist(phyglm_UN_Invert$residuals, breaks = 20) 

#lets get those values for our results table 
summary(phyglm_UN_Invert)
confint(phyglm_UN_Invert)

# z value = 1.080
# SE = 0.00880
# p value = 0.280
# 95% CI = (-0.00774, 0.0267)

alpha <- phyglm_UN_Invert$alpha 
alpha #0.0107
log(2)/alpha # Half-Life for the model = 65.057
#compared to T of 97.561, this is a small-moderate Half-Life 




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


######## method = logistic_MPLE ######## 


phyglm_UN_Vert <- phyloglm( Urban ~ Diet.Vert + Mass_log, 
                              data = DietTraitDat6, 
                              phy = Dietphy6, 
                              boot = 1000) # just add the boot = argument to the function

#Warning message:
 # In phyloglm(Urban ~ Diet.Vert + Mass_log, data = DietTraitDat6,  :
  #              the boundary of the linear predictor has been reached during the optimization procedure.
   #           You can increase this bound by increasing 'btol'.

# time to check out the model 
qqnorm(phyglm_UN_Vert$residuals)
qqline(phyglm_UN_Vert$residuals)
hist(phyglm_UN_Vert$residuals, breaks = 20) 

#lets get those values for our results table 
summary(phyglm_UN_Vert)
confint(phyglm_UN_Vert)

alpha <- phyglm_UN_Vert$alpha 
alpha #0.0103
log(2)/alpha # Half-Life for the model = 67.579
#compared to T of 97.561, this is a small-moderate Half-Life




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


######## method = logistic_MPLE ######## 


phyglm_UN_PS <- phyloglm( Urban ~ Diet.PS + Mass_log, 
                            data = DietTraitDat9, 
                            phy = Dietphy9, 
                            boot = 1000) # just add the boot = argument to the function

#Warning message:
 # In phyloglm(Urban ~ Diet.PS + Mass_log, data = DietTraitDat9, phy = Dietphy9,  :
  #              the boundary of the linear predictor has been reached during the optimization procedure.
   #           You can increase this bound by increasing 'btol'.



# time to check out the model 
qqnorm(phyglm_UN_PS$residuals)
qqline(phyglm_UN_PS$residuals)
hist(phyglm_UN_PS$residuals, breaks = 20) 

#lets get those values for our results table 
summary(phyglm_UN_PS)
confint(phyglm_UN_PS)

alpha <- phyglm_UN_PS$alpha 
alpha # 0.0108
log(2)/alpha # Half-Life for the model = 63.920
#compared to T of 97.561, this is a small-moderate Half-Life




######################## UN and % Diet Fruit/Nut ##########################


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


######## method = logistic_MPLE ######## 


phyglm_UN_FN <- phyloglm( Urban ~ Diet.FN + Mass_log, 
                          data = DietTraitDat12, 
                          phy = Dietphy12, 
                          boot = 1000) # just add the boot = argument to the function

#Warning message:
 # In phyloglm(Urban ~ Diet.FN + Mass_log, data = DietTraitDat12, phy = Dietphy12,  :
  #              the boundary of the linear predictor has been reached during the optimization procedure.
   #           You can increase this bound by increasing 'btol'.

# time to check out the model 
qqnorm(phyglm_UN_FN$residuals)
qqline(phyglm_UN_FN$residuals)
hist(phyglm_UN_FN$residuals, breaks = 20) 

#lets get those values for our results table 
summary(phyglm_UN_FN)
confint(phyglm_UN_FN)

alpha <- phyglm_UN_FN$alpha 
alpha # 0.5436619 
log(2)/alpha # Half-Life for the model = 1.27496
#compared to T of 97.561, this is a small  Half-Life






###################### LIFE HISTORY TRAITS ######################
#######################################################
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

######## method = logistic_MPLE --> the default ######## 

phyglm_UN_bv <- phyloglm( Urban ~ brood_value + Mass_log, 
                          data = LifehistTraitDat3, 
                          phy = Lifehistphy3, 
                          boot = 1000) # just add the boot = argument to the function

#Warning messages:
 # 1: In phyloglm(Urban ~ brood_value + Mass_log, data = LifehistTraitDat3,  :
  #                 the estimate of 'alpha' (0.559459062146366) reached the upper bound (0.559633190137694).
   #              This may simply reflect a flat likelihood at large alpha values,
    #             meaning that the phylogenetic correlation is estimated to be negligible.
     #            2: In phyloglm(Urban ~ brood_value + Mass_log, data = LifehistTraitDat3,  :
      #                            phyloglm failed to converge.

# time to check out the model 
qqnorm(phyglm_UN_bv$residuals)
qqline(phyglm_UN_bv$residuals)
hist(phyglm_UN_bv$residuals, breaks = 20) 

#lets get those values for our results table 
summary(phyglm_UN_bv)
confint(phyglm_UN_bv)

alpha <- phyglm_UN_bv$alpha 
alpha #  
log(2)/alpha # Half-Life for the model = 
#compared to T of 97.561, this is a ______  Half-Life




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

######## method = logistic_MPLE --> the default ######## 

phyglm_UN_clutch <- phyloglm( Urban ~ clutch_size + Mass_log, 
                          data = LifehistTraitDat6, 
                          phy = Lifehistphy6, 
                          boot = 1000) # just add the boot = argument to the function

#Warning messages:
 # 1: In phyloglm(Urban ~ clutch_size + Mass_log, data = LifehistTraitDat6,  :
  #                 the estimate of 'alpha' (0.559269614382704) reached the upper bound (0.559633190137426).
   #              This may simply reflect a flat likelihood at large alpha values,
    #             meaning that the phylogenetic correlation is estimated to be negligible.
     #            2: In phyloglm(Urban ~ clutch_size + Mass_log, data = LifehistTraitDat6,  :
      #                            phyloglm failed to converge.

# time to check out the model 
qqnorm(phyglm_UN_clutch$residuals)
qqline(phyglm_UN_clutch$residuals)
hist(phyglm_UN_clutch$residuals, breaks = 20) 

#lets get those values for our results table 
summary(phyglm_UN_clutch)
confint(phyglm_UN_clutch)

alpha <- phyglm_UN_clutch$alpha 
alpha #  
log(2)/alpha # Half-Life for the model = 
#compared to T of 97.561, this is a ______  Half-Life





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

######## method = logistic_MPLE --> the default ######## 

phyglm_UN_long <- phyloglm( Urban ~ longevity + Mass_log, 
                              data = LifehistTraitDat9, 
                              phy = Lifehistphy9, 
                              boot = 1000) # just add the boot = argument to the function

#Warning message:
 # In phyloglm(Urban ~ longevity + Mass_log, data = LifehistTraitDat9,  :
  #              the boundary of the linear predictor has been reached during the optimization procedure.
   #           You can increase this bound by increasing 'btol'.

# time to check out the model 
qqnorm(phyglm_UN_long$residuals)
qqline(phyglm_UN_long$residuals)
hist(phyglm_UN_long$residuals, breaks = 20) 

#lets get those values for our results table 
summary(phyglm_UN_long)
confint(phyglm_UN_long)

alpha <- phyglm_UN_long$alpha 
alpha #  0.01031167 
log(2)/alpha # Half-Life for the model = 67.21971 
#compared to T of 97.561, this is a small-moderate  Half-Life




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

######## method = logistic_MPLE --> the default ######## 

phyglm_UN_develop <- phyloglm( Urban ~ developmental_mode + Mass_log, 
                            data = LifehistTraitDat12, 
                            phy = Lifehistphy12, 
                            boot = 1000) # just add the boot = argument to the function

#Warning message:
 # In phyloglm(Urban ~ developmental_mode + Mass_log, data = LifehistTraitDat12,  :
  #              the estimate of 'alpha' (0.559597400576619) reached the upper bound (0.559633190138608).
   #           This may simply reflect a flat likelihood at large alpha values,
    #          meaning that the phylogenetic correlation is estimated to be negligible.

# time to check out the model 
qqnorm(phyglm_UN_develop$residuals)
qqline(phyglm_UN_develop$residuals)
hist(phyglm_UN_develop$residuals, breaks = 20) 

#lets get those values for our results table 
summary(phyglm_UN_develop)
confint(phyglm_UN_develop)

alpha <- phyglm_UN_develop$alpha 
alpha # 0.5595974
log(2)/alpha # Half-Life for the model =   1.238653 
#compared to T of 97.561, this is a small Half-Life






###################### NESTING TRAITS ######################
#######################################################
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


######## method = logistic_MPLE --> the default ######## 

phyglm_UN_nest_low <- phyloglm( Urban ~ NestSite_Low + Mass_log, 
                               data = NestTraitDat6, 
                               phy = Nestphy6, 
                               boot = 1000) # just add the boot = argument to the function

#Warning message:
#In phyloglm(Urban ~ NestSite_Low + Mass_log, data = NestTraitDat6,  :
 #               the estimate of 'alpha' (0.557838282933223) reached the upper bound (0.559633190138608).
  #            This may simply reflect a flat likelihood at large alpha values,
   #           meaning that the phylogenetic correlation is estimated to be negligible.

# time to check out the model 
qqnorm(phyglm_UN_nest_low$residuals)
qqline(phyglm_UN_nest_low$residuals)
hist(phyglm_UN_nest_low$residuals, breaks = 20) 

#lets get those values for our results table 
summary(phyglm_UN_nest_low)
confint(phyglm_UN_nest_low)

alpha <- phyglm_UN_nest_low$alpha 
alpha # 0.5578383 
log(2)/alpha # Half-Life for the model = 1.242559 
#compared to T of 97.561, this is a small Half-Life


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


######## method = logistic_MPLE --> the default ######## 

phyglm_UN_nest_low_only <- phyloglm( Urban ~ NestSite_Low + Mass_log, 
                                data = NestTraitDat6_only, 
                                phy = Nestphy6_only, 
                                boot = 1000) # just add the boot = argument to the function


# time to check out the model 
qqnorm(phyglm_UN_nest_low_only$residuals)
qqline(phyglm_UN_nest_low_only$residuals)
hist(phyglm_UN_nest_low_only$residuals, breaks = 20) 

#lets get those values for our results table 
summary(phyglm_UN_nest_low_only)
confint(phyglm_UN_nest_low_only)

alpha <- phyglm_UN_nest_low_only$alpha 
alpha # 0.04336139 
log(2)/alpha # Half-Life for the model = 15.98536 
#compared to T of 97.561, this is a small-moderate Half-Life





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


######## method = logistic_MPLE --> the default ######## 

phyglm_UN_nest_high <- phyloglm( Urban ~ NestSite_High + Mass_log, 
                                data = NestTraitDat9, 
                                phy = Nestphy9, 
                                boot = 1000) # just add the boot = argument to the function

#Warning message:
 # In phyloglm(Urban ~ NestSite_High + Mass_log, data = NestTraitDat9,  :
  #              the estimate of 'alpha' (0.559269401187482) reached the upper bound (0.559633190138608).
   #           This may simply reflect a flat likelihood at large alpha values,
    #          meaning that the phylogenetic correlation is estimated to be negligible.

# time to check out the model 
qqnorm(phyglm_UN_nest_high$residuals)
qqline(phyglm_UN_nest_high$residuals)
hist(phyglm_UN_nest_high$residuals, breaks = 20) 

#lets get those values for our results table 
summary(phyglm_UN_nest_high)
confint(phyglm_UN_nest_high)

alpha <- phyglm_UN_nest_high$alpha 
alpha #  1.23938 
log(2)/alpha # Half-Life for the model =  1.23938  
#compared to T of 97.561, this is a small Half-Life


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


######## method = logistic_MPLE --> the default ######## 

phyglm_UN_nest_high_only <- phyloglm( Urban ~ NestSite_High + Mass_log, 
                                     data = NestTraitDat9_only, 
                                     phy = Nestphy9_only, 
                                     boot = 1000) # just add the boot = argument to the function
#Warning message:
 # In phyloglm(Urban ~ NestSite_High + Mass_log, data = NestTraitDat9_only,  :
  #              the estimate of 'alpha' (0.558408354005091) reached the upper bound (0.55963319013772).
   #           This may simply reflect a flat likelihood at large alpha values,
    #          meaning that the phylogenetic correlation is estimated to be negligible.

# time to check out the model 
qqnorm(phyglm_UN_nest_high_only$residuals)
qqline(phyglm_UN_nest_high_only$residuals)
hist(phyglm_UN_nest_high_only$residuals, breaks = 20) 

#lets get those values for our results table 
summary(phyglm_UN_nest_high_only)
confint(phyglm_UN_nest_high_only)

alpha <- phyglm_UN_nest_high_only$alpha 
alpha # 0.5584084  
log(2)/alpha # Half-Life for the model =  1.241291 
#compared to T of 97.561, this is a small Half-Life





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


######## method = logistic_MPLE --> the default ######## 

phyglm_UN_nest_safety <- phyloglm( Urban ~ nest.safety + Mass_log, 
                                 data = NestTraitDat12, 
                                 phy = Nestphy12, 
                                 boot = 1000) # just add the boot = argument to the function

#Warning messages:
 # 1: In phyloglm(Urban ~ nest.safety + Mass_log, data = NestTraitDat12,  :
  #                 the estimate of 'alpha' (0.559398725309794) reached the upper bound (0.559633190138608).
   #              This may simply reflect a flat likelihood at large alpha values,
    #             meaning that the phylogenetic correlation is estimated to be negligible.
     #            2: In phyloglm(Urban ~ nest.safety + Mass_log, data = NestTraitDat12,  :
      #                            the boundary of the linear predictor has been reached during the optimization procedure.
       #                         You can increase this bound by increasing 'btol'.
        #                        3: In phyloglm(Urban ~ nest.safety + Mass_log, data = NestTraitDat12,  :
         #                                        phyloglm failed to converge.

# time to check out the model 
qqnorm(phyglm_UN_nest_safety$residuals)
qqline(phyglm_UN_nest_safety$residuals)
hist(phyglm_UN_nest_safety$residuals, breaks = 20) 

#lets get those values for our results table 
summary(phyglm_UN_nest_safety)
confint(phyglm_UN_nest_safety)

alpha <- phyglm_UN_nest_safety$alpha 
alpha #  ____ 
log(2)/alpha # Half-Life for the model =    
#compared to T of 97.561, this is a ____ Half-Life






###################### SEXUAL SELECTION TRAITS ######################
#######################################################
#######################################################



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


######## method = logistic_MPLE --> the default ######## 

phyglm_UN_brightness <- phyloglm( Urban ~ Dichrom_bright + Mass_log, 
                                   data = SSTraitDat3, 
                                   phy = SSphy3, 
                                   boot = 1000) # just add the boot = argument to the function

#Warning message:
 # In phyloglm(Urban ~ Dichrom_bright + Mass_log, data = SSTraitDat3,  :
  #              the boundary of the linear predictor has been reached during the optimization procedure.
   #           You can increase this bound by increasing 'btol'.

# time to check out the model 
qqnorm(phyglm_UN_brightness$residuals)
qqline(phyglm_UN_brightness$residuals)
hist(phyglm_UN_brightness$residuals, breaks = 20) 

#lets get those values for our results table 
summary(phyglm_UN_brightness)
confint(phyglm_UN_brightness)

alpha <- phyglm_UN_brightness$alpha 
alpha #  0.04296851  
log(2)/alpha # Half-Life for the model =  16.13151   
#compared to T of 97.561, this is a small/moderate Half-Life





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


######## method = logistic_MPLE --> the default ######## 

phyglm_UN_hue <- phyloglm( Urban ~ Dichrom_hue + Mass_log, 
                                  data = SSTraitDat6, 
                                  phy = SSphy6, 
                                  boot = 1000) # just add the boot = argument to the function

#Warning message:
# In phyloglm(Urban ~ NestSite_High + Mass_log, data = NestTraitDat9,  :
#              the estimate of 'alpha' (0.559269401187482) reached the upper bound (0.559633190138608).
#           This may simply reflect a flat likelihood at large alpha values,
#          meaning that the phylogenetic correlation is estimated to be negligible.

# time to check out the model 
qqnorm(phyglm_UN_hue$residuals)
qqline(phyglm_UN_hue$residuals)
hist(phyglm_UN_hue$residuals, breaks = 20) 

#lets get those values for our results table 
summary(phyglm_UN_hue)
confint(phyglm_UN_hue)

alpha <- phyglm_UN_hue$alpha 
alpha #  0.4830611  
log(2)/alpha # Half-Life for the model =   1.434906  
#compared to T of 97.561, this is a small Half-Life




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


######## method = logistic_MPLE --> the default ######## 

phyglm_UN_ssm <- phyloglm( Urban ~ sex.sel.m + Mass_log, 
                           data = SSTraitDat9, 
                           phy = SSphy9, 
                           boot = 1000) # just add the boot = argument to the function

#Warning messages:
 # 1: In phyloglm(Urban ~ sex.sel.m + Mass_log, data = SSTraitDat9, phy = SSphy9,  :
  #                 the estimate of 'alpha' (0.557824566925717) reached the upper bound (0.559633190138608).
   #              This may simply reflect a flat likelihood at large alpha values,
    #             meaning that the phylogenetic correlation is estimated to be negligible.
     #            2: In phyloglm(Urban ~ sex.sel.m + Mass_log, data = SSTraitDat9, phy = SSphy9,  :
      #                            phyloglm failed to converge.

# time to check out the model 
qqnorm(phyglm_UN_ssm$residuals)
qqline(phyglm_UN_ssm$residuals)
hist(phyglm_UN_ssm$residuals, breaks = 20) 

#lets get those values for our results table 
summary(phyglm_UN_ssm)
confint(phyglm_UN_ssm)

alpha <- phyglm_UN_ssm$alpha 
alpha #  ____  
log(2)/alpha # Half-Life for the model =    
#compared to T of 97.561, this is a ____ Half-Life



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


######## method = logistic_MPLE --> the default ######## 

phyglm_UN_ssf <- phyloglm( Urban ~ sex.sel.f + Mass_log, 
                           data = SSTraitDat12, 
                           phy = SSphy12, 
                           boot = 1000) # just add the boot = argument to the function

#Warning messages:
 # 1: In phyloglm(Urban ~ sex.sel.f + Mass_log, data = SSTraitDat12, phy = SSphy12,  :
  #                 the estimate of 'alpha' (0.559342431481888) reached the upper bound (0.559633190138608).
   #              This may simply reflect a flat likelihood at large alpha values,
    #             meaning that the phylogenetic correlation is estimated to be negligible.
     #            2: In phyloglm(Urban ~ sex.sel.f + Mass_log, data = SSTraitDat12, phy = SSphy12,  :
      #                            phyloglm failed to converge.

# time to check out the model 
qqnorm(phyglm_UN_ssf$residuals)
qqline(phyglm_UN_ssf$residuals)
hist(phyglm_UN_ssf$residuals, breaks = 20) 

#lets get those values for our results table 
summary(phyglm_UN_ssf)
confint(phyglm_UN_ssf)

alpha <- phyglm_UN_ssf$alpha 
alpha #  ____  
log(2)/alpha # Half-Life for the model =    
#compared to T of 97.561, this is a ____ Half-Life





###################### SEXUAL SELECTION TRAITS ######################
#######################################################
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


######## method = logistic_MPLE --> the default ######## 

phyglm_UN_territorial <- phyloglm( Urban ~ territoriality + Mass_log, 
                           data = SocialTraitDat3, 
                           phy = Socialphy3, 
                           boot = 1000) # just add the boot = argument to the function


# time to check out the model 
qqnorm(phyglm_UN_territorial$residuals)
qqline(phyglm_UN_territorial$residuals)
hist(phyglm_UN_territorial$residuals, breaks = 20) 

#lets get those values for our results table 
summary(phyglm_UN_territorial)
confint(phyglm_UN_territorial)

alpha <- phyglm_UN_territorial$alpha 
alpha #  0.07617762   
log(2)/alpha # Half-Life for the model = 9.099092 
#compared to T of 97.561, this is a small  Half-Life




######################## UN and cooperative ##########################

# do not do --> too small of N for Cooperative = 1 (only N = 8)


