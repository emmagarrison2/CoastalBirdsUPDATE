# The objective of this script is to run phylogenetic linear models for all 
# social traits. 

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



######################## UAI and territoriality ##########################

# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UAIDataUT <- C_Social_dat2 %>% filter(!is.na(aveUAI)) 
SocialData1 <- UAIDataUT %>% filter(!is.na(territoriality)) 
length(SocialData1$territoriality)
#766 species with UAI and territoriality

colnames(SocialData1)

###### add and pair tree

# add rownames to data
row.names(SocialData1) <- SocialData1$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Socialphydat1 <- treedata(tree_out,SocialData1, sort=T)

Socialphy1 <- Socialphydat1$phy
SocialTraitDat1 <- as.data.frame(Socialphydat1$data)

str(SocialTraitDat1)
length(SocialTraitDat1$territoriality)
#766


### convert traits of interest to numeric

SocialTraitDat1$aveUAI <- as.numeric(SocialTraitDat1$aveUAI)
SocialTraitDat1$Mass_log <- as.numeric(SocialTraitDat1$Mass_log)
SocialTraitDat1$territoriality <- as.numeric(SocialTraitDat1$territoriality)


#lets run the model!


UAI_GLS_territory <- gls(aveUAI~ territoriality + Mass_log, data = SocialTraitDat1, 
                      correlation = corPagel(0.5, phy=Socialphy1,fixed=F, form = ~Species_Jetz), 
                      method = "ML") 
#check out the model
check_model(UAI_GLS_territory) 
qqnorm(resid(UAI_GLS_territory)) 
qqline(resid(UAI_GLS_territory)) 
hist(resid(UAI_GLS_territory)) 


summary(UAI_GLS_territory) 
confint(UAI_GLS_territory)




######################## MUTI and territoriality ##########################

# lets first simplify a NEW database by removing records where we don't have an MUTI / brood_value
MUTIDataUT <- C_Social_dat2 %>% filter(!is.na(MUTIscore)) 
SocialData2 <- MUTIDataUT %>% filter(!is.na(territoriality)) 
length(SocialData2$territoriality)
#127 species with MUTI and territoriality

colnames(SocialData2)

###### add and pair tree

# add rownames to data
row.names(SocialData2) <- SocialData2$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Socialphydat2 <- treedata(tree_out,SocialData2, sort=T)

Socialphy2 <- Socialphydat2$phy
SocialTraitDat2 <- as.data.frame(Socialphydat2$data)

str(SocialTraitDat2)
length(SocialTraitDat2$territoriality)
#127


### convert traits of interest to numeric

SocialTraitDat2$MUTIscore <- as.numeric(SocialTraitDat2$MUTIscore)
SocialTraitDat2$Mass_log <- as.numeric(SocialTraitDat2$Mass_log)
SocialTraitDat2$territoriality <- as.numeric(SocialTraitDat2$territoriality)


#lets run the model!


MUTI_GLS_territory <- gls(MUTIscore~ territoriality + Mass_log, data = SocialTraitDat2, 
                         correlation = corPagel(0.5, phy=Socialphy2,fixed=F, form = ~Species_Jetz), 
                         method = "ML") 
#check out the model
check_model(MUTI_GLS_territory) 
qqnorm(resid(MUTI_GLS_territory)) 
qqline(resid(MUTI_GLS_territory)) 
hist(resid(MUTI_GLS_territory)) 


summary(MUTI_GLS_territory) 
confint(MUTI_GLS_territory)



######################## UN and territoriality ##########################

# lets first simplify a NEW database by removing records where we don't have an UN / brood_value
UNDataUT <- C_Social_dat2 %>% filter(!is.na(Urban)) 
SocialData3 <- UNDataUT %>% filter(!is.na(territoriality)) 
length(SocialData3$territoriality)
#129 species with UN and territoriality

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
#129


### convert traits of interest to numeric

SocialTraitDat3$Urban <- as.numeric(SocialTraitDat3$Urban)
SocialTraitDat3$Mass_log <- as.numeric(SocialTraitDat3$Mass_log)
SocialTraitDat3$territoriality <- as.numeric(SocialTraitDat3$territoriality)



#lets run the model using Phylolm!  

#(have to use lambda, until we figure out a way to make GLS work with binomial linear regression)
UN_M_territory <- phylolm(Urban~ territoriality + Mass_log, data=SocialTraitDat3,
                    phy=Socialphy3, model="lambda") 

# time to check out the model 
qqnorm(UN_M_territory$residuals)
qqline(UN_M_territory$residuals)
hist(UN_M_territory$residuals, breaks = 12) 

#lets get those values for our results table 
summary(UN_M_territory)
confint(UN_M_territory)





######################## UAI and cooperative ##########################

# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
colnames(C_Social_dat2)

UAIDataUT <- C_Social_dat2 %>% filter(!is.na(aveUAI)) 
SocialData4 <- UAIDataUT %>% filter(!is.na(cooperative)) 
length(SocialData4$cooperative)
#766 species with UAI and cooperative

colnames(SocialData4)

###### add and pair tree

# add rownames to data
row.names(SocialData4) <- SocialData4$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Socialphydat4 <- treedata(tree_out,SocialData4, sort=T)

Socialphy4 <- Socialphydat4$phy
SocialTraitDat4 <- as.data.frame(Socialphydat4$data)

str(SocialTraitDat4)
length(SocialTraitDat4$cooperative)
#766


### convert traits of interest to numeric

SocialTraitDat4$aveUAI <- as.numeric(SocialTraitDat4$aveUAI)
SocialTraitDat4$Mass_log <- as.numeric(SocialTraitDat4$Mass_log)
SocialTraitDat4$cooperative <- as.numeric(SocialTraitDat4$cooperative)



#lets run the model!


UAI_GLS_cooperative <- gls(aveUAI~ cooperative + Mass_log, data = SocialTraitDat4, 
                          correlation = corPagel(0.5, phy=Socialphy4,fixed=F, form = ~Species_Jetz), 
                          method = "ML") 
#check out the model
check_model(UAI_GLS_cooperative) 
qqnorm(resid(UAI_GLS_cooperative)) 
qqline(resid(UAI_GLS_cooperative)) 
hist(resid(UAI_GLS_cooperative)) 


summary(UAI_GLS_cooperative) 
confint(UAI_GLS_cooperative)




######################## MUTI and cooperative ##########################

# lets first simplify a NEW database by removing records where we don't have an MUTI / brood_value

MUTIDataUT <- C_Social_dat2 %>% filter(!is.na(MUTIscore)) 
SocialData5 <- MUTIDataUT %>% filter(!is.na(cooperative)) 
length(SocialData5$cooperative)
#127 species with MUTI and cooperative

colnames(SocialData5)

###### add and pair tree

# add rownames to data
row.names(SocialData5) <- SocialData5$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Socialphydat5 <- treedata(tree_out,SocialData5, sort=T)

Socialphy5 <- Socialphydat5$phy
SocialTraitDat5 <- as.data.frame(Socialphydat5$data)

str(SocialTraitDat5)
length(SocialTraitDat5$cooperative)
#127


### convert traits of interest to numeric

SocialTraitDat5$MUTIscore <- as.numeric(SocialTraitDat5$MUTIscore)
SocialTraitDat5$Mass_log <- as.numeric(SocialTraitDat5$Mass_log)
SocialTraitDat5$cooperative <- as.numeric(SocialTraitDat5$cooperative)



#lets run the model!


MUTI_GLS_cooperative <- gls(MUTIscore~ cooperative + Mass_log, data = SocialTraitDat5, 
                           correlation = corPagel(0.5, phy=Socialphy5,fixed=F, form = ~Species_Jetz), 
                           method = "ML") 
#check out the model
check_model(MUTI_GLS_cooperative) 
qqnorm(resid(MUTI_GLS_cooperative)) 
qqline(resid(MUTI_GLS_cooperative)) 
hist(resid(MUTI_GLS_cooperative)) 


summary(MUTI_GLS_cooperative) 
confint(MUTI_GLS_cooperative)






######################## UN and cooperative ##########################

# lets first simplify a NEW database by removing records where we don't have an UN / brood_value

UNDataUT <- C_Social_dat2 %>% filter(!is.na(Urban)) 
SocialData6 <- UNDataUT %>% filter(!is.na(cooperative)) 
length(SocialData6$cooperative)
#129 species with UN and cooperative

colnames(SocialData6)

###### add and pair tree

# add rownames to data
row.names(SocialData6) <- SocialData6$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Socialphydat6 <- treedata(tree_out,SocialData6, sort=T)

Socialphy6 <- Socialphydat6$phy
SocialTraitDat6 <- as.data.frame(Socialphydat6$data)

str(SocialTraitDat6)
length(SocialTraitDat6$cooperative)
#129


### convert traits of interest to numeric

SocialTraitDat6$Urban <- as.numeric(SocialTraitDat6$Urban)
SocialTraitDat6$Mass_log <- as.numeric(SocialTraitDat6$Mass_log)
SocialTraitDat6$cooperative <- as.numeric(SocialTraitDat6$cooperative)



#lets run the model using Phylolm!  

#(have to use lambda, until we figure out a way to make GLS work with binomial linear regression)
UN_M_cooperative <- phylolm(Urban~ cooperative + Mass_log, data=SocialTraitDat6,
                          phy=Socialphy6, model="lambda") 

# time to check out the model 
qqnorm(UN_M_cooperative$residuals)
qqline(UN_M_cooperative$residuals)
hist(UN_M_cooperative$residuals, breaks = 12) 

#lets get those values for our results table 
summary(UN_M_cooperative)
confint(UN_M_cooperative)
