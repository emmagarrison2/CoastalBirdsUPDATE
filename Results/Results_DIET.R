######Results - Diet Data Traits 

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


#load in "Coastal_Species_Diet.rds" - since this contains coastal species and all diet trait variables :) 

C_Diet_dat <- readRDS(here("Outputs", "Coastal_Species_Diet.rds"))
str(C_Diet_dat)

C_Diet_dat2 <- C_Diet_dat %>%
  mutate(Species_Jetz  = str_replace(Species_Jetz, " ", "_"))
str(C_Diet_dat2)

C_Diet_dat2$Urban <- ifelse(C_Diet_dat2$Urban == "U", 1, 0)
View(C_Diet_dat2)
colnames(C_Diet_dat2)


######################## UAI and % Diet Invertebrates ##########################


# lets first simplify a NEW database by removing records where we don't have an UAI / %dietinv
UAIDataUT <- C_Diet_dat2 %>% filter(!is.na(aveUAI)) 
DietData1 <- UAIDataUT %>% filter(!is.na(Diet.Inv)) 
length(DietData1$Diet.Inv)
#798 species with UAI and CT

###### add and pair tree

# add rownames to data
row.names(DietData1) <- DietData1$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Dietphydat1 <- treedata(tree_out,DietData1, sort=T)

Dietphy1 <- Dietphydat1$phy
DietTraitDat1 <- as.data.frame(Dietphydat1$data)

str(DietTraitDat1)
length(DietTraitDat1$Diet.Inv)
#798

### convert traits of interest to numeric

DietTraitDat1$aveUAI <- as.numeric(DietTraitDat1$aveUAI)
DietTraitDat1$Mass_log <- as.numeric(DietTraitDat1$Mass_log)
DietTraitDat1$Diet.Inv <- as.numeric(DietTraitDat1$Diet.Inv)


#lets run the model!


UAI_GLS_invert <- gls(aveUAI~ Diet.Inv + Mass_log, data = DietTraitDat1, 
                   correlation = corPagel(0.5, phy=Dietphy1,fixed=F, form = ~Species_Jetz), 
                   method = "ML") 
#check out the model
check_model(UAI_GLS_invert) ## low collinearity - which is good! - normality of residuals line does not fall on line, but is in a straight line 
qqnorm(resid(UAI_GLS_invert)) 
qqline(resid(UAI_GLS_invert)) #most points fall on the line 
hist(resid(UAI_GLS_invert)) #roughly normal dist of residuals


summary(UAI_GLS_invert) # strong phylogenetic relationship between traits and response, but no remaining influence of any predictor trait and response
confint(UAI_GLS_invert)





######################## MUTI and % Diet Invertebrates ##########################


# lets first simplify a NEW database by removing records where we don't have an UAI / %dietinv
MUTIDataUT <- C_Diet_dat2 %>% filter(!is.na(MUTIscore)) 
DietData2 <- MUTIDataUT %>% filter(!is.na(Diet.Inv)) 
length(DietData2$Diet.Inv)
#798 species with UAI and CT

###### add and pair tree

# add rownames to data
row.names(DietData2) <- DietData2$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Dietphydat2 <- treedata(tree_out,DietData2, sort=T)

Dietphy2 <- Dietphydat2$phy
DietTraitDat2 <- as.data.frame(Dietphydat2$data)

str(DietTraitDat2)
length(DietTraitDat2$Diet.Inv)
#130

### convert traits of interest to numeric

DietTraitDat2$MUTIscore <- as.numeric(DietTraitDat2$MUTIscore)
DietTraitDat2$Mass_log <- as.numeric(DietTraitDat2$Mass_log)
DietTraitDat2$Diet.Inv <- as.numeric(DietTraitDat2$Diet.Inv)



#lets run the model!


MUTI_GLS_invert <- gls(MUTIscore~ Diet.Inv + Mass_log, data = DietTraitDat2, 
                      correlation = corPagel(0.5, phy=Dietphy2,fixed=F, form = ~Species_Jetz), 
                      method = "ML") 
#check out the model
check_model(MUTI_GLS_invert) ## low collinearity - which is good! - normality of residuals line does not fall on line, but is in a straight line 
qqnorm(resid(MUTI_GLS_invert)) 
qqline(resid(MUTI_GLS_invert)) #most points fall on the line 
hist(resid(MUTI_GLS_invert)) #roughly normal dist of residuals


summary(MUTI_GLS_invert) # strong phylogenetic relationship between traits and response, but no remaining influence of any predictor trait and response
confint(MUTI_GLS_invert)







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



#lets run the model using Phylolm!  

#(have to use lambda, until we figure out a way to make GLS work with binomial linear regression)
UN_M_invert <- phylolm(Urban~ Diet.Inv + Mass_log, data=DietTraitDat3,
                   phy=Dietphy3, model="lambda") 

# time to check out the model 
qqnorm(UN_M_invert$residuals)
qqline(UN_M_invert$residuals) # what is happening? two separate lines bc of binomial... but is this the correct model check for binomial regression?
#the two lines do not have overlap... maybe this is good? 
hist(UN_M_invert$residuals, breaks = 20) 

#lets get those values for our results table 
summary(UN_M_invert)
confint(UN_M_invert)



######################## UAI and % Diet Vertebrates ##########################


# lets first simplify a NEW database by removing records where we don't have an UAI / %dietinv
UAIDataUT <- C798_Diet_dat2 %>% filter(!is.na(aveUAI)) 
DietData4 <- UAIDataUT %>% filter(!is.na(Diet.Vert)) 
length(DietData4$Diet.Vert)
#798 species with UAI and CT

###### add and pair tree

DietData4 <- as.data.frame(DietData4)
# add rownames to data
row.names(DietData4) <- DietData4$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Dietphydat4 <- treedata(tree_out,DietData4, sort=T)

Dietphy4 <- Dietphydat4$phy
DietTraitDat4 <- as.data.frame(Dietphydat4$data)

str(DietTraitDat4)
length(DietTraitDat4$Diet.Vert)
#798

### convert traits of interest to numeric

DietTraitDat4$aveUAI <- as.numeric(DietTraitDat4$aveUAI)
DietTraitDat4$Mass_log <- as.numeric(DietTraitDat4$Mass_log)
DietTraitDat4$Diet.Vert <- as.numeric(DietTraitDat4$Diet.Vert)



#lets run the model using GLS!


UAI_GLS_vert <- gls(aveUAI~ Diet.Vert + Mass_log, data = DietTraitDat4, 
                      correlation = corPagel(0.5, phy=Dietphy4,fixed=F, form = ~Species_Jetz), 
                      method = "ML") 
#check out the model
check_model(UAI_GLS_vert) ## low collinearity - which is good! - normality of residuals line does not fall on line, but is in a straight line 
qqnorm(resid(UAI_GLS_vert)) 
qqline(resid(UAI_GLS_vert)) #most points fall on the line 
hist(resid(UAI_GLS_vert)) #roughly normal dist of residuals


summary(UAI_GLS_vert) # strong phylogenetic relationship between traits and response, but no remaining influence of any predictor trait and response
confint(UAI_GLS_vert)





######################## MUTI and % Diet Vertebrates ##########################


# lets first simplify a NEW database by removing records where we don't have an UAI / %dietinv
MUTIDataUT <- C_Diet_dat2 %>% filter(!is.na(MUTIscore)) 
DietData5 <- MUTIDataUT %>% filter(!is.na(Diet.Vert)) 
length(DietData5$Diet.Vert)
#130 species with UAI and CT

###### add and pair tree

DietData5 <- as.data.frame(DietData5)
# add rownames to data
row.names(DietData5) <- DietData5$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Dietphydat5 <- treedata(tree_out,DietData5, sort=T)

Dietphy5 <- Dietphydat5$phy
DietTraitDat5 <- as.data.frame(Dietphydat5$data)

str(DietTraitDat5)
length(DietTraitDat5$Diet.Vert)
#130

### convert traits of interest to numeric

DietTraitDat5$MUTIscore <- as.numeric(DietTraitDat5$MUTIscore)
DietTraitDat5$Mass_log <- as.numeric(DietTraitDat5$Mass_log)
DietTraitDat5$Diet.Vert <- as.numeric(DietTraitDat5$Diet.Vert)



#lets run the model using GLS!


MUTI_GLS_vert <- gls(MUTIscore~ Diet.Vert + Mass_log, data = DietTraitDat5, 
                    correlation = corPagel(0.5, phy=Dietphy5,fixed=F, form = ~Species_Jetz), 
                    method = "ML") 
#check out the model
check_model(MUTI_GLS_vert) ## low collinearity - which is good! - normality of residuals line does not fall on line, but is in a straight line 
qqnorm(resid(MUTI_GLS_vert)) 
qqline(resid(MUTI_GLS_vert)) #most points fall on the line 
hist(resid(MUTI_GLS_vert)) #roughly normal dist of residuals


summary(MUTI_GLS_vert) # strong phylogenetic relationship between traits and response, but no remaining influence of any predictor trait and response
confint(MUTI_GLS_vert)







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


#lets run the model using Phylolm!  

#(have to use lambda, until we figure out a way to make GLS work with binomial linear regression)
UN_M_vert <- phylolm(Urban~ Diet.Vert + Mass_log, data=DietTraitDat6,
                       phy=Dietphy6, model="lambda") 

# time to check out the model 
qqnorm(UN_M_vert$residuals)
qqline(UN_M_vert$residuals) # what is happening? two separate lines bc of binomial... but is this the correct model check for binomial regression?
#the two lines do not have overlap... maybe this is good? 
hist(UN_M_vert$residuals, breaks = 20) 

#lets get those values for our results table 
summary(UN_M_vert)
confint(UN_M_vert)
