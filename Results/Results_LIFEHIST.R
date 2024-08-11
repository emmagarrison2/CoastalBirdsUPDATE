######Results - Life History Data Traits 

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


#load in "Coastal_Species_LifeHistory.rds" - since this contains coastal species and all diet trait variables :) 

C_LifeHist_dat <- readRDS(here("Outputs", "Coastal_Species_LifeHistory.rds"))
str(C_LifeHist_dat)

C_LifeHist_dat2 <- C_LifeHist_dat %>%
  mutate(Species_Jetz  = str_replace(Species_Jetz, " ", "_"))
str(C_LifeHist_dat2)

C_LifeHist_dat2$Urban <- ifelse(C_LifeHist_dat2$Urban == "U", 1, 0)
View(C_LifeHist_dat2)
colnames(C_LifeHist_dat2)



######################## UAI and % brood value ##########################


# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UAIDataUT <- C_LifeHist_dat2 %>% filter(!is.na(aveUAI)) 
LifehistData1 <- UAIDataUT %>% filter(!is.na(brood_value)) 
length(LifehistData1$brood_value)
#480 species with UAI and brood_value

###### add and pair tree

# add rownames to data
row.names(LifehistData1) <- LifehistData1$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Lifehistphydat1 <- treedata(tree_out,LifehistData1, sort=T)

Lifehistphy1 <- Lifehistphydat1$phy
LifehistTraitDat1 <- as.data.frame(Lifehistphydat1$data)

str(LifehistTraitDat1)
length(LifehistTraitDat1$brood_value)
#480

### convert traits of interest to numeric

LifehistTraitDat1$aveUAI <- as.numeric(LifehistTraitDat1$aveUAI)
LifehistTraitDat1$Mass_log <- as.numeric(LifehistTraitDat1$Mass_log)
LifehistTraitDat1$brood_value <- as.numeric(LifehistTraitDat1$brood_value)


#lets run the model!


UAI_GLS_bv <- gls(aveUAI~ brood_value + Mass_log, data = LifehistTraitDat1, 
                      correlation = corPagel(0.5, phy=Lifehistphy1,fixed=F, form = ~Species_Jetz), 
                      method = "ML") 
#check out the model
check_model(UAI_GLS_bv) ## low collinearity - which is good! - normality of residuals line does not fall on line, but is in a straight line 
qqnorm(resid(UAI_GLS_bv)) 
qqline(resid(UAI_GLS_bv)) #most points fall on the line 
hist(resid(UAI_GLS_bv)) #roughly normal dist of residuals


summary(UAI_GLS_bv) # strong phylogenetic relationship between traits and response, but no remaining influence of any predictor trait and response
confint(UAI_GLS_bv)




######################## MUTI and % brood value ##########################


# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
MUTIDataUT <- C_LifeHist_dat2 %>% filter(!is.na(MUTIscore)) 
LifehistData2 <- MUTIDataUT %>% filter(!is.na(brood_value)) 
length(LifehistData2$brood_value)
#121 species with UAI and brood_value

###### add and pair tree

# add rownames to data
row.names(LifehistData2) <- LifehistData2$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Lifehistphydat2 <- treedata(tree_out,LifehistData2, sort=T)

Lifehistphy2 <- Lifehistphydat2$phy
LifehistTraitDat2 <- as.data.frame(Lifehistphydat2$data)

str(LifehistTraitDat2)
length(LifehistTraitDat2$brood_value)
#121

### convert traits of interest to numeric

LifehistTraitDat2$MUTIscore <- as.numeric(LifehistTraitDat2$MUTIscore)
LifehistTraitDat2$Mass_log <- as.numeric(LifehistTraitDat2$Mass_log)
LifehistTraitDat2$brood_value <- as.numeric(LifehistTraitDat2$brood_value)


#lets run the model!


MUTI_GLS_bv <- gls(MUTIscore~ brood_value + Mass_log, data = LifehistTraitDat2, 
                  correlation = corPagel(0.5, phy=Lifehistphy2,fixed=F, form = ~Species_Jetz), 
                  method = "ML") 
#check out the model
check_model(MUTI_GLS_bv) ## low collinearity - which is good! - normality of residuals line does not fall on line, but is in a straight line 
qqnorm(resid(MUTI_GLS_bv)) 
qqline(resid(MUTI_GLS_bv)) #most points fall on the line 
hist(resid(MUTI_GLS_bv)) #roughly normal dist of residuals


summary(MUTI_GLS_bv) # strong phylogenetic relationship between traits and response, but no remaining influence of any predictor trait and response
confint(MUTI_GLS_bv)



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


#lets run the model using Phylolm!  

#(have to use lambda, until we figure out a way to make GLS work with binomial linear regression)
UN_M_bv <- phylolm(Urban~ brood_value + Mass_log, data=LifehistTraitDat3,
                   phy=Lifehistphy3, model="lambda") 

# time to check out the model 
qqnorm(UN_M_bv$residuals)
qqline(UN_M_bv$residuals) # what is happening? two separate lines bc of binomial... but is this the correct model check for binomial regression?
#the two lines do not have overlap... maybe this is good? 
hist(UN_M_bv$residuals, breaks = 20) 

#lets get those values for our results table 
summary(UN_M_bv)
confint(UN_M_bv)



######################## UAI and % clutch size ##########################


# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UAIDataUT <- C_LifeHist_dat2 %>% filter(!is.na(aveUAI)) 
LifehistData4 <- UAIDataUT %>% filter(!is.na(clutch_size)) 
length(LifehistData4$clutch_size)
#738 species with UAI and clutch_size

###### add and pair tree

# add rownames to data
LifehistData4 <- as.data.frame(LifehistData4)
row.names(LifehistData4) <- LifehistData4$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Lifehistphydat4 <- treedata(tree_out,LifehistData4, sort=T)

Lifehistphy4 <- Lifehistphydat4$phy
LifehistTraitDat4 <- as.data.frame(Lifehistphydat4$data)

str(LifehistTraitDat4)
length(LifehistTraitDat4$clutch_size)
#738

### convert traits of interest to numeric

LifehistTraitDat4$aveUAI <- as.numeric(LifehistTraitDat4$aveUAI)
LifehistTraitDat4$Mass_log <- as.numeric(LifehistTraitDat4$Mass_log)
LifehistTraitDat4$clutch_size <- as.numeric(LifehistTraitDat4$clutch_size)



#lets run the model!


UAI_GLS_clutch <- gls(aveUAI~ clutch_size + Mass_log, data = LifehistTraitDat4, 
                   correlation = corPagel(0.5, phy=Lifehistphy4,fixed=F, form = ~Species_Jetz), 
                   method = "ML") 
#check out the model
check_model(UAI_GLS_clutch) ## low collinearity - which is good! - normality of residuals line does not fall on line, but is in a straight line 
qqnorm(resid(UAI_GLS_clutch)) 
qqline(resid(UAI_GLS_clutch)) #most points fall on the line 
hist(resid(UAI_GLS_clutch)) #roughly normal dist of residuals


summary(UAI_GLS_clutch) # strong phylogenetic relationship between traits and response, but no remaining influence of any predictor trait and response
confint(UAI_GLS_clutch)



######################## MUTI and % clutch size ##########################


# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
MUTIDataUT <- C_LifeHist_dat2 %>% filter(!is.na(MUTIscore)) 
LifehistData5 <- MUTIDataUT %>% filter(!is.na(clutch_size)) 
length(LifehistData5$clutch_size)
#126 species with UAI and clutch_size

###### add and pair tree

# add rownames to data
LifehistData5 <- as.data.frame(LifehistData5)
row.names(LifehistData5) <- LifehistData5$Species_Jetz

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



#lets run the model!


MUTI_GLS_clutch <- gls(MUTIscore~ clutch_size + Mass_log, data = LifehistTraitDat5, 
                      correlation = corPagel(0.5, phy=Lifehistphy5,fixed=F, form = ~Species_Jetz), 
                      method = "ML") 
#check out the model
check_model(MUTI_GLS_clutch) ## low collinearity - which is good! - normality of residuals line does not fall on line, but is in a straight line 
qqnorm(resid(MUTI_GLS_clutch)) 
qqline(resid(MUTI_GLS_clutch)) #most points fall on the line 
hist(resid(MUTI_GLS_clutch)) #roughly normal dist of residuals


summary(MUTI_GLS_clutch) # strong phylogenetic relationship between traits and response, but no remaining influence of any predictor trait and response
confint(MUTI_GLS_clutch)



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



#lets run the model using Phylolm!  

#(have to use lambda, until we figure out a way to make GLS work with binomial linear regression)
UN_M_clutch <- phylolm(Urban~ clutch_size + Mass_log, data=LifehistTraitDat6,
                   phy=Lifehistphy6, model="lambda") 

# time to check out the model 
qqnorm(UN_M_clutch$residuals)
qqline(UN_M_clutch$residuals) # what is happening? two separate lines bc of binomial... but is this the correct model check for binomial regression?
#the two lines do not have overlap... maybe this is good? 
hist(UN_M_clutch$residuals, breaks = 20) 

#lets get those values for our results table 
summary(UN_M_clutch)
confint(UN_M_clutch)



######################## UAI and % longevity ##########################


# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UAIDataUT <- C_LifeHist_dat2 %>% filter(!is.na(aveUAI)) 
LifehistData7 <- UAIDataUT %>% filter(!is.na(longevity)) 
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



#lets run the model with GLS!


UAI_GLS_long <- gls(aveUAI~ longevity + Mass_log, data = LifehistTraitDat7, 
                      correlation = corPagel(0.5, phy=Lifehistphy7,fixed=F, form = ~Species_Jetz), 
                      method = "ML") 
#check out the model
check_model(UAI_GLS_long) 
qqnorm(resid(UAI_GLS_long)) 
qqline(resid(UAI_GLS_long)) 
hist(resid(UAI_GLS_long)) 


summary(UAI_GLS_long) 
confint(UAI_GLS_long)



######################## MUTI and % longevity ##########################


# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
MUTIDataUT <- C_LifeHist_dat2 %>% filter(!is.na(MUTIscore)) 
LifehistData8 <- MUTIDataUT %>% filter(!is.na(longevity)) 
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



#lets run the model with GLS!


MUTI_GLS_long <- gls(MUTIscore~ longevity + Mass_log, data = LifehistTraitDat8, 
                    correlation = corPagel(0.5, phy=Lifehistphy8,fixed=F, form = ~Species_Jetz), 
                    method = "ML") 
#check out the model
check_model(MUTI_GLS_long) 
qqnorm(resid(MUTI_GLS_long)) 
qqline(resid(MUTI_GLS_long)) 
hist(resid(MUTI_GLS_long)) 


summary(MUTI_GLS_long) 
confint(MUTI_GLS_long)




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





#lets run the model using Phylolm!  

#(have to use lambda, until we figure out a way to make GLS work with binomial linear regression)
UN_M_long <- phylolm(Urban~ longevity + Mass_log, data=LifehistTraitDat9,
                       phy=Lifehistphy9, model="lambda") 

# time to check out the model 
qqnorm(UN_M_long$residuals)
qqline(UN_M_long$residuals) # what is happening? two separate lines bc of binomial... but is this the correct model check for binomial regression?
#the two lines do not have overlap... maybe this is good? 
hist(UN_M_long$residuals, breaks = 20) 

#lets get those values for our results table 
summary(UN_M_long)
confint(UN_M_long)




######################## UAI and % developmental mode ##########################


# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UAIDataUT <- C_LifeHist_dat2 %>% filter(!is.na(aveUAI)) 
LifehistData10 <- UAIDataUT %>% filter(!is.na(developmental_mode)) 
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



#lets run the model with GLS!


UAI_GLS_develop <- gls(aveUAI~ developmental_mode + Mass_log, data = LifehistTraitDat10, 
                     correlation = corPagel(0.5, phy=Lifehistphy10,fixed=F, form = ~Species_Jetz), 
                     method = "ML") 
#check out the model
check_model(UAI_GLS_develop) 
qqnorm(resid(UAI_GLS_develop)) 
qqline(resid(UAI_GLS_develop)) 
hist(resid(UAI_GLS_develop)) 


summary(UAI_GLS_develop) 
confint(UAI_GLS_develop)



######################## MUTI and % developmental mode ##########################


# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
MUTIDataUT <- C_LifeHist_dat2 %>% filter(!is.na(MUTIscore)) 
LifehistData11 <- MUTIDataUT %>% filter(!is.na(developmental_mode)) 
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



#lets run the model with GLS!


MUTI_GLS_develop <- gls(MUTIscore~ developmental_mode + Mass_log, data = LifehistTraitDat11, 
                       correlation = corPagel(0.5, phy=Lifehistphy11,fixed=F, form = ~Species_Jetz), 
                       method = "ML") 
#check out the model
check_model(MUTI_GLS_develop) 
qqnorm(resid(MUTI_GLS_develop)) 
qqline(resid(MUTI_GLS_develop)) 
hist(resid(MUTI_GLS_develop)) 


summary(MUTI_GLS_develop) 
confint(MUTI_GLS_develop)



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




#lets run the model using Phylolm!  

#(have to use lambda, until we figure out a way to make GLS work with binomial linear regression)
UN_M_develop <- phylolm(Urban~ developmental_mode + Mass_log, data=LifehistTraitDat12,
                     phy=Lifehistphy12, model="lambda") 

# time to check out the model 
qqnorm(UN_M_develop$residuals)
qqline(UN_M_develop$residuals) # what is happening? two separate lines bc of binomial... but is this the correct model check for binomial regression?
#the two lines do not have overlap... maybe this is good? 
hist(UN_M_develop$residuals, breaks = 20) 

#lets get those values for our results table 
summary(UN_M_develop)
confint(UN_M_develop)


