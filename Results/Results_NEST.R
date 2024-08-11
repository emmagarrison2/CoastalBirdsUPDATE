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


#load in "Coastal_Species_Nest.rds" - since this contains coastal species and all diet trait variables :) 

C_Nest_dat <- readRDS(here("Outputs", "Coastal_Species_Nest.rds"))
str(C_Nest_dat)

C_Nest_dat2 <- C_Nest_dat %>%
  mutate(Species_Jetz  = str_replace(Species_Jetz, " ", "_"))
str(C_Nest_dat2)

C_Nest_dat2$Urban <- ifelse(C_Nest_dat2$Urban == "U", 1, 0)
View(C_Nest_dat2)
colnames(C_Nest_dat2)



######################## UAI and %  Nest Strategy ##########################
# 0 = enclosed
# 1 = open

# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UAIDataUT <- C_Nest_dat2 %>% filter(!is.na(aveUAI)) 
NestData1 <- UAIDataUT %>% filter(!is.na(NestStr)) 
length(NestData1$NestStr)
#733 species with UAI and NestStr

#there is no Log_mass column in this .rds ... let's put one in! 

# add column for log transformed body mass
NestData1 <- NestData1 %>%
  mutate(Mass_log = log(Mass))

colnames(NestData1)

###### add and pair tree

# add rownames to data
row.names(NestData1) <- NestData1$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Nestphydat1 <- treedata(tree_out,NestData1, sort=T)

Nestphy1 <- Nestphydat1$phy
NestTraitDat1 <- as.data.frame(Nestphydat1$data)

str(NestTraitDat1)
length(NestTraitDat1$NestStr)
View(NestTraitDat1)
#733


### convert traits of interest to numeric

NestTraitDat1$aveUAI <- as.numeric(NestTraitDat1$aveUAI)
NestTraitDat1$Mass_log <- as.numeric(NestTraitDat1$Mass_log)
NestTraitDat1$NestStr <- as.numeric(NestTraitDat1$NestStr)


#lets run the model!


UAI_GLS_neststr <- gls(aveUAI~ NestStr + Mass_log, data = NestTraitDat1, 
                  correlation = corPagel(0.5, phy=Nestphy1,fixed=F, form = ~Species_Jetz), 
                  method = "ML") 
#check out the model
check_model(UAI_GLS_neststr) ## low collinearity - which is good! - normality of residuals line does not fall on line, but is in a straight line 
qqnorm(resid(UAI_GLS_neststr)) 
qqline(resid(UAI_GLS_neststr)) #most points fall on the line 
hist(resid(UAI_GLS_neststr)) #roughly normal dist of residuals


summary(UAI_GLS_neststr) # strong phylogenetic relationship between traits and response, but no remaining influence of any predictor trait and response
confint(UAI_GLS_neststr)



######################## MUTI and %  Nest Strategy ##########################
# 0 = enclosed
# 1 = open

# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
MUTIDataUT <- C_Nest_dat2 %>% filter(!is.na(MUTIscore)) 
NestData2 <- MUTIDataUT %>% filter(!is.na(NestStr)) 
length(NestData2$NestStr)
#117 species with MUTIscore and NestStr

#there is no Log_mass column in this .rds ... let's put one in! 

# add column for log transformed body mass
NestData2 <- NestData2 %>%
  mutate(Mass_log = log(Mass))

colnames(NestData2)

###### add and pair tree

# add rownames to data
row.names(NestData2) <- NestData2$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Nestphydat2 <- treedata(tree_out,NestData2, sort=T)

Nestphy2 <- Nestphydat2$phy
NestTraitDat2 <- as.data.frame(Nestphydat2$data)

str(NestTraitDat2)
length(NestTraitDat2$NestStr)
#117


### convert traits of interest to numeric

NestTraitDat2$MUTIscore <- as.numeric(NestTraitDat2$MUTIscore)
NestTraitDat2$Mass_log <- as.numeric(NestTraitDat2$Mass_log)
NestTraitDat2$NestStr <- as.numeric(NestTraitDat2$NestStr)


#lets run the model!


MUTI_GLS_neststr <- gls(MUTIscore~ NestStr + Mass_log, data = NestTraitDat2, 
                       correlation = corPagel(0.5, phy=Nestphy2,fixed=F, form = ~Species_Jetz), 
                       method = "ML") 
#check out the model
check_model(MUTI_GLS_neststr) ## low collinearity - which is good! - normality of residuals line does not fall on line, but is in a straight line 
qqnorm(resid(MUTI_GLS_neststr)) 
qqline(resid(MUTI_GLS_neststr)) #most points fall on the line 
hist(resid(MUTI_GLS_neststr)) #roughly normal dist of residuals


summary(MUTI_GLS_neststr) # strong phylogenetic relationship between traits and response, but no remaining influence of any predictor trait and response
confint(MUTI_GLS_neststr)




######################## UN and %  Nest Strategy ##########################
# 0 = enclosed
# 1 = open

# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UNDataUT <- C_Nest_dat2 %>% filter(!is.na(Urban)) 
NestData3 <- UNDataUT %>% filter(!is.na(NestStr)) 
length(NestData3$NestStr)
#122 species with Urban and NestStr

#there is no Log_mass column in this .rds ... let's put one in! 

# add column for log transformed body mass
NestData3 <- NestData3 %>%
  mutate(Mass_log = log(Mass))

colnames(NestData3)

###### add and pair tree

# add rownames to data
row.names(NestData3) <- NestData3$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Nestphydat3 <- treedata(tree_out,NestData3, sort=T)

Nestphy3 <- Nestphydat3$phy
NestTraitDat3 <- as.data.frame(Nestphydat3$data)

str(NestTraitDat3)
length(NestTraitDat3$NestStr)
#122


### convert traits of interest to numeric

NestTraitDat3$Urban <- as.numeric(NestTraitDat3$Urban)
NestTraitDat3$Mass_log <- as.numeric(NestTraitDat3$Mass_log)
NestTraitDat3$NestStr <- as.numeric(NestTraitDat3$NestStr)


#lets run the model using Phylolm!  

#(have to use lambda, until we figure out a way to make GLS work with binomial linear regression)
UN_M_neststr <- phylolm(Urban~ NestStr + Mass_log, data=NestTraitDat3,
                   phy=Nestphy3, model="lambda") 

# time to check out the model 
qqnorm(UN_M_neststr$residuals)
qqline(UN_M_neststr$residuals) # what is happening? two separate lines bc of binomial... but is this the correct model check for binomial regression?
#the two lines do not have overlap... maybe this is good? 
hist(UN_M_neststr$residuals, breaks = 20) 

#lets get those values for our results table 
summary(UN_M_neststr)
confint(UN_M_neststr)


######################## UAI and %  Nest Site LOW ##########################
# 0 = not low
# 1 = LOW

# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UAIDataUT <- C_Nest_dat2 %>% filter(!is.na(aveUAI)) 
NestData4 <- UAIDataUT %>% filter(!is.na(NestSite_Low)) 
length(NestData4$NestStr)
#796 species with UAI and NestSite_Low

#there is no Log_mass column in this .rds ... let's put one in! 

# add column for log transformed body mass
NestData4 <- NestData4 %>%
  mutate(Mass_log = log(Mass))

colnames(NestData4)

###### add and pair tree

# add rownames to data
row.names(NestData4) <- NestData4$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Nestphydat4 <- treedata(tree_out,NestData4, sort=T)

Nestphy4 <- Nestphydat4$phy
NestTraitDat4 <- as.data.frame(Nestphydat4$data)

str(NestTraitDat4)
length(NestTraitDat4$NestSite_Low)
#796


### convert traits of interest to numeric

NestTraitDat4$aveUAI <- as.numeric(NestTraitDat4$aveUAI)
NestTraitDat4$Mass_log <- as.numeric(NestTraitDat4$Mass_log)
NestTraitDat4$NestSite_Low <- as.numeric(NestTraitDat4$NestSite_Low)


#lets run the model!


UAI_GLS_nest_low <- gls(aveUAI~ NestSite_Low + Mass_log, data = NestTraitDat4, 
                       correlation = corPagel(0.5, phy=Nestphy4,fixed=F, form = ~Species_Jetz), 
                       method = "ML") 
#check out the model
check_model(UAI_GLS_nest_low) 
qqnorm(resid(UAI_GLS_nest_low)) 
qqline(resid(UAI_GLS_nest_low))
hist(resid(UAI_GLS_nest_low))


summary(UAI_GLS_nest_low) 
confint(UAI_GLS_nest_low)




######################## MUTI and %  Nest Site LOW ##########################
# 0 = not low
# 1 = LOW

# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UAIDataUT <- C_Nest_dat2 %>% filter(!is.na(MUTIscore)) 
NestData5 <- UAIDataUT %>% filter(!is.na(NestSite_Low)) 
length(NestData5$NestSite_Low)
#130 species with MUTIscore and NestSite_Low

#there is no Log_mass column in this .rds ... let's put one in! 

# add column for log transformed body mass
NestData5 <- NestData5 %>%
  mutate(Mass_log = log(Mass))

colnames(NestData5)

###### add and pair tree

# add rownames to data
row.names(NestData5) <- NestData5$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Nestphydat5 <- treedata(tree_out,NestData5, sort=T)

Nestphy5 <- Nestphydat5$phy
NestTraitDat5 <- as.data.frame(Nestphydat5$data)

str(NestTraitDat5)
length(NestTraitDat5$NestSite_Low)
#130


### convert traits of interest to numeric

NestTraitDat5$MUTIscore <- as.numeric(NestTraitDat5$MUTIscore)
NestTraitDat5$Mass_log <- as.numeric(NestTraitDat5$Mass_log)
NestTraitDat5$NestSite_Low <- as.numeric(NestTraitDat5$NestSite_Low)


#lets run the model!


MUTI_GLS_nest_low <- gls(MUTIscore~ NestSite_Low + Mass_log, data = NestTraitDat5, 
                        correlation = corPagel(0.5, phy=Nestphy5,fixed=F, form = ~Species_Jetz), 
                        method = "ML") 
#check out the model
check_model(MUTI_GLS_nest_low) 
qqnorm(resid(MUTI_GLS_nest_low)) 
qqline(resid(MUTI_GLS_nest_low))
hist(resid(MUTI_GLS_nest_low))


summary(MUTI_GLS_nest_low) 
confint(MUTI_GLS_nest_low)




######################## UN and %  Nest Site LOW ##########################
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


#lets run the model using Phylolm!  

#(have to use lambda, until we figure out a way to make GLS work with binomial linear regression)
UN_M_nest_low <- phylolm(Urban~ NestStr + Mass_log, data=NestTraitDat6,
                        phy=Nestphy6, model="lambda") 

# time to check out the model 
qqnorm(UN_M_nest_low$residuals)
qqline(UN_M_nest_low$residuals) # what is happening? two separate lines bc of binomial... but is this the correct model check for binomial regression?
#the two lines do not have overlap... maybe this is good? 
hist(UN_M_nest_low$residuals, breaks = 20) 

#lets get those values for our results table 
summary(UN_M_nest_low)
confint(UN_M_nest_low)

