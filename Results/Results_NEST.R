# The objective of this script is to run phylogenetic linear models for all nesting traits. 

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

######################## UAI and %  Nest Strategy ##########################
# 0 = enclosed
# 1 = open

# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UAIDataUT <- C_Nest_dat2 %>% filter(!is.na(aveUAI)) 
NestData1 <- UAIDataUT %>% filter(!is.na(NestStr)) 
length(NestData1$NestStr)
#833 species with UAI and NestStr

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
#833


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
#118 species with MUTIscore and NestStr

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
#118


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
length(NestData4$NestSite_Low)
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

str(NestTraitDat4$aveUAI)


#lets run the model!


UAI_GLS_nest_low_2 <- gls(aveUAI~ NestSite_Low + Mass_log, data = NestTraitDat4, 
                       correlation = corPagel(0.5, phy=Nestphy4,fixed=F, form = ~Species_Jetz), 
                       method = "ML") 
#check out the model
check_model(UAI_GLS_nest_low_2) 
qqnorm(resid(UAI_GLS_nest_low_2)) 
qqline(resid(UAI_GLS_nest_low_2))
hist(resid(UAI_GLS_nest_low_2))


summary(UAI_GLS_nest_low_2) 
confint(UAI_GLS_nest_low_2)

#NOW, let's try filtering out all species that have both HIGH and LOW nesting scores of 1 
#### idea - are the "non-ambiguous" nesting species the ones driving this significant relationship 
#### between nesting and UAI index? 

#filter 
UAIDataUT <- C_Nest_dat2 %>% filter(!is.na(aveUAI)) 
NestData4_only <- UAIDataUT %>% filter(!is.na(NestSite_Low)) %>% 
  filter(!(NestSite_Low == 1 & NestSite_High == 1))
length(NestData4_only$NestSite_Low)
#575 species with UAI and ONLY NestSite_Low (not also high nesters)
796-575 # = 221 --> there are 221 species that were both High and Low nesters and had UAI scores 

#there is no Log_mass column in this .rds ... let's put one in! 

# add column for log transformed body mass
NestData4_only <- NestData4_only %>%
  mutate(Mass_log = log(Mass))

colnames(NestData4_only)

###### add and pair tree

# add rownames to data
row.names(NestData4_only) <- NestData4_only$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Nestphydat4_only <- treedata(tree_out,NestData4_only, sort=T)

Nestphy4_only <- Nestphydat4_only$phy
NestTraitDat4_only <- as.data.frame(Nestphydat4_only$data)

str(NestTraitDat4_only)
length(NestTraitDat4_only$NestSite_Low)
#575


### convert traits of interest to numeric
NestTraitDat4_only$aveUAI <- as.numeric(NestTraitDat4_only$aveUAI)
NestTraitDat4_only$Mass_log <- as.numeric(NestTraitDat4_only$Mass_log)
NestTraitDat4_only$NestSite_Low <- as.numeric(NestTraitDat4_only$NestSite_Low)


#lets run the model!


UAI_GLS_nest_low_only <- gls(aveUAI~ NestSite_Low + Mass_log, data = NestTraitDat4_only, 
                          correlation = corPagel(0.5, phy=Nestphy4_only,fixed=F, form = ~Species_Jetz), 
                          method = "ML") 
#check out the model
check_model(UAI_GLS_nest_low_only) 
qqnorm(resid(UAI_GLS_nest_low_only)) 
qqline(resid(UAI_GLS_nest_low_only))
hist(resid(UAI_GLS_nest_low_only))


summary(UAI_GLS_nest_low_only) 
confint(UAI_GLS_nest_low_only)





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




#NOW, let's try filtering out all species that have both HIGH and LOW nesting scores of 1 
#### idea - are the "non-ambiguous" nesting species the ones driving this significant relationship 
#### between nesting and UAI index? 

#filter 
MUTIDataUT <- C_Nest_dat2 %>% filter(!is.na(MUTIscore)) 
NestData5_only <- MUTIDataUT %>% filter(!is.na(NestSite_Low)) %>% 
  filter(!(NestSite_Low == 1 & NestSite_High == 1))
length(NestData5_only$NestSite_Low)
#91 species with UAI and ONLY NestSite_Low (not also high nesters)
130-91 # = 39 --> there are 39 species that were both High and Low nesters and had MUTI scores 

#there is no Log_mass column in this .rds ... let's put one in! 

# add column for log transformed body mass
NestData5_only <- NestData5_only %>%
  mutate(Mass_log = log(Mass))

colnames(NestData5_only)

###### add and pair tree

# add rownames to data
row.names(NestData5_only) <- NestData5_only$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Nestphydat5_only <- treedata(tree_out,NestData5_only, sort=T)

Nestphy5_only <- Nestphydat5_only$phy
NestTraitDat5_only <- as.data.frame(Nestphydat5_only$data)

str(NestTraitDat5_only)
length(NestTraitDat5_only$NestSite_Low)
#91


### convert traits of interest to numeric

NestTraitDat5_only$MUTIscore <- as.numeric(NestTraitDat5_only$MUTIscore)
NestTraitDat5_only$Mass_log <- as.numeric(NestTraitDat5_only$Mass_log)
NestTraitDat5_only$NestSite_Low <- as.numeric(NestTraitDat5_only$NestSite_Low)



########### since model not working as corPagel starting point = 0.5 and fixed = F... 
# let's find a
# lambda value to fix in the model 


# Create an empty vector to store AIC values
AIC_values <- numeric()

# Loop through different values of the parameter for corPagel
for (i in seq(0, 1, by = 0.1)) {
  # Fit the gls model with the current value of i
  model <- gls(MUTIscore ~ NestSite_Low + Mass_log, 
               data = NestTraitDat5_only, 
               correlation = corPagel(i, phy = Nestphy5_only, fixed = TRUE, form = ~Species_Jetz), 
               method = "ML")
  # Extract AIC value and store it in the vector
  AIC_values <- c(AIC_values, AIC(model))
}

# Print AIC values
print(AIC_values)
#0.3 = best AIC score 



#lets run the model!


MUTI_GLS_nest_low_only <- gls(MUTIscore~ NestSite_Low + Mass_log, data = NestTraitDat5_only, 
                         correlation = corPagel(0.3, phy=Nestphy5_only,fixed=T, form = ~Species_Jetz), 
                         method = "ML") 
#check out the model
check_model(MUTI_GLS_nest_low_only) 
qqnorm(resid(MUTI_GLS_nest_low_only)) 
qqline(resid(MUTI_GLS_nest_low_only))
hist(resid(MUTI_GLS_nest_low_only))


summary(MUTI_GLS_nest_low_only) 
confint(MUTI_GLS_nest_low_only)




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
UN_M_nest_low <- phylolm(Urban~ NestSite_Low + Mass_log, data=NestTraitDat6,
                        phy=Nestphy6, model="lambda") 

# time to check out the model 
qqnorm(UN_M_nest_low$residuals)
qqline(UN_M_nest_low$residuals) # what is happening? two separate lines bc of binomial... but is this the correct model check for binomial regression?
#the two lines do not have overlap... maybe this is good? 
hist(UN_M_nest_low$residuals, breaks = 20) 

#lets get those values for our results table 
summary(UN_M_nest_low)
confint(UN_M_nest_low)






#NOW, let's try filtering out all species that have both HIGH and LOW nesting scores of 1 
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



#lets run the model using Phylolm!  

#(have to use lambda, until we figure out a way to make GLS work with binomial linear regression)
UN_M_nest_low_only <- phylolm(Urban~ NestSite_Low + Mass_log, data=NestTraitDat6_only,
                         phy=Nestphy6_only, model="lambda") 

# time to check out the model 
qqnorm(UN_M_nest_low_only$residuals)
qqline(UN_M_nest_low_only$residuals) # what is happening? two separate lines bc of binomial... but is this the correct model check for binomial regression?
#the two lines do not have overlap... maybe this is good? 
hist(UN_M_nest_low_only$residuals, breaks = 20) 

#lets get those values for our results table 
summary(UN_M_nest_low_only)
confint(UN_M_nest_low_only)




######################## UAI and %  Nest Site HIGH ##########################
# 0 = not high
# 1 = HIGH

# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UAIDataUT <- C_Nest_dat2 %>% filter(!is.na(aveUAI)) 
NestData7 <- UAIDataUT %>% filter(!is.na(NestSite_High)) 
length(NestData7$NestSite_High)
#130 species with UAI and NestSite_High

#there is no Log_mass column in this .rds ... let's put one in! 

# add column for log transformed body mass
NestData7 <- NestData7 %>%
  mutate(Mass_log = log(Mass))

colnames(NestData7)

###### add and pair tree

# add rownames to data
row.names(NestData7) <- NestData7$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Nestphydat7 <- treedata(tree_out,NestData7, sort=T)

Nestphy7 <- Nestphydat7$phy
NestTraitDat7 <- as.data.frame(Nestphydat7$data)

str(NestTraitDat7)
length(NestTraitDat7$NestSite_High)
#130


### convert traits of interest to numeric

NestTraitDat7$aveUAI <- as.numeric(NestTraitDat7$aveUAI)
NestTraitDat7$Mass_log <- as.numeric(NestTraitDat7$Mass_log)
NestTraitDat7$NestSite_High <- as.numeric(NestTraitDat7$NestSite_High)


#lets run the model!


UAI_GLS_nest_safety <- gls(aveUAI~ NestSite_High + Mass_log, data = NestTraitDat7, 
                         correlation = corPagel(0.5, phy=Nestphy7,fixed=F, form = ~Species_Jetz), 
                         method = "ML") 
#check out the model
check_model(UAI_GLS_nest_safety) 
qqnorm(resid(UAI_GLS_nest_safety)) 
qqline(resid(UAI_GLS_nest_safety))
hist(resid(UAI_GLS_nest_safety))


summary(UAI_GLS_nest_safety) 
confint(UAI_GLS_nest_safety)



######################## MUTI and %  Nest Site HIGH ##########################
# 0 = not high
# 1 = HIGH

# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
MUTIDataUT <- C_Nest_dat2 %>% filter(!is.na(MUTIscore)) 
NestData8 <- MUTIDataUT %>% filter(!is.na(NestSite_High)) 
length(NestData8$NestSite_High)
#130 species with UAI and NestSite_High

#there is no Log_mass column in this .rds ... let's put one in! 

# add column for log transformed body mass
NestData8 <- NestData8 %>%
  mutate(Mass_log = log(Mass))

colnames(NestData8)

###### add and pair tree

# add rownames to data
row.names(NestData8) <- NestData8$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Nestphydat8 <- treedata(tree_out,NestData8, sort=T)

Nestphy8 <- Nestphydat8$phy
NestTraitDat8 <- as.data.frame(Nestphydat8$data)

str(NestTraitDat8)
length(NestTraitDat8$NestSite_High)
#130


### convert traits of interest to numeric

NestTraitDat8$MUTIscore <- as.numeric(NestTraitDat8$MUTIscore)
NestTraitDat8$Mass_log <- as.numeric(NestTraitDat8$Mass_log)
NestTraitDat8$NestSite_High <- as.numeric(NestTraitDat8$NestSite_High)


#lets run the model!


MUTI_GLS_nest_safety <- gls(MUTIscore~ NestSite_High + Mass_log, data = NestTraitDat8, 
                          correlation = corPagel(0.5, phy=Nestphy8,fixed=F, form = ~Species_Jetz), 
                          method = "ML") 
#check out the model
check_model(MUTI_GLS_nest_safety) 
qqnorm(resid(MUTI_GLS_nest_safety)) 
qqline(resid(MUTI_GLS_nest_safety))
hist(resid(MUTI_GLS_nest_safety))


summary(MUTI_GLS_nest_safety) 
confint(MUTI_GLS_nest_safety)



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



#lets run the model using Phylolm!  

#(have to use lambda, until we figure out a way to make GLS work with binomial linear regression)
UN_M_nest_safety <- phylolm(Urban~ NestSite_High + Mass_log, data=NestTraitDat9,
                         phy=Nestphy6, model="lambda") 

# time to check out the model 
qqnorm(UN_M_nest_safety$residuals)
qqline(UN_M_nest_safety$residuals) # what is happening? two separate lines bc of binomial... but is this the correct model check for binomial regression?
#the two lines do not have overlap... maybe this is good? 
hist(UN_M_nest_safety$residuals, breaks = 20) 

#lets get those values for our results table 
summary(UN_M_nest_safety)
confint(UN_M_nest_safety)




######################## UAI and %  Nest Safety ##########################


# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UAIDataUT <- C_Nest_dat2 %>% filter(!is.na(aveUAI)) 
NestData10 <- UAIDataUT %>% filter(!is.na(nest.safety)) 
length(NestData10$nest.safety)
#766 species with aveUAI and nest.safety

#there is no Log_mass column in this .rds ... let's put one in! 

# add column for log transformed body mass
NestData10 <- NestData10 %>%
  mutate(Mass_log = log(Mass))

colnames(NestData10)

###### add and pair tree

# add rownames to data
row.names(NestData10) <- NestData10$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Nestphydat10 <- treedata(tree_out,NestData10, sort=T)

Nestphy10 <- Nestphydat10$phy
NestTraitDat10 <- as.data.frame(Nestphydat10$data)

str(NestTraitDat10)
length(NestTraitDat10$nest.safety)
#766


### convert traits of interest to numeric

NestTraitDat10$aveUAI <- as.numeric(NestTraitDat10$aveUAI)
NestTraitDat10$Mass_log <- as.numeric(NestTraitDat10$Mass_log)
NestTraitDat10$nest.safety <- as.numeric(NestTraitDat10$nest.safety)



#lets run the model!


UAI_GLS_nest_safety <- gls(aveUAI~ nest.safety + Mass_log, data = NestTraitDat10, 
                            correlation = corPagel(0.5, phy=Nestphy10,fixed=F, form = ~Species_Jetz), 
                            method = "ML") 
#check out the model
check_model(UAI_GLS_nest_safety) 
qqnorm(resid(UAI_GLS_nest_safety)) 
qqline(resid(UAI_GLS_nest_safety))
hist(resid(UAI_GLS_nest_safety))


summary(UAI_GLS_nest_safety) 
confint(UAI_GLS_nest_safety)



######################## MUTI and %  Nest Safety ##########################


# lets first simplify a NEW database by removing records where we don't have an MUTI / brood_value
MUTIDataUT <- C_Nest_dat2 %>% filter(!is.na(MUTIscore)) 
NestData11 <- MUTIDataUT %>% filter(!is.na(nest.safety)) 
length(NestData11$nest.safety)
#127 species with MUTIscore and nest.safety

#there is no Log_mass column in this .rds ... let's put one in! 

# add column for log transformed body mass
NestData11 <- NestData11 %>%
  mutate(Mass_log = log(Mass))

colnames(NestData11)

###### add and pair tree

# add rownames to data
row.names(NestData11) <- NestData11$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Nestphydat11 <- treedata(tree_out,NestData11, sort=T)

Nestphy11 <- Nestphydat11$phy
NestTraitDat11 <- as.data.frame(Nestphydat11$data)

str(NestTraitDat11)
length(NestTraitDat11$nest.safety)
#127


### convert traits of interest to numeric

NestTraitDat11$MUTIscore <- as.numeric(NestTraitDat11$MUTIscore)
NestTraitDat11$Mass_log <- as.numeric(NestTraitDat11$Mass_log)
NestTraitDat11$nest.safety <- as.numeric(NestTraitDat11$nest.safety)



#lets run the model!


MUTI_GLS_nest_safety <- gls(MUTIscore~ nest.safety + Mass_log, data = NestTraitDat11, 
                           correlation = corPagel(0.5, phy=Nestphy11,fixed=F, form = ~Species_Jetz), 
                           method = "ML") 
#check out the model
check_model(MUTI_GLS_nest_safety) 
qqnorm(resid(MUTI_GLS_nest_safety)) 
qqline(resid(MUTI_GLS_nest_safety))
hist(resid(MUTI_GLS_nest_safety))


summary(MUTI_GLS_nest_safety) 
confint(MUTI_GLS_nest_safety)




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



#lets run the model using Phylolm!  

#(have to use lambda, until we figure out a way to make GLS work with binomial linear regression)
UN_M_nest_safety <- phylolm(Urban~ nest.safety + Mass_log, data=NestTraitDat12,
                            phy=Nestphy12, model="lambda") 

# time to check out the model 
qqnorm(UN_M_nest_safety$residuals)
qqline(UN_M_nest_safety$residuals) # what is happening? two separate lines bc of binomial... but is this the correct model check for binomial regression?
#the two lines do not have overlap... maybe this is good? 
hist(UN_M_nest_safety$residuals, breaks = 20) 

#lets get those values for our results table 
summary(UN_M_nest_safety)
confint(UN_M_nest_safety)
