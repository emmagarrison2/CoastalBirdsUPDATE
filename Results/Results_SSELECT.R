# The objective of this script is to run phylogenetic linear models for all 
# sexual selection traits. 

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

C_SSelect_dat <- readRDS(here("Outputs", "Coastal_Species_SSelect.rds"))
str(C_SSelect_dat)

C_SSelect_dat2 <- C_SSelect_dat %>%
  mutate(Species_Jetz  = str_replace(Species_Jetz, " ", "_"))
str(C_SSelect_dat2)

C_SSelect_dat2$Urban <- ifelse(C_SSelect_dat2$Urban == "U", 1, 0)
View(C_SSelect_dat2)
colnames(C_SSelect_dat2)



######################## UAI and Dichromatism - BRIGHTNESS ##########################

# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UAIDataUT <- C_SSelect_dat2 %>% filter(!is.na(aveUAI)) 
SSData1 <- UAIDataUT %>% filter(!is.na(Dichrom_bright)) 
length(SSData1$Dichrom_bright)
#201 species with UAI and Dichrom_bright

colnames(SSData1)

###### add and pair tree

# add rownames to data
row.names(SSData1) <- SSData1$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

SSphydat1 <- treedata(tree_out,SSData1, sort=T)

SSphy1 <- SSphydat1$phy
SSTraitDat1 <- as.data.frame(SSphydat1$data)

str(SSTraitDat1)
length(SSTraitDat1$Dichrom_bright)
#201


### convert traits of interest to numeric

SSTraitDat1$aveUAI <- as.numeric(SSTraitDat1$aveUAI)
SSTraitDat1$Mass_log <- as.numeric(SSTraitDat1$Mass_log)
SSTraitDat1$Dichrom_bright <- as.numeric(SSTraitDat1$Dichrom_bright)


#lets run the model!


UAI_GLS_bright <- gls(aveUAI~ Dichrom_bright + Mass_log, data = SSTraitDat1, 
                       correlation = corPagel(0.5, phy=SSphy1,fixed=F, form = ~Species_Jetz), 
                       method = "ML") 
#check out the model
check_model(UAI_GLS_bright) ## low collinearity - which is good! - normality of residuals line does not fall on line, but is in a straight line 
qqnorm(resid(UAI_GLS_bright)) 
qqline(resid(UAI_GLS_bright)) #most points fall on the line 
hist(resid(UAI_GLS_bright)) #roughly normal dist of residuals


summary(UAI_GLS_bright) # strong phylogenetic relationship between traits and response, but no remaining influence of any predictor trait and response
confint(UAI_GLS_bright)



######################## MUTI and Dichromatism - BRIGHTNESS ##########################

# lets first simplify a NEW database by removing records where we don't have an MUTI / brood_value
MUTIDataUT <- C_SSelect_dat2 %>% filter(!is.na(MUTIscore)) 
SSData2 <- MUTIDataUT %>% filter(!is.na(Dichrom_bright)) 
length(SSData2$Dichrom_bright)
#67 species with MUTI and Dichrom_bright

colnames(SSData2)

###### add and pair tree

# add rownames to data
row.names(SSData2) <- SSData2$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

SSphydat2 <- treedata(tree_out,SSData2, sort=T)

SSphy2 <- SSphydat2$phy
SSTraitDat2 <- as.data.frame(SSphydat2$data)

str(SSTraitDat2)
length(SSTraitDat2$Dichrom_bright)
#67


### convert traits of interest to numeric

SSTraitDat2$MUTIscore <- as.numeric(SSTraitDat2$MUTIscore)
SSTraitDat2$Mass_log <- as.numeric(SSTraitDat2$Mass_log)
SSTraitDat2$Dichrom_bright <- as.numeric(SSTraitDat2$Dichrom_bright)


########### since model not working as corPagel starting point = 0.5 and fixed = F... 
# let's find a lambda value to fix in the model 


# Create an empty vector to store AIC values
AIC_values <- numeric()

# Loop through different values of the parameter for corPagel
for (i in seq(0, 1, by = 0.1)) {
  # Fit the gls model with the current value of i
  model <- gls(MUTIscore ~ Dichrom_bright + Mass_log, 
               data = SSTraitDat2, 
               correlation = corPagel(i, phy = SSphy2, fixed = TRUE, form = ~Species_Jetz), 
               method = "ML")
  # Extract AIC value and store it in the vector
  AIC_values <- c(AIC_values, AIC(model))
}

# Print AIC values
print(AIC_values)
#0.4 = best AIC score 




#lets run the model!


MUTI_GLS_bright <- gls(MUTIscore~ Dichrom_bright + Mass_log, data = SSTraitDat2, 
                      correlation = corPagel(0.4, phy=SSphy2,fixed=T, form = ~Species_Jetz), 
                      method = "ML") 
#check out the model
check_model(MUTI_GLS_bright) ## low collinearity - which is good! - normality of residuals line does not fall on line, but is in a straight line 
qqnorm(resid(MUTI_GLS_bright)) 
qqline(resid(MUTI_GLS_bright)) #most points fall on the line 
hist(resid(MUTI_GLS_bright)) #roughly normal dist of residuals


summary(MUTI_GLS_bright) # strong phylogenetic relationship between traits and response, but no remaining influence of any predictor trait and response
confint(MUTI_GLS_bright)




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




#lets run the model using Phylolm!  
#(have to use lambda, until we figure out a way to make GLS work with binomial linear regression)
UN_M_bright <- phylolm(Urban~ Dichrom_bright + Mass_log, data=SSTraitDat3,
                   phy=SSphy3, model="lambda") 

# time to check out the model 
qqnorm(UN_M_bright$residuals)
qqline(UN_M_bright$residuals) # what is happening? two separate lines bc of binomial... but is this the correct model check for binomial regression?
#the two lines do not have overlap... maybe this is good? 
hist(UN_M_bright$residuals, breaks = 12) 

#lets get those values for our results table 
summary(UN_M_bright)
confint(UN_M_bright)



######################## UAI and Dichromatism - HUE ##########################

# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UAIDataUT <- C_SSelect_dat2 %>% filter(!is.na(aveUAI)) 
SSData4 <- UAIDataUT %>% filter(!is.na(Dichrom_hue)) 
length(SSData4$Dichrom_hue)
#201 species with UAI and Dichrom_hue

colnames(SSData4)

###### add and pair tree

# add rownames to data
row.names(SSData4) <- SSData4$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

SSphydat4 <- treedata(tree_out,SSData4, sort=T)

SSphy4 <- SSphydat4$phy
SSTraitDat4 <- as.data.frame(SSphydat4$data)

str(SSTraitDat4)
length(SSTraitDat4$Dichrom_hue)
#201


### convert traits of interest to numeric

SSTraitDat4$aveUAI <- as.numeric(SSTraitDat4$aveUAI)
SSTraitDat4$Mass_log <- as.numeric(SSTraitDat4$Mass_log)
SSTraitDat4$Dichrom_hue <- as.numeric(SSTraitDat4$Dichrom_hue)


#lets run the model!


UAI_GLS_hue <- gls(aveUAI~ Dichrom_hue + Mass_log, data = SSTraitDat4, 
                      correlation = corPagel(0.5, phy=SSphy4,fixed=F, form = ~Species_Jetz), 
                      method = "ML") 
#check out the model
check_model(UAI_GLS_hue) ## low collinearity - which is good! - normality of residuals line does not fall on line, but is in a straight line 
qqnorm(resid(UAI_GLS_hue)) 
qqline(resid(UAI_GLS_hue)) #most points fall on the line 
hist(resid(UAI_GLS_hue)) #roughly normal dist of residuals


summary(UAI_GLS_hue) # strong phylogenetic relationship between traits and response, but no remaining influence of any predictor trait and response
confint(UAI_GLS_hue)





######################## MUTI and Dichromatism - HUE ##########################

# lets first simplify a NEW database by removing records where we don't have an MUTI / brood_value
MUTIDataUT <- C_SSelect_dat2 %>% filter(!is.na(MUTIscore)) 
SSData5 <- MUTIDataUT %>% filter(!is.na(Dichrom_hue)) 
length(SSData5$Dichrom_hue)
#67 species with MUTI and Dichrom_hue

colnames(SSData5)

###### add and pair tree

# add rownames to data
row.names(SSData5) <- SSData5$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

SSphydat5 <- treedata(tree_out,SSData5, sort=T)

SSphy5 <- SSphydat5$phy
SSTraitDat5 <- as.data.frame(SSphydat5$data)

str(SSTraitDat5)
length(SSTraitDat5$Dichrom_hue)
#67


### convert traits of interest to numeric

SSTraitDat5$MUTIscore <- as.numeric(SSTraitDat5$MUTIscore)
SSTraitDat5$Mass_log <- as.numeric(SSTraitDat5$Mass_log)
SSTraitDat5$Dichrom_hue <- as.numeric(SSTraitDat5$Dichrom_hue)


#lets run the model!


MUTI_GLS_hue <- gls(MUTIscore~ Dichrom_hue + Mass_log, data = SSTraitDat5, 
                   correlation = corPagel(0.5, phy=SSphy5,fixed=F, form = ~Species_Jetz), 
                   method = "ML") 
#check out the model
check_model(MUTI_GLS_hue) 
qqnorm(resid(MUTI_GLS_hue)) 
qqline(resid(MUTI_GLS_hue)) 
hist(resid(MUTI_GLS_hue)) 


summary(MUTI_GLS_hue)
confint(MUTI_GLS_hue)





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



#lets run the model using Phylolm!  

#(have to use lambda, until we figure out a way to make GLS work with binomial linear regression)
UN_M_hue <- phylolm(Urban~ Dichrom_hue + Mass_log, data=SSTraitDat6,
                       phy=SSphy6, model="lambda") 

# time to check out the model 
qqnorm(UN_M_hue$residuals)
qqline(UN_M_hue$residuals) # what is happening? two separate lines bc of binomial... but is this the correct model check for binomial regression?
#the two lines do not have overlap... maybe this is good? 
hist(UN_M_hue$residuals, breaks = 12) 

#lets get those values for our results table 
summary(UN_M_hue)
confint(UN_M_hue)




######################## UAI and Dichromatism - HUE ##########################

# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UAIDataUT <- C_SSelect_dat2 %>% filter(!is.na(aveUAI)) 
SSData7 <- UAIDataUT %>% filter(!is.na(sex.sel.m)) 
length(SSData7$sex.sel.m)
#766 species with UAI and sex.sel.m

colnames(SSData7)

###### add and pair tree

# add rownames to data
row.names(SSData7) <- SSData7$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

SSphydat7 <- treedata(tree_out,SSData7, sort=T)

SSphy7 <- SSphydat7$phy
SSTraitDat7 <- as.data.frame(SSphydat7$data)

str(SSTraitDat7)
length(SSTraitDat7$sex.sel.m)
#766


### convert traits of interest to numeric

SSTraitDat7$aveUAI <- as.numeric(SSTraitDat7$aveUAI)
SSTraitDat7$Mass_log <- as.numeric(SSTraitDat7$Mass_log)
SSTraitDat7$sex.sel.m <- as.numeric(SSTraitDat7$sex.sel.m)


#lets run the model!


UAI_GLS_ssm <- gls(aveUAI~ sex.sel.m + Mass_log, data = SSTraitDat7, 
                   correlation = corPagel(0.5, phy=SSphy7,fixed=F, form = ~Species_Jetz), 
                   method = "ML") 
#check out the model
check_model(UAI_GLS_ssm)
qqnorm(resid(UAI_GLS_ssm)) 
qqline(resid(UAI_GLS_ssm))
hist(resid(UAI_GLS_ssm)) 


summary(UAI_GLS_ssm) 
confint(UAI_GLS_ssm)



######################## MUTI and Dichromatism - HUE ##########################

# lets first simplify a NEW database by removing records where we don't have an MUTI / brood_value
MUTIDataUT <- C_SSelect_dat2 %>% filter(!is.na(MUTIscore)) 
SSData8 <- MUTIDataUT %>% filter(!is.na(sex.sel.m)) 
length(SSData8$sex.sel.m)
#127 species with MUTI and sex.sel.m

colnames(SSData8)

###### add and pair tree

# add rownames to data
row.names(SSData8) <- SSData8$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

SSphydat8 <- treedata(tree_out,SSData8, sort=T)

SSphy8 <- SSphydat8$phy
SSTraitDat8 <- as.data.frame(SSphydat8$data)

str(SSTraitDat8)
length(SSTraitDat8$sex.sel.m)
#127


### convert traits of interest to numeric

SSTraitDat8$MUTIscore <- as.numeric(SSTraitDat8$MUTIscore)
SSTraitDat8$Mass_log <- as.numeric(SSTraitDat8$Mass_log)
SSTraitDat8$sex.sel.m <- as.numeric(SSTraitDat8$sex.sel.m)


#lets run the model!


MUTI_GLS_ssm <- gls(MUTIscore~ sex.sel.m + Mass_log, data = SSTraitDat8, 
                   correlation = corPagel(0.5, phy=SSphy8,fixed=F, form = ~Species_Jetz), 
                   method = "ML") 
#check out the model
check_model(MUTI_GLS_ssm)
qqnorm(resid(MUTI_GLS_ssm)) 
qqline(resid(MUTI_GLS_ssm))
hist(resid(MUTI_GLS_ssm)) 


summary(MUTI_GLS_ssm) 
confint(MUTI_GLS_ssm)




######################## UN and Dichromatism - HUE ##########################

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



#lets run the model using Phylolm!  

#(have to use lambda, until we figure out a way to make GLS work with binomial linear regression)
UN_M_ssm <- phylolm(Urban~ sex.sel.m + Mass_log, data=SSTraitDat9,
                    phy=SSphy9, model="lambda") 

# time to check out the model 
qqnorm(UN_M_ssm$residuals)
qqline(UN_M_ssm$residuals)
hist(UN_M_ssm$residuals, breaks = 12) 

#lets get those values for our results table 
summary(UN_M_ssm)
confint(UN_M_ssm)




######################## UAI and Sexual Selection - Female ##########################

# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UAIDataUT <- C_SSelect_dat2 %>% filter(!is.na(aveUAI)) 
SSData10 <- UAIDataUT %>% filter(!is.na(sex.sel.f)) 
length(SSData10$sex.sel.f)
#766 species with UAI and sex.sel.f

colnames(SSData10)

###### add and pair tree

# add rownames to data
row.names(SSData10) <- SSData10$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

SSphydat10 <- treedata(tree_out,SSData10, sort=T)

SSphy10 <- SSphydat10$phy
SSTraitDat10 <- as.data.frame(SSphydat10$data)

str(SSTraitDat10)
length(SSTraitDat10$sex.sel.f)
#766


### convert traits of interest to numeric

SSTraitDat10$aveUAI <- as.numeric(SSTraitDat10$aveUAI)
SSTraitDat10$Mass_log <- as.numeric(SSTraitDat10$Mass_log)
SSTraitDat10$sex.sel.f <- as.numeric(SSTraitDat10$sex.sel.f)



#lets run the model!


UAI_GLS_ssf <- gls(aveUAI~ sex.sel.f + Mass_log, data = SSTraitDat10, 
                    correlation = corPagel(0.5, phy=SSphy10,fixed=F, form = ~Species_Jetz), 
                    method = "ML") 
#check out the model
check_model(UAI_GLS_ssf)
qqnorm(resid(UAI_GLS_ssf)) 
qqline(resid(UAI_GLS_ssf))
hist(resid(UAI_GLS_ssf)) 


summary(UAI_GLS_ssf) 
confint(UAI_GLS_ssf)




######################## MUTI and Sexual Selection - Female ##########################

# lets first simplify a NEW database by removing records where we don't have an MUTI / brood_value
MUTIDataUT <- C_SSelect_dat2 %>% filter(!is.na(MUTIscore)) 
SSData11 <- MUTIDataUT %>% filter(!is.na(sex.sel.f)) 
length(SSData11$sex.sel.f)
#127 species with MUTI and sex.sel.f

colnames(SSData11)

###### add and pair tree

# add rownames to data
row.names(SSData11) <- SSData11$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

SSphydat11 <- treedata(tree_out,SSData11, sort=T)

SSphy11 <- SSphydat11$phy
SSTraitDat11 <- as.data.frame(SSphydat11$data)

str(SSTraitDat11)
length(SSTraitDat11$sex.sel.f)
#127


### convert traits of interest to numeric

SSTraitDat11$MUTIscore <- as.numeric(SSTraitDat11$MUTIscore)
SSTraitDat11$Mass_log <- as.numeric(SSTraitDat11$Mass_log)
SSTraitDat11$sex.sel.f <- as.numeric(SSTraitDat11$sex.sel.f)



#lets run the model!


MUTI_GLS_ssf <- gls(MUTIscore~ sex.sel.f + Mass_log, data = SSTraitDat11, 
                   correlation = corPagel(0.5, phy=SSphy11,fixed=F, form = ~Species_Jetz), 
                   method = "ML") 
#check out the model
check_model(MUTI_GLS_ssf)
qqnorm(resid(MUTI_GLS_ssf)) 
qqline(resid(MUTI_GLS_ssf))
hist(resid(MUTI_GLS_ssf)) 


summary(MUTI_GLS_ssf) 
confint(MUTI_GLS_ssf)




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



#lets run the model using Phylolm!  

#(have to use lambda, until we figure out a way to make GLS work with binomial linear regression)
UN_M_ssf <- phylolm(Urban~ sex.sel.f + Mass_log, data=SSTraitDat12,
                    phy=SSphy12, model="lambda") 

# time to check out the model 
qqnorm(UN_M_ssf$residuals)
qqline(UN_M_ssf$residuals)
hist(UN_M_ssf$residuals, breaks = 12) 

#lets get those values for our results table 
summary(UN_M_ssf)
confint(UN_M_ssf)
