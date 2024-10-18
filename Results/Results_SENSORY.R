######Results - Sensory Data Traits 

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
if(!require(phylolm)){
  install.packages("phylolm")
  require(phylolm)
}
library(phylolm)

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


#############################################################
###########
#C.T and UAI YAYAYAY!

# lets first simplify a NEW database by removing records where we don't have an UAI value from Neate-Clegg
NeateTraitDataUT <- C_Sensory_dat2 %>% filter(!is.na(aveUAI)) 
NeateTraitData1 <- NeateTraitDataUT %>% filter(!is.na(C.T)) 
length(NeateTraitData1$C.T)
#237 species with UAI and CT

###### add and pair tree

# add rownames to data
row.names(NeateTraitData1) <- NeateTraitData1$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Neatephydat1 <- treedata(tree_out,NeateTraitData1, sort=T)

Neatephy1 <- Neatephydat1$phy
NeateTraitDat1 <- as.data.frame(Neatephydat1$data)

str(NeateTraitDat1)
length(NeateTraitDat1$C.T)
#237

### convert traits of interest to numeric

NeateTraitDat1$aveUAI <- as.numeric(NeateTraitDat1$aveUAI)
NeateTraitDat1$Mass_log <- as.numeric(NeateTraitDat1$Mass_log)
NeateTraitDat1$C.T <- as.numeric(NeateTraitDat1$C.T)

########### since model not working as corPagel starting point = 0.5 and fixed = F... 
# let's find a
# lambda value to fix in the model 


# Create an empty vector to store AIC values
AIC_values <- numeric()

# Loop through different values of the parameter for corPagel
for (i in seq(0, 1, by = 0.1)) {
  # Fit the gls model with the current value of i
  model <- gls(aveUAI ~ C.T + Mass_log, 
               data = NeateTraitDat1,
               correlation = corPagel(i, phy = Neatephy1, fixed = TRUE, form = ~Species_Jetz), 
               method = "ML")
  # Extract AIC value and store it in the vector
  AIC_values <- c(AIC_values, AIC(model))
}

# Print AIC values
print(AIC_values)
#0 = best AIC score 




#lets run the model!


UAI_GLS_C.T <- gls(aveUAI~ C.T + Mass_log, data = NeateTraitDat1, 
                   correlation = corPagel(0, phy=Neatephy1,fixed=T, form = ~Species_Jetz), 
                   method = "ML") 
#check out the model
check_model(UAI_GLS_C.T) ## low collinearity - which is good! - normality of residuals line does not fall on line, but is in a straight line 
qqnorm(resid(UAI_GLS_C.T)) 
qqline(resid(UAI_GLS_C.T)) #most points fall on the line 
hist(resid(UAI_GLS_C.T)) #roughly normal dist of residuals


summary(UAI_GLS_C.T) # strong phylogenetic relationship between traits and response, but no remaining influence of any predictor trait and response
confint(UAI_GLS_C.T)


########################### C.T and MUTI #############################

# lets first simplify a NEW database by removing records where we don't have an UAI value from Neate-Clegg
TraitDataUTSensory <- C_Sensory_dat2 %>% filter(!is.na(MUTIscore)) 
SensoryTraitData1 <- TraitDataUTSensory %>% filter(!is.na(C.T)) 
length(SensoryTraitData1$C.T)
#69 species with UAI and CT

###### add and pair tree

# add rownames to data
row.names(SensoryTraitData1) <- SensoryTraitData1$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Sensoryphydat1 <- treedata(tree_out,SensoryTraitData1, sort=T)

Sensoryphy1 <- Sensoryphydat1$phy
SensoryTraitDat1 <- as.data.frame(Sensoryphydat1$data)

str(SensoryTraitDat1)
length(SensoryTraitDat1$C.T)
#69

### convert traits of interest to numeric

SensoryTraitDat1$MUTIscore <- as.numeric(SensoryTraitDat1$MUTIscore)
SensoryTraitDat1$Mass_log <- as.numeric(SensoryTraitDat1$Mass_log)
SensoryTraitDat1$C.T <- as.numeric(SensoryTraitDat1$C.T)

########### since model not working as corPagel starting point = 0.5 and fixed = F... 
# let's find a
# lambda value to fix in the model 


# Create an empty vector to store AIC values
AIC_values <- numeric()

# Loop through different values of the parameter for corPagel
for (i in seq(0, 1, by = 0.1)) {
  # Fit the gls model with the current value of i
  model <- gls(MUTIscore ~ C.T + Mass_log, 
               data = SensoryTraitDat1, 
               correlation = corPagel(i, phy = Sensoryphy1, fixed = TRUE, form = ~Species_Jetz), 
               method = "ML")
  # Extract AIC value and store it in the vector
  AIC_values <- c(AIC_values, AIC(model))
}

# Print AIC values
print(AIC_values)
#0.1 = best AIC score 


#lets run the model!


MUTI_GLS_C.T <- gls(MUTIscore~ C.T + Mass_log, data = SensoryTraitDat1, 
                   correlation = corPagel(0.1, phy=Sensoryphy1,fixed=T, form = ~Species_Jetz), 
                   method = "ML") 
#check out the model
check_model(MUTI_GLS_C.T) ## low collinearity - which is good! - normality of residuals line does not fall on line, but is in a straight line 
qqnorm(resid(MUTI_GLS_C.T)) 
qqline(resid(MUTI_GLS_C.T)) #most points fall on the line 
hist(resid(MUTI_GLS_C.T)) #roughly normal dist of residuals


summary(MUTI_GLS_C.T) # strong phylogenetic relationship between traits and response, but no remaining influence of any predictor trait and response
confint(MUTI_GLS_C.T)

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

########### since model not working as corPagel starting point = 0.5 and fixed = F... 
# let's find a
# lambda value to fix in the model 


# Create an empty vector to store AIC values
AIC_values <- numeric()

# Print AIC values
print(AIC_values)
#0.1 = best AIC score 


#lets run the model!(have to use lambda, until we figure out a way to make GLS work with binomial linear regression)
UN_M_C.T <- phylolm(Urban~ C.T + Mass_log, data=SensoryTraitDat2,
                    phy=Sensoryphy2, model="lambda") 

# time to check out the model 
qqnorm(UN_M_C.T$residuals)
qqline(UN_M_C.T$residuals) # what is happening? two separate lines bc of binomial... but is this the correct model check for binomial regression?
#the two lines do not have overlap... maybe this is good? 
hist(UN_M_C.T$residuals, breaks = 12) 

#lets get those values for our results table 
summary(UN_M_C.T)
confint(UN_M_C.T)






############################# Peak Frequency and UAI #################################

# lets first simplify a NEW database by removing records where we don't have an UAI value from Neate-Clegg
UAISensory <- C_Sensory_dat2 %>% filter(!is.na(aveUAI)) 
SensoryTraitData3 <- UAISensory %>% filter(!is.na(peak_freq)) 
length(SensoryTraitData3$peak_freq)
#202 species with UAI and CT

###### add and pair tree

# add rownames to data
row.names(SensoryTraitData3) <- SensoryTraitData3$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Sensoryphydat3 <- treedata(tree_out,SensoryTraitData3, sort=T)

Sensoryphy3 <- Sensoryphydat3$phy
SensoryTraitDat3 <- as.data.frame(Sensoryphydat3$data)

str(SensoryTraitDat3)
length(SensoryTraitDat3$peak_freq)
#202

### convert traits of interest to numeric

SensoryTraitDat3$aveUAI <- as.numeric(SensoryTraitDat3$aveUAI)
SensoryTraitDat3$Mass_log <- as.numeric(SensoryTraitDat3$Mass_log)
SensoryTraitDat3$peak_freq <- as.numeric(SensoryTraitDat3$peak_freq)

########### since model not working as corPagel starting point = 0.5 and fixed = F... 
# let's find a lambda value to fix in the model 


# Create an empty vector to store AIC values
AIC_values <- numeric()

# Loop through different values of the parameter for corPagel
for (i in seq(0, 1, by = 0.1)) {
  # Fit the gls model with the current value of i
  model <- gls(aveUAI ~ peak_freq + Mass_log, 
               data = SensoryTraitDat3, 
               correlation = corPagel(i, phy = Sensoryphy3, fixed = TRUE, form = ~Species_Jetz), 
               method = "ML")
  # Extract AIC value and store it in the vector
  AIC_values <- c(AIC_values, AIC(model))
}

# Print AIC values
print(AIC_values)
#0.7 = best AIC score 


#lets run the model!using GLS 

UAI_GLS_pf <- gls(aveUAI~ peak_freq + Mass_log, data = SensoryTraitDat3, 
                    correlation = corPagel(0.5, phy=Sensoryphy3,fixed=F, form = ~Species_Jetz), 
                    method = "ML") 

#check out the model
check_model(UAI_GLS_pf) ## low collinearity - which is good! - normality of residuals line does not fall on line, but is in a straight line 
qqnorm(resid(UAI_GLS_pf)) 
qqline(resid(UAI_GLS_pf)) #most points fall on the line 
hist(resid(UAI_GLS_pf)) #roughly normal dist of residuals


summary(UAI_GLS_pf) # strong phylogenetic relationship between traits and response, but no remaining influence of any predictor trait and response
confint(UAI_GLS_pf)




############################# Peak Frequency and MUTI #################################

# lets first simplify a NEW database by removing records where we don't have an UAI value from Neate-Clegg
MUTISensory <- C_Sensory_dat2 %>% filter(!is.na(MUTIscore)) 
SensoryTraitData4 <- MUTISensory %>% filter(!is.na(peak_freq)) 
length(SensoryTraitData4$peak_freq)
#202 species with UAI and peak frequency

###### add and pair tree

# add rownames to data
row.names(SensoryTraitData4) <- SensoryTraitData4$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Sensoryphydat4 <- treedata(tree_out,SensoryTraitData4, sort=T)

Sensoryphy4 <- Sensoryphydat4$phy
SensoryTraitDat4 <- as.data.frame(Sensoryphydat4$data)

str(SensoryTraitDat4)
length(SensoryTraitDat4$peak_freq)
#68

### convert traits of interest to numeric

SensoryTraitDat4$MUTIscore <- as.numeric(SensoryTraitDat4$MUTIscore)
SensoryTraitDat4$Mass_log <- as.numeric(SensoryTraitDat4$Mass_log)
SensoryTraitDat4$peak_freq <- as.numeric(SensoryTraitDat4$peak_freq)

########### since model not working as corPagel starting point = 0.5 and fixed = F... 
# let's find a lambda value to fix in the model 


# Create an empty vector to store AIC values
AIC_values <- numeric()

# Loop through different values of the parameter for corPagel
for (i in seq(0, 1, by = 0.1)) {
  # Fit the gls model with the current value of i
  model <- gls(MUTIscore ~ peak_freq + Mass_log, 
               data = SensoryTraitDat4, 
               correlation = corPagel(i, phy = Sensoryphy4, fixed = TRUE, form = ~Species_Jetz), 
               method = "ML")
  # Extract AIC value and store it in the vector
  AIC_values <- c(AIC_values, AIC(model))
}

# Print AIC values
print(AIC_values)
#0.1 = best AIC score 


#lets run the model!using GLS 

MUTI_GLS_pf <- gls(MUTIscore~ peak_freq + Mass_log, data = SensoryTraitDat4, 
                  correlation = corPagel(0.1, phy=Sensoryphy4,fixed=T, form = ~Species_Jetz), 
                  method = "ML") 

#check out the model
check_model(MUTI_GLS_pf) ## low collinearity - which is good! - normality of residuals line does not fall on line, but is in a straight line 
qqnorm(resid(MUTI_GLS_pf)) 
qqline(resid(MUTI_GLS_pf)) #most points fall on the line 
hist(resid(MUTI_GLS_pf)) #roughly normal dist of residuals


summary(MUTI_GLS_pf) # strong phylogenetic relationship between traits and response, but no remaining influence of any predictor trait and response
confint(MUTI_GLS_pf)




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


#lets run the model using Phylolm!  
#(have to use lambda, until we figure out a way to make GLS work with binomial linear regression)
UN_M_pf <- phylolm(Urban~ peak_freq + Mass_log, data=SensoryTraitDat5,
                    phy=Sensoryphy5, model="lambda") 

# time to check out the model 
qqnorm(UN_M_pf$residuals)
qqline(UN_M_pf$residuals) # what is happening? two separate lines bc of binomial... but is this the correct model check for binomial regression?
#the two lines do not have overlap... maybe this is good? 
hist(UN_M_pf$residuals, breaks = 12) 

#lets get those values for our results table 
summary(UN_M_pf)
confint(UN_M_pf)


