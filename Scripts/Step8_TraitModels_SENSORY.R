##########################################################
######### Coastal Birds Urban Tolerance Project ##########
##########################################################
# Step 8 Phylogenetic Trait Models - SENSORY TRAITS
# Authors: Sarah L. Jennings, Emma M. Garrison
##########################################################
# The objective of this script is to run phylogenetic models using sensory traits
# 2 sensory traits: dim light vision (C.T. ratio) and 


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


# load in "Coastal_Species_Sensory.rds" - this contains coastal species and sensory traits

C_Sensory_dat <- readRDS(here("Outputs", "Coastal_Species_Sensory.rds"))
str(C_Sensory_dat)

C_Sensory_dat2 <- C_Sensory_dat %>%
  mutate(Species_Jetz  = str_replace(Species_Jetz, " ", "_")) %>%
  rowwise() %>%
  mutate(peak_freq_kHz = peak_freq/1000) # find peak vocal frequency in kHz as this will provide easier to interpret model parameters
         
str(C_Sensory_dat2)

C_Sensory_dat2$Urban <- ifelse(C_Sensory_dat2$Urban == "U", 1, 0)
View(C_Sensory_dat2)
colnames(C_Sensory_dat2)


#################### UAI and Dim Light Vision #########################################

# we use C.T ratio as a proxy for dim light vision 

# create a new data frame by removing species with no UAI value or that are missing C.T ratio
CT_UAI_Data <- C_Sensory_dat2 %>% filter(!is.na(aveUAI))  %>% filter(!is.na(C.T)) 
length(CT_UAI_Data$C.T)
#237 species with UAI and CT

###### add and pair tree

# add rownames to data
row.names(CT_UAI_Data) <- CT_UAI_Data$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

CT_UAI_phydat <- treedata(tree_out, CT_UAI_Data, sort=T)

CT_UAI_phy <- CT_UAI_phydat$phy
CT_UAI_dat <- as.data.frame(CT_UAI_phydat$data)

str(CT_UAI_dat)
length(CT_UAI_dat$C.T)
#237

### convert traits of interest to numeric
CT_UAI_dat$aveUAI <- as.numeric(CT_UAI_dat$aveUAI)
CT_UAI_dat$Mass_log <- as.numeric(CT_UAI_dat$Mass_log)
CT_UAI_dat$C.T <- as.numeric(CT_UAI_dat$C.T)


# Run phylogenetic linear model for UAI
# model does not as corPagel with starting point = 0.5 and fixed = F... 
# we need to find a fixed lambda value for the model based on AIC

# create an empty vector to store AIC values
AIC_values <- numeric()

# loop through different values of the parameter for corPagel
for (i in seq(0, 1, by = 0.1)) {
  # Fit the gls model with the current value of i
  model <- gls(aveUAI ~ C.T + Mass_log, 
               data = CT_UAI_dat,
               correlation = corPagel(i, phy = CT_UAI_phy, fixed = TRUE, form = ~Species_Jetz), 
               method = "ML")
  # Extract AIC value and store it in the vector
  AIC_values <- c(AIC_values, AIC(model))
}

# Print AIC values
print(AIC_values)
# 0 = best AIC score 


# run the model with lambda fixed at zero
UAI_GLS_C.T <- gls(aveUAI~ C.T + Mass_log, data = CT_UAI_dat, 
                   correlation = corPagel(0, phy = CT_UAI_phy, fixed=T, form = ~Species_Jetz), 
                   method = "ML") 

# model results and summary
summary(UAI_GLS_C.T) 
confint(UAI_GLS_C.T)

# model diagnostics
check_model(UAI_GLS_C.T) 
qqnorm(resid(UAI_GLS_C.T)) 
qqline(resid(UAI_GLS_C.T)) 
hist(resid(UAI_GLS_C.T)) 

# save model for easy retrieval 
saveRDS(UAI_GLS_C.T, here("Models/UAI", "UAI_GLS_C.T.rds"))


########################### MUTI and Dim Light Vision #############################

# create a new data frame by removing species with no MUTI value or that are missing C.T ratio
CT_MUTI_Data <- C_Sensory_dat2 %>% filter(!is.na(MUTIscore)) %>% 
  filter(!is.na(C.T))
length(CT_MUTI_Data$C.T)
# 69 species with MUTI and CT


# add rownames to data
row.names(CT_MUTI_Data) <- CT_MUTI_Data$Species_Jetz

###### add and pair tree
tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

CT_MUTI_phydat <- treedata(tree_out, CT_MUTI_Data, sort=T)

CT_MUTI_phy <- CT_MUTI_phydat$phy
CT_MUTI_dat <- as.data.frame(CT_MUTI_phydat$data)

str(CT_MUTI_dat)
length(CT_MUTI_dat$C.T)
#69

### convert traits of interest to numeric
CT_MUTI_dat$MUTIscore <- as.numeric(CT_MUTI_dat$MUTIscore)
CT_MUTI_dat$Mass_log <- as.numeric(CT_MUTI_dat$Mass_log)
CT_MUTI_dat$C.T <- as.numeric(CT_MUTI_dat$C.T)

# Run phylogenetic linear model for MUTI
# model does not as corPagel with starting point = 0.5 and fixed = F... 
# we need to find a fixed lambda value for the model based on AIC

# Create an empty vector to store AIC values
AIC_values <- numeric()

# Loop through different values of the parameter for corPagel
for (i in seq(0, 1, by = 0.1)) {
  # Fit the gls model with the current value of i
  model <- gls(MUTIscore ~ C.T + Mass_log, 
               data = CT_MUTI_dat, 
               correlation = corPagel(i, phy = CT_MUTI_phy, fixed = TRUE, form = ~Species_Jetz), 
               method = "ML")
  # Extract AIC value and store it in the vector
  AIC_values <- c(AIC_values, AIC(model))
}

# Print AIC values
print(AIC_values)
# 0.1 = best AIC score 

# run model with lambda fixed at 0.1
MUTI_GLS_C.T <- gls(MUTIscore~ C.T + Mass_log, data = CT_MUTI_dat, 
                   correlation = corPagel(0.1, phy = CT_MUTI_phy, fixed=T, form = ~Species_Jetz), 
                   method = "ML") 

# model summary and results
summary(MUTI_GLS_C.T) 
confint(MUTI_GLS_C.T)

# model diagnostics
check_model(MUTI_GLS_C.T) 
qqnorm(resid(MUTI_GLS_C.T)) 
qqline(resid(MUTI_GLS_C.T)) 
hist(resid(MUTI_GLS_C.T)) 

# save model for easy retrieval 
saveRDS(MUTI_GLS_C.T, here("Models/MUTI", "MUTI_GLS_C.T.rds"))

############################## UN and Dim Light Vision ####################################

# create a new data frame by removing species with no UN value or that are missing C.T ratio
# add rownames to data frame
CT_UN_Data <- C_Sensory_dat2 %>% filter(!is.na(Urban)) %>% filter(!is.na(C.T)) %>%
  column_to_rownames(., var="Species_Jetz")
length(CT_UN_Data$C.T)
# 38 species with UN and CT

###### add and pair tree

# import tree
tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

# combine tree with data
CT_UN_phydat <- treedata(tree_out, CT_UN_Data, sort=T)

CT_UN_phy <- CT_UN_phydat$phy
CT_UN_dat <- as.data.frame(CT_UN_phydat$data)

str(CT_UN_dat)
length(CT_UN_dat$C.T)
#38

### convert traits of interest to numeric
CT_UN_dat$Urban <- as.numeric(CT_UN_dat$Urban)
CT_UN_dat$Mass_log <- as.numeric(CT_UN_dat$Mass_log)
CT_UN_dat$C.T <- as.numeric(CT_UN_dat$C.T)


# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
set.seed(584)
phyglm_UN_CT_scale <- phyloglm(Urban ~ scale(C.T) + scale(Mass_log), 
                               data = CT_UN_dat, 
                               phy = CT_UN_phy,
                               boot = 1000) 

# this successfully converges
summary(phyglm_UN_CT_scale)

# save model for easy retrieval
saveRDS(phyglm_UN_CT_scale, here("Models/UN", "phyglm_UN_CT_scale.rds"))
# load model
phyglm_UN_CT_scale <- readRDS(here("Models/UN", "phyglm_UN_CT_scale.rds"))

# get alpha, t, and half life for the model
(alpha_CT <- phyglm_UN_CT_scale$alpha) # alpha
(phyglm_UN_CT_scale$mean.tip.height) # t (aka mean tip height)
(hl_CT<- log(2)/alpha_CT) # half-life
# compared to t, this is a small Half-Life -> low phylogenetic signal

#############################################################################
#############################################################################

############################# UAI and Peak Frequency ########################

# create a new data frame by removing species with no UAI value or that are missing peak frequency
UAI_PeakFreq <- C_Sensory_dat2 %>% filter(!is.na(aveUAI)) %>% 
  filter(!is.na(peak_freq)) 
length(UAI_PeakFreq$peak_freq_kHz)
#202 species with UAI and peak freq

###### add and pair tree

# add rownames to data
row.names(UAI_PeakFreq) <- UAI_PeakFreq$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_PeakFreq_phydat <- treedata(tree_out, UAI_PeakFreq, sort=T)

UAI_PeakFreq_phy <- UAI_PeakFreq_phydat$phy
UAI_PeakFreq_dat <- as.data.frame(UAI_PeakFreq_phydat$data)

str(UAI_PeakFreq_dat)
length(UAI_PeakFreq_dat$peak_freq_kHz)
#202

### convert traits of interest to numeric
UAI_PeakFreq_dat$aveUAI <- as.numeric(UAI_PeakFreq_dat$aveUAI)
UAI_PeakFreq_dat$Mass_log <- as.numeric(UAI_PeakFreq_dat$Mass_log)
UAI_PeakFreq_dat$peak_freq_kHz <- as.numeric(UAI_PeakFreq_dat$peak_freq_kHz)

# Run phylogenetic linear model for UAI
UAI_GLS_pf <- gls(aveUAI ~ peak_freq_kHz + Mass_log, data = UAI_PeakFreq_dat, 
                    correlation = corPagel(0.5, phy = UAI_PeakFreq_phy, fixed=F, form = ~Species_Jetz), 
                    method = "ML") 

# model summary and results
summary(UAI_GLS_pf) 
confint(UAI_GLS_pf)

# model diagnostics
check_model(UAI_GLS_pf) 
qqnorm(resid(UAI_GLS_pf)) 
qqline(resid(UAI_GLS_pf))
hist(resid(UAI_GLS_pf)) 

# save model for easy retrieval 
saveRDS(UAI_GLS_pf, here("Models/UAI", "UAI_GLS_pf.rds"))

############################# MUTI and Peak Frequency ########################

# create a new data frame by removing species with no MUTI value or that are missing peak frequency
MUTI_PeakFreq <- C_Sensory_dat2 %>% filter(!is.na(MUTIscore)) %>% filter(!is.na(peak_freq)) 
length(MUTI_PeakFreq$peak_freq_kHz)
# 68 species with MUTI and peak frequency

###### add and pair tree

# add rownames to data
row.names(MUTI_PeakFreq) <- MUTI_PeakFreq$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_PeakFreq_phydat <- treedata(tree_out, MUTI_PeakFreq, sort=T)

MUTI_PeakFreq_phy <- MUTI_PeakFreq_phydat$phy
MUTI_PeakFreq_dat <- as.data.frame(MUTI_PeakFreq_phydat$data)

str(MUTI_PeakFreq_dat)
length(MUTI_PeakFreq_dat$peak_freq_kHz)
#68

### convert traits of interest to numeric
MUTI_PeakFreq_dat$MUTIscore <- as.numeric(MUTI_PeakFreq_dat$MUTIscore)
MUTI_PeakFreq_dat$Mass_log <- as.numeric(MUTI_PeakFreq_dat$Mass_log)
MUTI_PeakFreq_dat$peak_freq_kHz <- as.numeric(MUTI_PeakFreq_dat$peak_freq_kHz)


# Run phylogenetic linear model for MUTI
# model does not as corPagel with starting point = 0.5 and fixed = F... 
# we need to find a fixed lambda value for the model based on AIC

# Create an empty vector to store AIC values
AIC_values <- numeric()

# Loop through different values of the parameter for corPagel
for (i in seq(0, 1, by = 0.1)) {
  # Fit the gls model with the current value of i
  model <- gls(MUTIscore ~ peak_freq_kHz + Mass_log, 
               data = MUTI_PeakFreq_dat, 
               correlation = corPagel(i, phy = MUTI_PeakFreq_phy, fixed = TRUE, form = ~Species_Jetz), 
               method = "ML")
  # Extract AIC value and store it in the vector
  AIC_values <- c(AIC_values, AIC(model))
}

# Print AIC values
print(AIC_values)
# 0.1 = best AIC score 

# run the model with lambda fixed at 0.1
MUTI_GLS_pf <- gls(MUTIscore ~ peak_freq_kHz + Mass_log, data = MUTI_PeakFreq_dat, 
                  correlation = corPagel(0.1, phy = MUTI_PeakFreq_phy, fixed=T, form = ~Species_Jetz), 
                  method = "ML") 

# model summary and results
summary(MUTI_GLS_pf) 
confint(MUTI_GLS_pf)

# model diagnostics
check_model(MUTI_GLS_pf) 
qqnorm(resid(MUTI_GLS_pf)) 
qqline(resid(MUTI_GLS_pf))  
hist(resid(MUTI_GLS_pf)) 

# save model for easy retrieval 
saveRDS(MUTI_GLS_pf, here("Models/MUTI", "MUTI_GLS_pf.rds"))


############################# UN and Peak Frequency ########################

# create a new data frame by removing species with no UN value or that are missing peak frequency
UN_PeakFreq <- C_Sensory_dat2 %>% filter(!is.na(Urban))  %>% 
  filter(!is.na(peak_freq)) %>%
  column_to_rownames(., var="Species_Jetz")
length(UN_PeakFreq$peak_freq_kHz)
#129 species with UN and peak frequency

###### add and pair tree

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UN_PeakFreq_phydat <- treedata(tree_out, UN_PeakFreq, sort=T)

UN_PeakFreq_phy <- UN_PeakFreq_phydat$phy
UN_PeakFreq_dat <- as.data.frame(UN_PeakFreq_phydat$data)

str(UN_PeakFreq_dat)
length(UN_PeakFreq_dat$peak_freq_kHz)
#129

### convert traits of interest to numeric
UN_PeakFreq_dat$Urban <- as.numeric(UN_PeakFreq_dat$Urban)
UN_PeakFreq_dat$Mass_log <- as.numeric(UN_PeakFreq_dat$Mass_log)
UN_PeakFreq_dat$peak_freq_kHz <- as.numeric(UN_PeakFreq_dat$peak_freq_kHz)

# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
set.seed(649)
phyglm_UN_pf_scale <- phyloglm( Urban ~ scale(peak_freq_kHz) + scale(Mass_log), 
                                data = UN_PeakFreq_dat, 
                                phy = UN_PeakFreq_phy, 
                                boot = 1000)
# alpha reached upper bound but model does converge
summary(phyglm_UN_pf_scale)

# save model
saveRDS(phyglm_UN_pf_scale, here("Models/UN", "phyglm_UN_pf_scale.rds"))
# load model
phyglm_UN_pf_scale <- readRDS(here("Models/UN", "phyglm_UN_pf_scale.rds"))

# as alpha is at upper bound, also look at regular logistic model
glm_UN_pf_scale <- logistf(Urban ~ scale(peak_freq_kHz) + scale(Mass_log), 
                           data = UN_PeakFreq)
summary(glm_UN_pf_scale)
# coefficient for peak freq is quite different, but we reach the same conclusions

# get alpha, t, and half life for the model
(phyglm_UN_pf_scale$mean.tip.height) # t
(alpha_pfreq <- phyglm_UN_pf_scale$alpha) # alpha
(hl_pfreq<- log(2)/alpha_pfreq) # half life
# low half life relative to mean tip height -> low phylogenetic signal
