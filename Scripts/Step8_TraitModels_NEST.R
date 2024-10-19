##########################################################
######### Coastal Birds Urban Tolerance Project ##########
##########################################################
# Step 8 Phylogenetic Trait Models - NEST TRAITS
# Authors: Sarah L. Jennings, Emma M. Garrison
##########################################################
# The objective of this script is to run phylogenetic models using nest traits
# 4 life history traits: nest strategy, nest site, and nest safety


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

###################### Prep

# load in "Coastal_Species_Nest.rds"
C_Nest_dat <- readRDS(here("Outputs", "Coastal_Species_Nest.rds"))
str(C_Nest_dat)
nrow(C_Nest_dat)
#807

C_Nest_dat2 <- C_Nest_dat %>%
  mutate(Species_Jetz  = str_replace(Species_Jetz, " ", "_"))
str(C_Nest_dat2)


C_Nest_dat2$Urban <- ifelse(C_Nest_dat2$Urban == "U", 1, 0)
View(C_Nest_dat2)
colnames(C_Nest_dat2)


######################## UAI and Nest Strategy ##########################
# 0 = enclosed
# 1 = open

# create a new data frame that contains only species with both UAI and nest strategy
UAI_NestStr <- C_Nest_dat2 %>% filter(!is.na(aveUAI)) %>% 
  filter(!is.na(NestStr)) %>% as.data.frame()
length(UAI_NestStr$NestStr)
# 833 species with UAI and NestStr

###### add and pair tree

# add rownames to data
row.names(UAI_NestStr) <- UAI_NestStr$Species_Jetz

tree_out <- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_NestStr_phydat <- treedata(tree_out, UAI_NestStr, sort=T)

UAI_NestStr_phy <- UAI_NestStr_phydat$phy
UAI_NestStr_dat <- as.data.frame(UAI_NestStr_phydat$data)

str(UAI_NestStr_dat)
length(UAI_NestStr_dat$NestStr)
View(UAI_NestStr_dat)
# 833

# convert traits of interest to numeric
UAI_NestStr_dat$aveUAI <- as.numeric(UAI_NestStr_dat$aveUAI)
UAI_NestStr_dat$Mass_log <- as.numeric(UAI_NestStr_dat$Mass_log)
UAI_NestStr_dat$NestStr <- as.numeric(UAI_NestStr_dat$NestStr)

# run phylogenetic linear model
UAI_GLS_neststr <- gls(aveUAI~ NestStr + Mass_log, data = UAI_NestStr_dat, 
                  correlation = corPagel(0.5, phy= UAI_NestStr_phy, fixed=F, form = ~Species_Jetz), 
                  method = "ML") 

# model summary and results
summary(UAI_GLS_neststr) 
confint(UAI_GLS_neststr)

# model diagnostics
check_model(UAI_GLS_neststr) 
qqnorm(resid(UAI_GLS_neststr)) 
qqline(resid(UAI_GLS_neststr)) 
hist(resid(UAI_GLS_neststr)) 

# save model for easy retrieval 
saveRDS(UAI_GLS_neststr, here("Models/UAI", "UAI_GLS_neststr.rds"))


######################## MUTI and Nest Strategy ##########################
# 0 = enclosed
# 1 = open

# create a new data frame that contains only species with both MUTI and nest strategy
MUTI_NestStr <- C_Nest_dat2 %>% filter(!is.na(MUTIscore)) %>% 
  filter(!is.na(NestStr)) %>% as.data.frame()
length(MUTI_NestStr$NestStr)
# 117 species with MUTIscore and NestStr

###### add and pair tree

# add rownames to data
row.names(MUTI_NestStr) <- MUTI_NestStr$Species_Jetz

tree_out <- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_NestStr_phydat <- treedata(tree_out, MUTI_NestStr, sort=T)

MUTI_NestStr_phy <- MUTI_NestStr_phydat$phy
MUTI_NestStr_dat <- as.data.frame(MUTI_NestStr_phydat$data)

str(MUTI_NestStr_dat)
length(MUTI_NestStr_dat$NestStr)
#117


# convert traits of interest to numeric
MUTI_NestStr_dat$MUTIscore <- as.numeric(MUTI_NestStr_dat$MUTIscore)
MUTI_NestStr_dat$Mass_log <- as.numeric(MUTI_NestStr_dat$Mass_log)
MUTI_NestStr_dat$NestStr <- as.numeric(MUTI_NestStr_dat$NestStr)


# run phylogenetic linear model
MUTI_GLS_neststr <- gls(MUTIscore~ NestStr + Mass_log, data = MUTI_NestStr_dat, 
                       correlation = corPagel(0.5, phy= MUTI_NestStr_phy, fixed = F, form = ~Species_Jetz), 
                       method = "ML") 

# model summary and results
summary(MUTI_GLS_neststr) 
confint(MUTI_GLS_neststr)

# model diagnostic
check_model(MUTI_GLS_neststr) 
qqnorm(resid(MUTI_GLS_neststr)) 
qqline(resid(MUTI_GLS_neststr)) 
hist(resid(MUTI_GLS_neststr)) 


# save model for easy retrieval 
saveRDS(MUTI_GLS_neststr, here("Models/MUTI", "MUTI_GLS_neststr.rds"))


######################## UN and Nest Strategy ##########################
# 0 = enclosed
# 1 = open

# create a new data frame that contains only species with both UN and nest strategy
UN_NestStr <- C_Nest_dat2 %>% filter(!is.na(Urban)) %>%
  filter(!is.na(NestStr)) %>% column_to_rownames(., var="Species_Jetz")
length(UN_NestStr$NestStr)
# 122 species with UN and NestStr

UN_NestStr %>% filter(!is.na(NestStr)) %>% 
  filter(!is.na(Urban)) %>% 
  group_by(NestStr) %>% count()

# do not run --> only 8 species for UN with "enclosed" = 0 Nest Strategy 
# too small of an N for this category 

##########################################################################
##########################################################################

######################## UAI and Nest Site LOW ##########################
# 0 = not low
# 1 = low

# create a new data frame that contains only species with both UAI and nest site low
UAI_NestLow <- C_Nest_dat2 %>% filter(!is.na(aveUAI))  %>% 
  filter(!is.na(NestSite_Low)) %>% as.data.frame()
length(UAI_NestLow$NestSite_Low)
#796 species with UAI and NestSite_Low

###### add and pair tree

# add rownames to data
row.names(UAI_NestLow) <- UAI_NestLow$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_NestLow_phydat <- treedata(tree_out, UAI_NestLow, sort=T)

UAI_NestLow_phy <- UAI_NestLow_phydat$phy
UAI_NestLow_dat <- as.data.frame(UAI_NestLow_phydat$data)

str(UAI_NestLow_dat)
length(UAI_NestLow_dat$NestSite_Low)
#796

# convert traits of interest to numeric
UAI_NestLow_dat$aveUAI <- as.numeric(UAI_NestLow_dat$aveUAI)
UAI_NestLow_dat$Mass_log <- as.numeric(UAI_NestLow_dat$Mass_log)
UAI_NestLow_dat$NestSite_Low <- as.numeric(UAI_NestLow_dat$NestSite_Low)

# run a phylogenetic linear model
UAI_GLS_nest_low <- gls(aveUAI~ NestSite_Low + Mass_log, data = UAI_NestLow_dat, 
                       correlation = corPagel(0.5, phy = UAI_NestLow_phy, fixed = F, form = ~Species_Jetz), 
                       method = "ML") 

# model summary and results
summary(UAI_GLS_nest_low) 
confint(UAI_GLS_nest_low)

# model diagnostics
check_model(UAI_GLS_nest_low) 
qqnorm(resid(UAI_GLS_nest_low)) 
qqline(resid(UAI_GLS_nest_low))
hist(resid(UAI_GLS_nest_low))

# save model
saveRDS(UAI_GLS_nest_low, here("Models/UAI", "UAI_GLS_nest_low.rds"))

#########################################
# Filter out all species that use both HIGH and LOW nest sites 
# we are doing this to test whether the "non-ambiguous" nesting species the ones driving the significant relationship in the previous section


# apply filter 
UAI_NestLow_only <- C_Nest_dat2 %>% filter(!is.na(aveUAI)) %>% 
  filter(!is.na(NestSite_Low)) %>% 
  filter(!(NestSite_Low == 1 & NestSite_High == 1)) %>% as.data.frame()
length(UAI_NestLow_only$NestSite_Low)
# 575 species with UAI and ONLY NestSite_Low (not also high nesters)
796-575 # = 221 --> there are 221 species that were both High and Low nesters and had UAI scores 

###### add and pair tree

# add rownames to data
row.names(UAI_NestLow_only) <- UAI_NestLow_only$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_NestLow_only_phydat <- treedata(tree_out, UAI_NestLow_only, sort=T)

UAI_NestLow_only_phy <- UAI_NestLow_only_phydat$phy
UAI_NestLow_only_dat <- as.data.frame(UAI_NestLow_only_phydat$data)

str(UAI_NestLow_only_dat)
length(UAI_NestLow_only_dat$NestSite_Low)
#575

# convert traits of interest to numeric
UAI_NestLow_only_dat$aveUAI <- as.numeric(UAI_NestLow_only_dat$aveUAI)
UAI_NestLow_only_dat$Mass_log <- as.numeric(UAI_NestLow_only_dat$Mass_log)
UAI_NestLow_only_dat$NestSite_Low <- as.numeric(UAI_NestLow_only_dat$NestSite_Low)

# run a phylogenetic linear model
UAI_GLS_nest_low_only <- gls(aveUAI~ NestSite_Low + Mass_log, data = UAI_NestLow_only_dat, 
                          correlation = corPagel(0.5, phy = UAI_NestLow_only_phy, fixed = F, form = ~Species_Jetz), 
                          method = "ML") 


# model summary and results
summary(UAI_GLS_nest_low_only) 
confint(UAI_GLS_nest_low_only, level = 0.95)
confint(UAI_GLS_nest_low_only, level = 0.85)

# model diagnostics
check_model(UAI_GLS_nest_low_only) 
qqnorm(resid(UAI_GLS_nest_low_only)) 
qqline(resid(UAI_GLS_nest_low_only))
hist(resid(UAI_GLS_nest_low_only))

# save model
saveRDS(UAI_GLS_nest_low_only, here("Models/UAI", "UAI_GLS_nest_low_only.rds"))


######################## MUTI and Nest Site LOW ##########################
# 0 = not low
# 1 = low

# create a new data frame that contains only species with both MUTI and nest site low
MUTI_NestLow <- C_Nest_dat2 %>% filter(!is.na(MUTIscore)) %>%
  filter(!is.na(NestSite_Low)) %>% as.data.frame()
length(MUTI_NestLow$NestSite_Low)
#130 species with MUTIscore and NestSite_Low

###### add and pair tree

# add rownames to data
row.names(MUTI_NestLow) <- MUTI_NestLow$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_NestLow_phydat <- treedata(tree_out, MUTI_NestLow, sort=T)

MUTI_NestLow_phy <- MUTI_NestLow_phydat$phy
MUTI_NestLow_dat <- as.data.frame(MUTI_NestLow_phydat$data)

str(MUTI_NestLow_dat)
length(MUTI_NestLow_dat$NestSite_Low)
#130

# convert traits of interest to numeric
MUTI_NestLow_dat$MUTIscore <- as.numeric(MUTI_NestLow_dat$MUTIscore)
MUTI_NestLow_dat$Mass_log <- as.numeric(MUTI_NestLow_dat$Mass_log)
MUTI_NestLow_dat$NestSite_Low <- as.numeric(MUTI_NestLow_dat$NestSite_Low)

# run a phylogenetic linear model
MUTI_GLS_nest_low <- gls(MUTIscore~ NestSite_Low + Mass_log, data = MUTI_NestLow_dat, 
                        correlation = corPagel(0.5, phy = MUTI_NestLow_phy, fixed = F, form = ~Species_Jetz), 
                        method = "ML") 

# model summary and results
summary(MUTI_GLS_nest_low) 
confint(MUTI_GLS_nest_low, level = 0.95)
confint(MUTI_GLS_nest_low, level = 0.85)

# model diagnostics
check_model(MUTI_GLS_nest_low) 
qqnorm(resid(MUTI_GLS_nest_low)) 
qqline(resid(MUTI_GLS_nest_low))
hist(resid(MUTI_GLS_nest_low))

# save model for easy retrieval 
saveRDS(MUTI_GLS_nest_low, here("Models/MUTI", "MUTI_GLS_nest_low.rds"))


#########################################
# Filter out all species that use both HIGH and LOW nest sites 
# we are doing this to test whether the "non-ambiguous" nesting species the ones driving the significant relationship in the previous section

# apply filter 
MUTI_NestLow_only <- C_Nest_dat2 %>% filter(!is.na(MUTIscore)) %>% 
  filter(!is.na(NestSite_Low)) %>% 
  filter(!(NestSite_Low == 1 & NestSite_High == 1)) %>% as.data.frame()
length(MUTI_NestLow_only$NestSite_Low)
# 91 species 
130-91 # = 39 --> there are 39 species that were both High and Low nesters and had MUTI scores 

###### add and pair tree

# add rownames to data
row.names(MUTI_NestLow_only) <- MUTI_NestLow_only$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_NestLow_only_phydat <- treedata(tree_out, MUTI_NestLow_only, sort=T)

MUTI_NestLow_only_phy <- MUTI_NestLow_only_phydat$phy
MUTI_NestLow_only_dat <- as.data.frame(MUTI_NestLow_only_phydat$data)

str(MUTI_NestLow_only_dat)
length(MUTI_NestLow_only_dat$NestSite_Low)
#91


# convert traits of interest to numeric
MUTI_NestLow_only_dat$MUTIscore <- as.numeric(MUTI_NestLow_only_dat$MUTIscore)
MUTI_NestLow_only_dat$Mass_log <- as.numeric(MUTI_NestLow_only_dat$Mass_log)
MUTI_NestLow_only_dat$NestSite_Low <- as.numeric(MUTI_NestLow_only_dat$NestSite_Low)

# Model is not working with corPagel starting point = 0.5 and fixed = F
# we need to find a value of lambda value to fix in the model 

# Create an empty vector to store AIC values
AIC_values <- numeric()

# Loop through different values of the parameter for corPagel
for (i in seq(0, 1, by = 0.1)) {
  # Fit the gls model with the current value of i
  model <- gls(MUTIscore ~ NestSite_Low + Mass_log, 
               data = MUTI_NestLow_only_dat, 
               correlation = corPagel(i, phy = MUTI_NestLow_only_phy, fixed = TRUE, form = ~Species_Jetz), 
               method = "ML")
  # Extract AIC value and store it in the vector
  AIC_values <- c(AIC_values, AIC(model))
}

# Print AIC values
print(AIC_values)
#0.3 = best AIC score 


# run a phylogenetic linear model
MUTI_GLS_nest_low_only <- gls(MUTIscore~ NestSite_Low + Mass_log, data = MUTI_NestLow_only_dat, 
                         correlation = corPagel(0.3, phy = MUTI_NestLow_only_phy, fixed = T, form = ~Species_Jetz), 
                         method = "ML") 

# model summary and results
summary(MUTI_GLS_nest_low_only) 
confint(MUTI_GLS_nest_low_only, level = 0.95)

# model diagnostics
check_model(MUTI_GLS_nest_low_only) 
qqnorm(resid(MUTI_GLS_nest_low_only)) 
qqline(resid(MUTI_GLS_nest_low_only))
hist(resid(MUTI_GLS_nest_low_only))

# save model for easy retrieval 
saveRDS(MUTI_GLS_nest_low_only, here("Models/MUTI", "MUTI_GLS_nest_low_only.rds"))


######################## UN and Nest Site LOW ##########################
# 0 = not low
# 1 = low

# create a new data frame that contains only species with both UN and nest site low
UN_NestLow <- C_Nest_dat2 %>% filter(!is.na(Urban)) %>%
  filter(!is.na(NestSite_Low)) %>% column_to_rownames(., var = "Species_Jetz")
length(UN_NestLow$NestSite_Low)
# 129 species with UN and NestSite_Low

###### add and pair tree

tree_out <- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UN_NestLow_phydat <- treedata(tree_out, UN_NestLow, sort=T)

UN_NestLow_phy <- UN_NestLow_phydat$phy
UN_NestLow_dat <- as.data.frame(UN_NestLow_phydat$data)

str(UN_NestLow_dat)
length(UN_NestLow_dat$NestSite_Low)
#129


# convert traits of interest to numeric
UN_NestLow_dat$Urban <- as.numeric(UN_NestLow_dat$Urban)
UN_NestLow_dat$Mass_log <- as.numeric(UN_NestLow_dat$Mass_log)
UN_NestLow_dat$NestSite_Low <- as.numeric(UN_NestLow_dat$NestSite_Low)

# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
phyglm_UN_nest_low_scale <- phyloglm( Urban ~ NestSite_Low + scale(Mass_log), 
                                      data = UN_NestLow_dat, 
                                      phy = UN_NestLow_phy, 
                                      boot = 1000) 

summary(phyglm_UN_nest_low_scale)
# fails to converge
# alpha at upper bounds


# print AIC values for models with different upper bounds
# intervals of 0.1 from 0 up to 4
for (i in seq(0, 4, by = 0.1)) {
  print(phyloglm(Urban ~  NestSite_Low + scale(Mass_log), 
                 data = UN_NestLow_dat, 
                 phy = UN_NestLow_phy,
                 log.alpha.bound = i)$aic)
}
# AIC values support models with larger values of alpha (low phylo signal)


# try fixing alpha at upper bounds
# giving the model a little more searching space for alpha because AIC is actually lowest slightly below log.alpha.bounds = 4
exp(3.7)/(phyglm_UN_nest_low_scale$mean.tip.height) # equals 0.41. Using this as start.alpha

set.seed(701)
phyglm_UN_nest_low_fix <- phyloglm( Urban ~ NestSite_Low + scale(Mass_log), 
                                    data = UN_NestLow_dat, 
                                    phy = UN_NestLow_phy,
                                    log.alpha.bound = 4,
                                    start.alpha = 0.41,
                                    boot = 1000)
summary(phyglm_UN_nest_low_fix)
# model converges
# alpha at upper bounds (= 0.559)


# save model
saveRDS(phyglm_UN_nest_low_fix, here("Models/UN", "phyglm_UN_nest_low_fix.rds"))
# load model
phyglm_UN_nest_low_fix <- readRDS(here("Models/UN", "phyglm_UN_nest_low_fix.rds"))


# compare with non-phylogenetic model
glm_UN_nest_low <- logistf(Urban ~ NestSite_Low + scale(Mass_log), 
                           data = UN_NestLow)
summary(glm_UN_nest_low)
# some differences for coefficients but same conclusions reached


# get alpha, t, and half life for the model
(phyglm_UN_nest_low_fix$mean.tip.height) # t
(alpha_Nlow <- phyglm_UN_nest_low_fix$alpha) # alpha
(hl_NLow <- log(2)/alpha_Nlow) # half life
# compared to t, this is a small half life

#########################################
# Filter out all species that use both HIGH and LOW nest sites 
# we are doing this to test whether the "non-ambiguous" nesting species the ones driving the significant relationship in the previous section

UN_NestLow_only <- C_Nest_dat2 %>% filter(!is.na(Urban)) %>% 
  filter(!is.na(NestSite_Low)) %>% 
  filter(!(NestSite_Low == 1 & NestSite_High == 1)) %>%
  column_to_rownames(., var = "Species_Jetz")
length(UN_NestLow_only$NestSite_Low)
#104 species 
129-104 # = 25 --> there are 25 species that were both High and Low nesters and had UN scores 

###### add and pair tree

tree_out <- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UN_NestLow_only_phydat <- treedata(tree_out, UN_NestLow_only, sort=T)

UN_NestLow_only_phy <- UN_NestLow_only_phydat$phy
UN_NestLow_only_dat <- as.data.frame(UN_NestLow_only_phydat$data)

str(UN_NestLow_only_dat)
length(UN_NestLow_only_dat$NestSite_Low)
#104


# convert traits of interest to numeric
UN_NestLow_only_dat$Urban <- as.numeric(UN_NestLow_only_dat$Urban)
UN_NestLow_only_dat$Mass_log <- as.numeric(UN_NestLow_only_dat$Mass_log)
UN_NestLow_only_dat$NestSite_Low <- as.numeric(UN_NestLow_only_dat$NestSite_Low)


# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
set.seed(438)
phyglm_UN_nest_low_only_scale <- phyloglm( Urban ~ NestSite_Low + scale(Mass_log), 
                                           data = UN_NestLow_only_dat, 
                                           phy = UN_NestLow_only_phy,
                                           boot = 1000)


summary(phyglm_UN_nest_low_only_scale)
# alpha at upper bounds
# the p-value for NestSite_low flagged as significant
# but the bootstrapped CI interval strongly overlaps zero
# we will use the bootstrapped 95% CI here over the p-value


# save model
saveRDS(phyglm_UN_nest_low_only_scale, here("Models/UN", "phyglm_UN_nest_low_only_scale.rds"))
# load model
phyglm_UN_nest_low_only_scale <- readRDS(here("Models/UN", "phyglm_UN_nest_low_only_scale.rds"))


# get alpha, t, and half life for the model
(phyglm_UN_nest_low_only_scale$mean.tip.height) # t
(alpha_Nlowonly <- phyglm_UN_nest_low_only_scale$alpha) # alpha
(hl_Nlowonly <- log(2)/alpha_Nlowonly) # half life
# compared to t, this is a small half-life -> low phylogenetic signal


##########################################################################
##########################################################################

######################## UAI and Nest Site HIGH ##########################
# 0 = not high
# 1 = HIGH

# create a new data frame that contains only species with both UAI and nest site high
UAI_NestHigh <- C_Nest_dat2 %>% filter(!is.na(aveUAI)) %>% 
  filter(!is.na(NestSite_High)) %>% as.data.frame()
length(UAI_NestHigh$NestSite_High)
#796 species with UAI and NestSite_High

###### add and pair tree

# add rownames to data
row.names(UAI_NestHigh) <- UAI_NestHigh$Species_Jetz

tree_out <- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_NestHigh_phydat <- treedata(tree_out, UAI_NestHigh, sort=T)

UAI_NestHigh_phy <- UAI_NestHigh_phydat$phy
UAI_NestHigh_dat <- as.data.frame(UAI_NestHigh_phydat$data)

str(UAI_NestHigh_dat)
length(UAI_NestHigh_dat$NestSite_High)
#796


# convert traits of interest to numeric
UAI_NestHigh_dat$aveUAI <- as.numeric(UAI_NestHigh_dat$aveUAI)
UAI_NestHigh_dat$Mass_log <- as.numeric(UAI_NestHigh_dat$Mass_log)
UAI_NestHigh_dat$NestSite_High <- as.numeric(UAI_NestHigh_dat$NestSite_High)

# run phylogenetic linear model
UAI_GLS_nest_high <- gls(aveUAI~ NestSite_High + Mass_log, data = UAI_NestHigh_dat, 
                         correlation = corPagel(0.5, phy = UAI_NestHigh_phy, fixed = F, form = ~Species_Jetz), 
                         method = "ML") 

# model summary and results
summary(UAI_GLS_nest_high) 
confint(UAI_GLS_nest_high)

# model diagnostics
check_model(UAI_GLS_nest_high) 
qqnorm(resid(UAI_GLS_nest_high)) 
qqline(resid(UAI_GLS_nest_high))
hist(resid(UAI_GLS_nest_high))

# save model
saveRDS(UAI_GLS_nest_high, here("Models/UAI", "UAI_GLS_nest_high.rds"))


#########################################
# Filter out all species that use both HIGH and LOW nest sites 
# we are doing this to test whether the "non-ambiguous" nesting species the ones driving the significant relationship in the previous section

# apply filter 
UAI_NestHigh_only <- C_Nest_dat2 %>% filter(!is.na(aveUAI)) %>%
  filter(!is.na(NestSite_High)) %>% 
  filter(!(NestSite_Low == 1 & NestSite_High == 1)) %>% as.data.frame()
length(UAI_NestHigh_only$NestSite_High)
#575 species 
796-575 # = 221 --> there are 221 species that were both High and Low nesters and had UAI scores 

###### add and pair tree

# add rownames to data
row.names(UAI_NestHigh_only) <- UAI_NestHigh_only$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_NestHigh_only_phydat <- treedata(tree_out, UAI_NestHigh_only, sort=T)

UAI_NestHigh_only_phy <- UAI_NestHigh_only_phydat$phy
UAI_NestHigh_only_dat <- as.data.frame(UAI_NestHigh_only_phydat$data)

str(UAI_NestHigh_only_dat)
length(UAI_NestHigh_only_dat$NestSite_High)
#575

### convert traits of interest to numeric
UAI_NestHigh_only_dat$aveUAI <- as.numeric(UAI_NestHigh_only_dat$aveUAI)
UAI_NestHigh_only_dat$Mass_log <- as.numeric(UAI_NestHigh_only_dat$Mass_log)
UAI_NestHigh_only_dat$NestSite_High <- as.numeric(UAI_NestHigh_only_dat$NestSite_High)

# run phylogenetic linear model
UAI_GLS_nest_high_only <- gls(aveUAI~ NestSite_High + Mass_log, data = UAI_NestHigh_only_dat, 
                              correlation = corPagel(0.5, phy = UAI_NestHigh_only_phy, fixed = F, form = ~Species_Jetz), 
                              method = "ML") 

# model summary and results
summary(UAI_GLS_nest_high_only) 
confint(UAI_GLS_nest_high_only)

# model diagnostics
check_model(UAI_GLS_nest_high_only) 
qqnorm(resid(UAI_GLS_nest_high_only)) 
qqline(resid(UAI_GLS_nest_high_only))
hist(resid(UAI_GLS_nest_high_only))

# save model
saveRDS(UAI_GLS_nest_high_only, here("Models/UAI", "UAI_GLS_nest_high_only.rds"))


######################## MUTI and Nest Site HIGH ##########################
# 0 = not high
# 1 = HIGH

# create a new data frame that contains only species with both MUTI and nest site high
MUTI_NestHigh <- C_Nest_dat2 %>% filter(!is.na(MUTIscore)) %>% 
  filter(!is.na(NestSite_High)) %>% as.data.frame()
length(MUTI_NestHigh$NestSite_High)
# 130 species with UAI and NestSite_High

###### add and pair tree

# add rownames to data
row.names(MUTI_NestHigh) <- MUTI_NestHigh$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_NestHigh_phydat <- treedata(tree_out, MUTI_NestHigh, sort=T)

MUTI_NestHigh_phy <- MUTI_NestHigh_phydat$phy
MUTI_NestHigh_dat <- as.data.frame(MUTI_NestHigh_phydat$data)

str(MUTI_NestHigh_dat)
length(MUTI_NestHigh_dat$NestSite_High)
# 130


# convert traits of interest to numeric
MUTI_NestHigh_dat$MUTIscore <- as.numeric(MUTI_NestHigh_dat$MUTIscore)
MUTI_NestHigh_dat$Mass_log <- as.numeric(MUTI_NestHigh_dat$Mass_log)
MUTI_NestHigh_dat$NestSite_High <- as.numeric(MUTI_NestHigh_dat$NestSite_High)


# run phylogenetic linear model
MUTI_GLS_nest_high <- gls(MUTIscore~ NestSite_High + Mass_log, data = MUTI_NestHigh_dat, 
                          correlation = corPagel(0.5, phy = MUTI_NestHigh_phy, fixed = F, form = ~Species_Jetz), 
                          method = "ML") 

# model summary and results
summary(MUTI_GLS_nest_high) 
confint(MUTI_GLS_nest_high)
confint(MUTI_GLS_nest_high, level = 0.85)

# model diagnostics
check_model(MUTI_GLS_nest_high) 
qqnorm(resid(MUTI_GLS_nest_high)) 
qqline(resid(MUTI_GLS_nest_high))
hist(resid(MUTI_GLS_nest_high))

# save model
saveRDS(MUTI_GLS_nest_high, here("Models/MUTI", "MUTI_GLS_nest_high.rds"))

#########################################
# Filter out all species that use both HIGH and LOW nest sites 
# we are doing this to test whether the "non-ambiguous" nesting species the ones driving the significant relationship in the previous section

# apply filter 
MUTI_NestHigh_only <- C_Nest_dat2 %>% filter(!is.na(MUTIscore)) %>%
  filter(!is.na(NestSite_High)) %>% 
  filter(!(NestSite_Low == 1 & NestSite_High == 1)) %>% as.data.frame()
length(MUTI_NestHigh_only$NestSite_High)
# 91 species
130-91 # = 39 --> there are 39 species that were both High and Low nesters and had MUTI scores 

###### add and pair tree

# add rownames to data
row.names(MUTI_NestHigh_only) <- MUTI_NestHigh_only$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_NestHigh_only_phydat <- treedata(tree_out, MUTI_NestHigh_only, sort=T)

MUTI_NestHigh_only_phy <- MUTI_NestHigh_only_phydat$phy
MUTI_NestHigh_only_dat <- as.data.frame(MUTI_NestHigh_only_phydat$data)

str(MUTI_NestHigh_only_dat)
length(MUTI_NestHigh_only_dat$NestSite_High)
# 91


# convert traits of interest to numeric
MUTI_NestHigh_only_dat$MUTIscore <- as.numeric(MUTI_NestHigh_only_dat$MUTIscore)
MUTI_NestHigh_only_dat$Mass_log <- as.numeric(MUTI_NestHigh_only_dat$Mass_log)
MUTI_NestHigh_only_dat$NestSite_High <- as.numeric(MUTI_NestHigh_only_dat$NestSite_High)

# model is not working with corPagel starting point = 0.5 and fixed = F
# we need to find a lambda value to fix in the model 

# Create an empty vector to store AIC values
AIC_values <- numeric()

# Loop through different values of the parameter for corPagel
for (i in seq(0, 1, by = 0.1)) {
  # Fit the gls model with the current value of i
  model <- gls(MUTIscore ~ NestSite_High + Mass_log, 
               data = MUTI_NestHigh_only_dat, 
               correlation = corPagel(i, phy = MUTI_NestHigh_only_phy, fixed = TRUE, form = ~Species_Jetz), 
               method = "ML")
  # Extract AIC value and store it in the vector
  AIC_values <- c(AIC_values, AIC(model))
}

# Print AIC values
print(AIC_values)
# 0.3 = best AIC score 


# run phylogenetic linear model
MUTI_GLS_nest_high_only <- gls(MUTIscore~ NestSite_High + Mass_log, data = MUTI_NestHigh_only_dat, 
                              correlation = corPagel(0.3, phy= MUTI_NestHigh_only_phy, fixed = T, form = ~Species_Jetz), 
                              method = "ML") 

# model summary and results
summary(MUTI_GLS_nest_high_only) 
confint(MUTI_GLS_nest_high_only)

# model diagnostics
check_model(MUTI_GLS_nest_high_only) 
qqnorm(resid(MUTI_GLS_nest_high_only)) 
qqline(resid(MUTI_GLS_nest_high_only))
hist(resid(MUTI_GLS_nest_high_only))

# save model
saveRDS(MUTI_GLS_nest_high_only, here("Models/MUTI", "MUTI_GLS_nest_high_only.rds"))


######################## UN and Nest Site HIGH ##########################
# 0 = not high
# 1 = HIGH

# create a new data frame that contains only species with both UN and nest site high
UN_NestHigh <- C_Nest_dat2 %>% filter(!is.na(Urban)) %>% 
  filter(!is.na(NestSite_High)) %>% column_to_rownames(., var ="Species_Jetz")
length(UN_NestHigh$NestSite_High)
#129 species with UN and NestSite_High


###### add and pair tree

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UN_NestHigh_phydat <- treedata(tree_out, UN_NestHigh, sort=T)

UN_NestHigh_phy <- UN_NestHigh_phydat$phy
UN_NestHigh_dat <- as.data.frame(UN_NestHigh_phydat$data)

str(UN_NestHigh_dat)
length(UN_NestHigh_dat$NestSite_High)
# 129


### convert traits of interest to numeric
UN_NestHigh_dat$Urban <- as.numeric(UN_NestHigh_dat$Urban)
UN_NestHigh_dat$Mass_log <- as.numeric(UN_NestHigh_dat$Mass_log)
UN_NestHigh_dat$NestSite_High <- as.numeric(UN_NestHigh_dat$NestSite_High)


# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
phyglm_UN_nest_high_scale <- phyloglm( Urban ~ NestSite_High + scale(Mass_log), 
                                       data = UN_NestHigh_dat, 
                                       phy = UN_NestHigh_phy, 
                                       boot = 1000)
summary(phyglm_UN_nest_high_scale) # this fails to converge


# print AIC values for models with different upper bounds
# intervals of 0.1 from 0 up to 4
for (i in seq(0, 4, by = 0.1)) {
  print(phyloglm(Urban ~  NestSite_High + scale(Mass_log), 
                 data = UN_NestHigh_dat, 
                 phy = UN_NestHigh_phy, 
                 log.alpha.bound = i)$aic)
}
# AIC values support models with larger values of alpha (low phylo signal)


# try fixing alpha at upper bounds
# give the model a little more searching space for alpha because AIC is actually lowest slightly below log.alpha.bounds = 4
exp(3.8)/(phyglm_UN_nest_high_scale$mean.tip.height) # equals 0.45. Using this as start.alpha

set.seed(782)
phyglm_UN_nest_high_fix <- phyloglm( Urban ~ NestSite_High + scale(Mass_log), 
                                     data = UN_NestHigh_dat, 
                                     phy = UN_NestHigh_phy,  
                                     log.alpha.bound = 4,
                                     start.alpha = 0.45,
                                     boot = 1000)
summary(phyglm_UN_nest_high_fix)
# model converges


# save model
saveRDS(phyglm_UN_nest_high_fix, here("Models/UN", "phyglm_UN_nest_high_fix.rds"))
# load model
phyglm_UN_nest_high_fix <- readRDS(here("Models/UN", "phyglm_UN_nest_high_fix.rds"))


# look at a non-phylogenetic logistic model
glm_UN_nest_high <- logistf(Urban ~ NestSite_High + scale(Mass_log), 
                            data = UN_NestHigh)

summary(glm_UN_nest_high) # we reach same conclusions


# get alpha, t, and half life for the model
(phyglm_UN_nest_high_fix$mean.tip.height) # t
(alpha_Nhigh <- phyglm_UN_nest_high_fix$alpha) # alpha
(hl_Nhigh <- log(2)/alpha_Nhigh) # half life
#compared to t, this is a small half life


#########################################
# Filter out all species that use both HIGH and LOW nest sites 
# we are doing this to test whether the "non-ambiguous" nesting species the ones driving the significant relationship in the previous section

# apply filter 
UN_NestHigh_only <- C_Nest_dat2 %>% filter(!is.na(Urban))  %>% 
  filter(!is.na(NestSite_High)) %>% 
  filter(!(NestSite_Low == 1 & NestSite_High == 1)) %>%
  column_to_rownames(., var="Species_Jetz")
length(UN_NestHigh_only$NestSite_High)
# 104 species
129-104 # = 25 --> there are 25 species that were both High and Low nesters and had UN scores 

###### add and pair tree

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UN_NestHigh_only_phydat <- treedata(tree_out, UN_NestHigh_only, sort=T)

UN_NestHigh_only_phy <- UN_NestHigh_only_phydat$phy
UN_NestHigh_only_dat <- as.data.frame(UN_NestHigh_only_phydat$data)

str(UN_NestHigh_only_dat)
length(UN_NestHigh_only_dat$NestSite_High)
# 104

# convert traits of interest to numeric
UN_NestHigh_only_dat$Urban <- as.numeric(UN_NestHigh_only_dat$Urban)
UN_NestHigh_only_dat$Mass_log <- as.numeric(UN_NestHigh_only_dat$Mass_log)
UN_NestHigh_only_dat$NestSite_High <- as.numeric(UN_NestHigh_only_dat$NestSite_High)


# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
set.seed(107)
phyglm_UN_nest_high_only_scale <- phyloglm( Urban ~ NestSite_High + scale(Mass_log), 
                                            data = UN_NestHigh_only_dat, 
                                            phy = UN_NestHigh_only_phy, 
                                            boot = 1000) 
# model converges
# alpha at upper bounds
summary(phyglm_UN_nest_high_only_scale)
# nest site high is no longer important


# save model
saveRDS(phyglm_UN_nest_high_only_scale, here("Models/UN", "phyglm_UN_nest_high_only_scale.rds"))
# load model
phyglm_UN_nest_high_only_scale <- readRDS(here("Models/UN", "phyglm_UN_nest_high_only_scale.rds"))


# get alpha, t, and half life for the model
(phyglm_UN_nest_high_only_scale$mean.tip.height) # t
(alpha_Nhighonly <- phyglm_UN_nest_high_only_scale$alpha) # alpha
(hl_Nhighonly <- log(2)/alpha_Nhighonly) # half life
# compared to t, this is a small half life -> low phylogenetic signal



##########################################################################
##########################################################################

######################## UAI and Nest Safety ##########################

# create a new data frame that contains only species with both UAI and nest safety
UAI_NestSafety <- C_Nest_dat2 %>% filter(!is.na(aveUAI)) %>% 
  filter(!is.na(nest.safety)) %>% as.data.frame()
length(UAI_NestSafety$nest.safety)
#766 species with UAI and nest.safety

###### add and pair tree

# add rownames to data
row.names(UAI_NestSafety) <- UAI_NestSafety$Species_Jetz

tree_out <- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_NestSafety_phydat <- treedata(tree_out, UAI_NestSafety, sort=T)

UAI_NestSafety_phy <- UAI_NestSafety_phydat$phy
UAI_NestSafety_dat <- as.data.frame(UAI_NestSafety_phydat$data)

str(UAI_NestSafety_dat)
length(UAI_NestSafety_dat$nest.safety)
#766


### convert traits of interest to numeric
UAI_NestSafety_dat$aveUAI <- as.numeric(UAI_NestSafety_dat$aveUAI)
UAI_NestSafety_dat$Mass_log <- as.numeric(UAI_NestSafety_dat$Mass_log)
UAI_NestSafety_dat$nest.safety <- as.numeric(UAI_NestSafety_dat$nest.safety)


# run phylogenetic linear model
UAI_GLS_nest_safety <- gls(aveUAI~ nest.safety + Mass_log, data = UAI_NestSafety_dat, 
                            correlation = corPagel(0.5, phy = UAI_NestSafety_phy, fixed=F, form = ~Species_Jetz), 
                            method = "ML") 

# model summary and results
summary(UAI_GLS_nest_safety) 
confint(UAI_GLS_nest_safety)

# model diagnostics
check_model(UAI_GLS_nest_safety) 
qqnorm(resid(UAI_GLS_nest_safety)) 
qqline(resid(UAI_GLS_nest_safety))
hist(resid(UAI_GLS_nest_safety))

# save model
saveRDS(UAI_GLS_nest_safety, here("Models/UAI", "UAI_GLS_nest_safety.rds"))


######################## MUTI and Nest Safety ##########################

# create a new data frame that contains only species with both MUTI and nest safety
MUTI_NestSafety <- C_Nest_dat2 %>% filter(!is.na(MUTIscore))  %>%
  filter(!is.na(nest.safety)) %>% as.data.frame()
length(MUTI_NestSafety$nest.safety)
# 127 species with MUTIscore and nest.safety

###### add and pair tree

# add rownames to data
row.names(MUTI_NestSafety) <- MUTI_NestSafety$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_NestSafety_phydat <- treedata(tree_out, MUTI_NestSafety, sort=T)

MUTI_NestSafety_phy <- MUTI_NestSafety_phydat$phy
MUTI_NestSafety_dat <- as.data.frame(MUTI_NestSafety_phydat$data)

str(MUTI_NestSafety_dat)
length(MUTI_NestSafety_dat$nest.safety)
#127

# convert traits of interest to numeric
MUTI_NestSafety_dat$MUTIscore <- as.numeric(MUTI_NestSafety_dat$MUTIscore)
MUTI_NestSafety_dat$Mass_log <- as.numeric(MUTI_NestSafety_dat$Mass_log)
MUTI_NestSafety_dat$nest.safety <- as.numeric(MUTI_NestSafety_dat$nest.safety)


# run a phylogenetic linear model
MUTI_GLS_nest_safety <- gls(MUTIscore~ nest.safety + Mass_log, data = MUTI_NestSafety_dat, 
                           correlation = corPagel(0.5, phy = MUTI_NestSafety_phy, fixed = F, form = ~Species_Jetz), 
                           method = "ML") 

# model summary and results
summary(MUTI_GLS_nest_safety) 
confint(MUTI_GLS_nest_safety)

# model diagnostics
check_model(MUTI_GLS_nest_safety) 
qqnorm(resid(MUTI_GLS_nest_safety)) 
qqline(resid(MUTI_GLS_nest_safety))
hist(resid(MUTI_GLS_nest_safety))

# save model
saveRDS(MUTI_GLS_nest_safety, here("Models/MUTI", "MUTI_GLS_nest_safety.rds"))


######################## UN and Nest Safety ##########################


# create a new data frame that contains only species with both UN and nest safety
UN_NestSafety <- C_Nest_dat2 %>% filter(!is.na(Urban))  %>% 
  filter(!is.na(nest.safety)) %>% column_to_rownames(., var="Species_Jetz")
length(UN_NestSafety$nest.safety)
#129 species with Urban and nest.safety

###### add and pair tree
tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UN_NestSafety_phydat <- treedata(tree_out, UN_NestSafety, sort=T)

UN_NestSafety_phy <- UN_NestSafety_phydat$phy
UN_NestSafety_dat <- as.data.frame(UN_NestSafety_phydat$data)

str(UN_NestSafety_dat)
length(UN_NestSafety_dat$nest.safety)
#129

### convert traits of interest to numeric
UN_NestSafety_dat$Urban <- as.numeric(UN_NestSafety_dat$Urban)
UN_NestSafety_dat$Mass_log <- as.numeric(UN_NestSafety_dat$Mass_log)
UN_NestSafety_dat$nest.safety <- as.numeric(UN_NestSafety_dat$nest.safety)

# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
set.seed(753)
phyglm_UN_nest_safety_scale <- phyloglm( Urban ~ scale(nest.safety) + scale(Mass_log), 
                                         data = UN_NestSafety_dat, 
                                         phy = UN_NestSafety_phy, 
                                         boot = 1000) 
summary(phyglm_UN_nest_safety_scale)
# this model converges
# this is a positive trend for nest.safety based on bootstrapped 95% CI

# save model
saveRDS(phyglm_UN_nest_safety_scale, here("Models/UN", "phyglm_UN_nest_safety_scale.rds"))
# load model
phyglm_UN_nest_safety_scale <- readRDS(here("Models/UN", "phyglm_UN_nest_safety_scale.rds"))


# get alpha, t, and half life for the model
(phyglm_UN_nest_safety_scale$mean.tip.height) # t
(alpha_Nsafe <- phyglm_UN_nest_safety_scale$alpha) # alpha
(hl_Nsafe <- log(2)/alpha_Nsafe) # half life
# compared to t, this is small half life
