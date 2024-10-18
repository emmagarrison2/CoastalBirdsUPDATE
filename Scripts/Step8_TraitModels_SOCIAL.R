##########################################################
######### Coastal Birds Urban Tolerance Project ##########
##########################################################
# Step 8 Phylogenetic Trait Models - SOCIAL TRAITS
# Authors: Sarah L. Jennings, Emma M. Garrison
##########################################################
# The objective of this script is to run phylogenetic models using social traits
# 2 social traits: territoriality and cooperative breeding


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

# load  "Coastal_Species_Social.rds"
C_Social_dat <- readRDS(here("Outputs", "Coastal_Species_Social.rds"))
str(C_Social_dat)

C_Social_dat2 <- C_Social_dat %>%
  mutate(Species_Jetz  = str_replace(Species_Jetz, " ", "_"))
str(C_Social_dat2)

C_Social_dat2$Urban <- ifelse(C_Social_dat2$Urban == "U", 1, 0)
View(C_Social_dat2)
colnames(C_Social_dat2)


######################## UAI and Territoriality ##########################

# create a new data frame that contains only species with both UAI and territoriality values 
UAI_Terr <- C_Social_dat2 %>% filter(!is.na(aveUAI)) %>% filter(!is.na(territoriality))
length(UAI_Terr$territoriality)
# 766 species with UAI and territoriality


###### add and pair tree

# add rownames to data
row.names(UAI_Terr) <- UAI_Terr$Species_Jetz

tree_out <- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_Terr_phydat <- treedata(tree_out, UAI_Terr, sort=T)

UAI_Terr_phy <- UAI_Terr_phydat$phy
UAI_Terr_dat <- as.data.frame(UAI_Terr_phydat$data)

str(UAI_Terr_dat)
length(UAI_Terr_dat$territoriality)
# 766


### convert traits of interest to numeric
UAI_Terr_dat$aveUAI <- as.numeric(UAI_Terr_dat$aveUAI)
UAI_Terr_dat$Mass_log <- as.numeric(UAI_Terr_dat$Mass_log)
UAI_Terr_dat$territoriality <- as.numeric(UAI_Terr_dat$territoriality)


# Run phylogenetic linear model
UAI_GLS_territory <- gls(aveUAI ~ territoriality + Mass_log, data = UAI_Terr_dat, 
                      correlation = corPagel(0.5, phy=UAI_Terr_phy,fixed=F, form = ~Species_Jetz), 
                      method = "ML") 

# model summary and results
summary(UAI_GLS_territory) 
confint(UAI_GLS_territory)

# model diagnostics
check_model(UAI_GLS_territory) 
qqnorm(resid(UAI_GLS_territory)) 
qqline(resid(UAI_GLS_territory)) 
hist(resid(UAI_GLS_territory)) 

# save model for easy retrieval 
saveRDS(UAI_GLS_territory, here("Models/UAI", "UAI_GLS_territory.rds"))

### make a factor version of territoriality 
UAI_Terr_dat$territoriality_f <- as.factor(UAI_Terr_dat$territoriality)


# Run phylogenetic linear model
UAI_GLS_territory_f <- gls(aveUAI~ territoriality_f + Mass_log, data = UAI_Terr_dat, 
                         correlation = corPagel(0.5, phy = UAI_Terr_phy, fixed=F, form = ~Species_Jetz), 
                         method = "ML") 

summary(UAI_GLS_territory_f)

# compare levels of territoriality
pred.terr <- predict_response(UAI_GLS_territory_f, "territoriality_f", margin = "marginalmeans") # get CIs for different levels of territoriality
pred.terr # this gives the predicted mean for each level of territoriality and the 95% CI

# run comparison of each pair of levels
test_predictions(pred.terr) # 0 and 1 are the levels that are significantly different
# territoriality = 0 has higher UAI scores than territoriality = 1
# there is a trend that 0 is higher than 2, but it's not significant using alpha = 0.05

# make a quick plot to visualize the relationship
territoriality.plot <- plot(pred.terr, show_data=T, jitter=0.5, dot_size = 1.5, line_size= 1.5, colors =  "#21918c", alpha = 0.4) + 
  theme_classic(base_size = 12, base_family = "") +
  labs(x = "Territoriality" , y = "UAI Score")
territoriality.plot

 ######################## MUTI and Territoriality ##########################

# create a new data frame that contains only species with both MUTI and territoriality values
MUTI_Terr <- C_Social_dat2 %>% filter(!is.na(MUTIscore)) %>% 
  filter(!is.na(territoriality)) %>% 
  filter(territoriality != 2) # drop species with scores = 2 as there are very few
length(MUTI_Terr$territoriality)
View(MUTI_Terr)
#123 species with MUTI and territoriality


###### add and pair tree

# add rownames to data
row.names(MUTI_Terr) <- MUTI_Terr$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_Terr_phydat <- treedata(tree_out, MUTI_Terr, sort=T)

MUTI_Terr_phy <- MUTI_Terr_phydat$phy
MUTI_Terr_dat <- as.data.frame(MUTI_Terr_phydat$data)

str(MUTI_Terr_dat)
length(MUTI_Terr_dat$territoriality)
#123


### convert traits of interest to numeric
MUTI_Terr_dat$MUTIscore <- as.numeric(MUTI_Terr_dat$MUTIscore)
MUTI_Terr_dat$Mass_log <- as.numeric(MUTI_Terr_dat$Mass_log)
MUTI_Terr_dat$territoriality <- as.numeric(MUTI_Terr_dat$territoriality)


# Run phylogenetic linear model
MUTI_GLS_territory <- gls(MUTIscore~ territoriality + Mass_log, data = MUTI_Terr_dat, 
                         correlation = corPagel(0.5, phy = MUTI_Terr_phy, fixed=F, form = ~Species_Jetz), 
                         method = "ML") 

# model summary and results
summary(MUTI_GLS_territory) 
confint(MUTI_GLS_territory)

# model diagnostics
check_model(MUTI_GLS_territory) 
qqnorm(resid(MUTI_GLS_territory)) 
qqline(resid(MUTI_GLS_territory)) 
hist(resid(MUTI_GLS_territory)) 

# save model for easy retrieval 
saveRDS(MUTI_GLS_territory, here("Models/MUTI", "MUTI_GLS_territory.rds"))


######################## UN and Territoriality ##########################

# create a new data frame that contains only species with both UN and territoriality values
UN_Terr <- C_Social_dat2 %>% filter(!is.na(Urban)) %>%
  filter(!is.na(territoriality)) %>% 
  filter(territoriality != 2) %>% # drop species with scores = 2 as there are very few
  column_to_rownames(., var="Species_Jetz")
length(UN_Terr$territoriality)
#122 species with UN and territoriality

###### add and pair tree
tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UN_Terr_phydat <- treedata(tree_out, UN_Terr, sort=T)

UN_Terr_phy <- UN_Terr_phydat$phy
UN_Terr_dat <- as.data.frame(UN_Terr_phydat$data)

str(UN_Terr_dat)
length(UN_Terr_dat$territoriality)
# 122


### convert traits of interest to numeric
UN_Terr_dat$Urban <- as.numeric(UN_Terr_dat$Urban)
UN_Terr_dat$Mass_log <- as.numeric(UN_Terr_dat$Mass_log)
UN_Terr_dat$territoriality <- as.numeric(UN_Terr_dat$territoriality)

# Run the model using phyloglm(), which performs a logistic phylogenetic model to account for binary UN index
# default method ="logistic_MPLE"
# we will also scale and center the response variable to help with convergence
set.seed(124)
phyglm_UN_territorial_scale <- phyloglm( Urban ~ scale(territoriality) + scale(Mass_log), 
                                         data = UN_Terr_dat, 
                                         phy = UN_Terr_phy, 
                                         boot = 1000) 

summary(phyglm_UN_territorial_scale)

# save model for easy retrieval
saveRDS(phyglm_UN_territorial_scale, here("Models/UN", "phyglm_UN_territorial_scale.rds"))
# load model
phyglm_UN_territorial_scale <- readRDS(here("Models/UN", "phyglm_UN_territorial_scale.rds"))

# get alpha, t, and half life for the model
(phyglm_UN_territorial_scale$mean.tip.height) # t
(alpha_terr <- phyglm_UN_territorial_scale$alpha) # alpha
(hl_terr <- log(2)/alpha_terr) # half life
# small half life compared to t -> low phylogenetic signal



#############################################################################
#############################################################################

######################## UAI and Cooperative Breeding #######################

# create a new data frame that contains only species with both UAI and cooperative breeding values
UAI_CoopB <- C_Social_dat2 %>% filter(!is.na(aveUAI)) %>% filter(!is.na(cooperative)) 
length(UAI_CoopB$cooperative)
#766 species with UAI and cooperative


###### add and pair tree

# add rownames to data
row.names(UAI_CoopB) <- UAI_CoopB$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

UAI_CoopB_phydat <- treedata(tree_out, UAI_CoopB, sort=T)

UAI_CoopB_phy <- UAI_CoopB_phydat$phy
UAI_CoopB_dat <- as.data.frame(UAI_CoopB_phydat$data)

str(UAI_CoopB_dat)
length(UAI_CoopB_dat$cooperative)
#766


### convert traits of interest to numeric
UAI_CoopB_dat$aveUAI <- as.numeric(UAI_CoopB_dat$aveUAI)
UAI_CoopB_dat$Mass_log <- as.numeric(UAI_CoopB_dat$Mass_log)
UAI_CoopB_dat$cooperative <- as.numeric(UAI_CoopB_dat$cooperative)


# Run phylogenetic linear model
UAI_GLS_cooperative <- gls(aveUAI~ cooperative + Mass_log, data = UAI_CoopB_dat, 
                          correlation = corPagel(0.5, phy = UAI_CoopB_phy, fixed=F, form = ~Species_Jetz), 
                          method = "ML") 

# model summary and results
summary(UAI_GLS_cooperative) 
confint(UAI_GLS_cooperative)

# model diagnostics
check_model(UAI_GLS_cooperative) 
qqnorm(resid(UAI_GLS_cooperative)) 
qqline(resid(UAI_GLS_cooperative)) 
hist(resid(UAI_GLS_cooperative)) 

# save model for easy retrieval 
saveRDS(UAI_GLS_cooperative, here("Models/UAI", "UAI_GLS_cooperative.rds"))


######################## MUTI and Cooperative Breeding ##########################

# create a new data frame that contains only species with both MUTI and cooperative breeding values
MUTI_CoopB <- C_Social_dat2 %>% filter(!is.na(MUTIscore)) %>% filter(!is.na(cooperative)) 
length(MUTI_CoopB$cooperative)
#127 species with MUTI and cooperative breeding

###### add and pair tree

# add rownames to data
row.names(MUTI_CoopB) <- MUTI_CoopB$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

MUTI_CoopB_phydat <- treedata(tree_out, MUTI_CoopB, sort=T)

MUTI_CoopB_phy <- MUTI_CoopB_phydat$phy
MUTI_CoopB_dat <- as.data.frame(MUTI_CoopB_phydat$data)

str(MUTI_CoopB_dat)
length(MUTI_CoopB_dat$cooperative)
#127


### convert traits of interest to numeric
MUTI_CoopB_dat$MUTIscore <- as.numeric(MUTI_CoopB_dat$MUTIscore)
MUTI_CoopB_dat$Mass_log <- as.numeric(MUTI_CoopB_dat$Mass_log)
MUTI_CoopB_dat$cooperative <- as.numeric(MUTI_CoopB_dat$cooperative)


# Run phylogenetic linear model
MUTI_GLS_cooperative <- gls(MUTIscore~ cooperative + Mass_log, data = MUTI_CoopB_dat, 
                           correlation = corPagel(0.5, phy = MUTI_CoopB_phy, fixed=F, form = ~Species_Jetz), 
                           method = "ML") 

# model summary and results
summary(MUTI_GLS_cooperative) 
confint(MUTI_GLS_cooperative)

# model diagnostics
check_model(MUTI_GLS_cooperative) 
qqnorm(resid(MUTI_GLS_cooperative)) 
qqline(resid(MUTI_GLS_cooperative)) 
hist(resid(MUTI_GLS_cooperative))
 
# save model for easy retrieval 
saveRDS(MUTI_GLS_cooperative, here("Models/MUTI", "MUTI_GLS_cooperative.rds"))


######################## UN and Cooperative Breeding ##########################

# create a new data frame that contains only species with both UN and cooperative breeding values
UN_CoopB <- C_Social_dat2 %>% filter(!is.na(Urban)) %>% 
  filter(!is.na(cooperative)) %>%
  column_to_rownames(., var="Species_Jetz")
length(UN_CoopB$cooperative)
#129 species with UN and cooperative breeding

UN_CoopB %>% filter(cooperative==0) %>% nrow()
UN_CoopB %>% filter(cooperative==1) %>% nrow()
# There are only 8 species with cooperative = 1, all others have cooperative = 0
# We will not run a model here as this is highly unbalanced
