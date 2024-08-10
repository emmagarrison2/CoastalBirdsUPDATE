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



