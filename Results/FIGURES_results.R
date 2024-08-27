#the purpose of this script is to create results figures for significant (or notable) models. 
#priotizing models where results for a predictor trait are supported across multiple indexes, or differ across multiple indexes



library(ape)
library(geiger)
library(nlme)
library(effects)
library(ggplot2)
library(ggeffects)
library(sjPlot)
library(dplyr)
if(!require(prediction)){
  install.packages("prediction")
  require(prediction)
}
library(prediction)


if(!require(effects)){
  install.packages("effects")
  require(effects)
}

if(!require(ggeffects)){
  install.packages("ggeffects")
  require(ggeffects)
}

if(!require(sjPlot)){
  install.packages("sjPlot")
  require(sjPlot)
}

if(!require(performance)){
  install.packages("performance")
  require(performance)
}

if(!require(here)){
  install.packages("here")
  require(here)
}

if(!require(phylolm)){
  install.packages("phylolm")
  require(phylolm)
}

if(!require(gridExtra)){
  install.packages("gridExtra")
  require(gridExtra)
}
library(gridExtra)

if(!require(see)){
  install.packages("see")
  require(see)
}



########################## BODY MASS ########################## 


###################### Prep

#load in "Coastal_Species_w_Mass.rds" - since this contains coastal species and all diet trait variables :) 

C_mass_dat <- readRDS(here("Outputs", "Coastal_Species_w_Mass.rds"))
str(C_mass_dat)

C_mass_dat2 <- C_mass_dat %>%
  mutate(Species_Jetz  = str_replace(Species_Jetz, " ", "_"))
str(C_mass_dat2)

C_mass_dat2$Urban <- ifelse(C_mass_dat2$Urban == "U", 1, 0)
View(C_mass_dat2)
colnames(C_mass_dat2)



######################## UAI and body mass ##########################

# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UAIDataUT <- C_mass_dat2 %>% filter(!is.na(aveUAI)) 
MassData1 <- UAIDataUT %>% filter(!is.na(Mass_log)) 
length(MassData1$Mass_log)
#798 species with UAI and Mass_log

colnames(MassData1)

###### add and pair tree

# add rownames to data
row.names(MassData1) <- MassData1$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Massphydat1 <- treedata(tree_out,MassData1, sort=T)

Massphy1 <- Massphydat1$phy
MassTraitDat1 <- as.data.frame(Massphydat1$data)

str(MassTraitDat1)
length(MassTraitDat1$Mass_log)
#798


### convert traits of interest to numeric

MassTraitDat1$aveUAI <- as.numeric(MassTraitDat1$aveUAI)
MassTraitDat1$Mass_log <- as.numeric(MassTraitDat1$Mass_log)


#lets run the model!


UAI_GLS_mass <- gls(aveUAI~ Mass_log, data = MassTraitDat1, 
                    correlation = corPagel(0.5, phy=Massphy1,fixed=F, form = ~Species_Jetz), 
                    method = "ML") 
#check out the model
check_model(UAI_GLS_mass) 
qqnorm(resid(UAI_GLS_mass)) 
qqline(resid(UAI_GLS_mass)) 
hist(resid(UAI_GLS_mass)) 


summary(UAI_GLS_mass) 
confint(UAI_GLS_mass)

############
#plot it 

library(effects)
Mass.UAIdb <- predictorEffect("Mass_log" , UAI_GLS_mass)


plot(Mass.UAIdb, ask =FALSE, xlab = "Body Mass", ylab = "UAI", main ="", lines=list(multiline=TRUE, col=c("indianred4")), confint=list(style="auto"), auto.key=FALSE, ylim=c(-20,10))

MASS.UAI.DF<-data.frame(Mass.UAIdb)

BodyMass_UAI_plot<-ggplot(data=MASS.UAI.DF,aes(x=Mass_log, y=fit)) +
  geom_line(color="indianred3",lwd=1.5) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="indianred3",alpha=.2, lwd=.1)+xlim(-1000,1000) +
  
  geom_point(data=MassTraitDat1,aes(x=jitter(Mass_log, 1), y =aveUAI),color="blue", bg="blue",alpha=.2, size=2,pch=21) +
  coord_cartesian(ylim = c(-0, 4.2), xlim =c(0,10.5)) + theme_classic() +
  theme(axis.text.x =  element_text(color="black", size = 13),axis.text.y =  element_text(color="black", size = 13), axis.title.x = element_text(size =14), axis.title.y = element_text(size =14)) +
  xlab("log(Body Mass)") + 
  ylab("Average UAI") 

BodyMass_UAI_plot

saveRDS(BodyMass_MUTI_PLOT, here("Results", "BodyMass_UAI.rds"))




######################## MUTI and body mass ##########################

# lets first simplify a NEW database by removing records where we don't have an MUTI / brood_value
MUTIDataUT <- C_mass_dat2 %>% filter(!is.na(MUTIscore)) 
MassData2 <- MUTIDataUT %>% filter(!is.na(Mass_log)) 
length(MassData2$Mass_log)
#130 species with MUTI and Mass_log

colnames(MassData2)

###### add and pair tree

# add rownames to data
row.names(MassData2) <- MassData2$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Massphydat2 <- treedata(tree_out,MassData2, sort=T)

Massphy2 <- Massphydat2$phy
MassTraitDat2 <- as.data.frame(Massphydat2$data)

str(MassTraitDat2)
length(MassTraitDat2$Mass_log)
#130


### convert traits of interest to numeric

MassTraitDat2$MUTIscore <- as.numeric(MassTraitDat2$MUTIscore)
MassTraitDat2$Mass_log <- as.numeric(MassTraitDat2$Mass_log)


#lets run the model!


MUTI_GLS_mass <- gls(MUTIscore~ Mass_log, data = MassTraitDat2, 
                     correlation = corPagel(0.5, phy=Massphy2,fixed=F, form = ~Species_Jetz), 
                     method = "ML") 
#check out the model
check_model(MUTI_GLS_mass) 
qqnorm(resid(MUTI_GLS_mass)) 
qqline(resid(MUTI_GLS_mass)) 
hist(resid(MUTI_GLS_mass)) 


summary(MUTI_GLS_mass) 
confint(MUTI_GLS_mass)


############
#plot it 

library(effects)
Mass.MUTIdb <- predictorEffect("Mass_log" , MUTI_GLS_mass)


plot(Mass.MUTIdb, ask =FALSE, xlab = "Body Mass", ylab = "MUTI", main ="", lines=list(multiline=TRUE, col=c("indianred4")), confint=list(style="auto"), auto.key=FALSE, ylim=c(-20,10))

MASS.MUTI.DF<-data.frame(Mass.MUTIdb)

BodyMass_MUTI_PLOT <-ggplot(data=MASS.MUTI.DF,aes(x=Mass_log, y=fit)) +
  geom_line(color="indianred3",lwd=1.5) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="indianred3",alpha=.2, lwd=.1)+xlim(-1000,1000) +
  
  geom_point(data= MassTraitDat2, aes(x=jitter(Mass_log, 1), y =MUTIscore),color="blue", bg="blue",alpha=.2, size=2,pch=21) +
  coord_cartesian(ylim = c(-4, 4), xlim =c(0,10.5)) + theme_classic() +
  theme(axis.text.x =  element_text(color="black", size = 13),axis.text.y =  element_text(color="black", size = 13), axis.title.x = element_text(size =14), axis.title.y = element_text(size =14)) +
  xlab("log(Body Mass)") + 
  ylab("MUTI") 

BodyMass_MUTI_PLOT

saveRDS(BodyMass_MUTI_PLOT, here("Results", "BodyMass_MUTI.rds"))









##################################################################
########################## LIFE HISTORY ########################## 


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




###############This is the correct code for boxplot 



if(!require(crayon)){
  install.packages("crayon")
  require(crayon)
}
library(crayon)
#plot it 

# Extract predicted values
#predictions <- predict(MUTI_GLS_develop, newdata = LifehistTraitDat11, type = "response")

# Combine predictions with original data
#LifehistTraitDat11$predicted_MUTI <- predictions

# Filter data for category 0 and category 1
category0_data <- subset(LifehistTraitDat11, developmental_mode == 0)
category1_data <- subset(LifehistTraitDat11, developmental_mode == 1)

View(LifehistTraitDat11)

# Plotting boxplots side by side

LifehistTraitDat11$developmental_mode <- as.numeric(as.character(LifehistTraitDat11$developmental_mode))

# now, plot! 

correct_developMUTI.plot <- ggplot() +
  geom_boxplot(data = category0_data, aes(x = "Precocial", y = MUTIscore), fill = "lightblue") +
  geom_boxplot(data = category1_data, aes(x = "Altricial", y = MUTIscore), fill = "lightgreen") +
  geom_point(data = LifehistTraitDat11, aes(x = ifelse(developmental_mode == 0, "Precocial", "Altricial"), y = MUTIscore), color = "deepskyblue4", size = 2, shape = 21, fill = "deepskyblue4", alpha = 0.3, position = position_jitter(width = 0.2)) +
  theme_classic() +
  xlab("Developmental Mode") +
  ylab("MUTIscore") +
  scale_x_discrete(labels = c("Precocial", "Altricial"))

print(correct_developMUTI.plot)

saveRDS(correct_developMUTI.plot, here("Results", "MUTI_developmental.rds"))





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




###############This is the correct code for boxplot 



if(!require(crayon)){
  install.packages("crayon")
  require(crayon)
}
library(crayon)
#plot it 



# Filter data for category 0 and category 1
category0_data <- subset(LifehistTraitDat12, developmental_mode == 0)
category1_data <- subset(LifehistTraitDat12, developmental_mode == 1)

# Plotting boxplots side by side

LifehistTraitDat12$developmental_mode <- as.numeric(as.character(LifehistTraitDat12$developmental_mode))

# now, plot! 

correct_develop_UN_plot <- ggplot() +
  geom_boxplot(data = category0_data, aes(x = "Precocial", y = Urban), fill = "lightblue") +
  geom_boxplot(data = category1_data, aes(x = "Altricial", y = Urban), fill = "lightgreen") +
  geom_point(data = LifehistTraitDat12, aes(x = ifelse(developmental_mode == 0, "Precocial", "Altricial"), y = Urban), color = "deepskyblue4", size = 2, shape = 21, fill = "deepskyblue4", alpha = 0.3, position = position_jitter(width = 0.2)) +
  theme_classic() +
  xlab("Developmental Mode") +
  ylab("MUTIscore") +
  scale_x_discrete(labels = c("Precocial", "Altricial"))

print(correct_develop_UN_plot)

saveRDS(correct_develop_UN_plot, here("Results", "UN_developmental.rds"))







