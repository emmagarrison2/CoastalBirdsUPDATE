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

saveRDS(BodyMass_UAI_plot, here("Results", "BodyMass_UAI.rds"))




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



##############################################################################
#arrange all body mass in a grid (besides two opposing clutch size)

UAI_bodymass_plot <- readRDS(here("Results", "BodyMass_UAI.rds"))
MUTI_bodymass_plot <- readRDS(here("Results", "BodyMass_MUTI.rds"))

bodymass_plots <- grid.arrange(UAI_bodymass_plot, MUTI_bodymass_plot, ncol=2, nrow = 1)




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
  ylab("Urban") +
  scale_x_discrete(labels = c("Precocial", "Altricial"))

print(correct_develop_UN_plot)

saveRDS(correct_develop_UN_plot, here("Results", "UN_developmental.rds"))



##############################################################################
#arrange all developmental mode plots in a grid

MUTI_develop_plot <- readRDS(here("Results", "MUTI_developmental.rds"))
UN_develop_plot <- readRDS(here("Results", "UN_developmental.rds"))

devepop_mode_plots <- grid.arrange(MUTI_develop_plot, UN_develop_plot, ncol=2, nrow = 1)





##################################################################
########################## NESTING TRAITS ########################## 




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


#let's plot the UAI & Nest site Low model that DOESN'T exclude ambiguous-nesting species for now. We can redo [Nest Site Low (ONLY)] plot later 

###############This is the correct code for boxplot 



if(!require(crayon)){
  install.packages("crayon")
  require(crayon)
}
library(crayon)
#plot it 



# Filter data for category 0 and category 1
category0_data <- subset(NestTraitDat4, NestSite_Low == 0)
#0 = not low
category1_data <- subset(NestTraitDat4, NestSite_Low == 1)
#1 = low 

# Plotting boxplots side by side

NestTraitDat4$NestSite_Low <- as.numeric(as.character(NestTraitDat4$NestSite_Low))

# now, plot! 

correct_low_UAI_plot <- ggplot() +
  geom_boxplot(data = category0_data, aes(x = "Not Low", y = aveUAI), fill = "lightblue") +
  geom_boxplot(data = category1_data, aes(x = "Low", y = aveUAI), fill = "lightgreen") +
  geom_point(data = NestTraitDat4, aes(x = ifelse(NestSite_Low == 0, "Not Low", "Low"), y = aveUAI), color = "deepskyblue4", size = 2, shape = 21, fill = "deepskyblue4", alpha = 0.3, position = position_jitter(width = 0.2)) +
  theme_classic() +
  xlab("Nest Site") +
  ylab("aveUAI") +
  scale_x_discrete(labels = c("Not Low", "Low"))

print(correct_low_UAI_plot)

saveRDS(correct_low_UAI_plot, here("Results", "UAI_low_nest.rds"))





######################## MUTI and %  Nest Site LOW ##########################
# 0 = not low
# 1 = LOW

# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
MUTIDataUT <- C_Nest_dat2 %>% filter(!is.na(MUTIscore)) 
NestData5 <- MUTIDataUT %>% filter(!is.na(NestSite_Low)) 
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



#let's plot the MUTI & Nest site Low model that DOESN'T exclude ambiguous-nesting species for now. We can redo [Nest Site Low (ONLY)] plot later 

###############This is the correct code for boxplot 



if(!require(crayon)){
  install.packages("crayon")
  require(crayon)
}
library(crayon)
#plot it 



# Filter data for category 0 and category 1
category0_data <- subset(NestTraitDat5, NestSite_Low == 0)
#0 = not low
category1_data <- subset(NestTraitDat5, NestSite_Low == 1)
#1 = low 

# Plotting boxplots side by side

NestTraitDat5$NestSite_Low <- as.numeric(as.character(NestTraitDat5$NestSite_Low))

# now, plot! 

correct_low_MUTI_plot <- ggplot() +
  geom_boxplot(data = category0_data, aes(x = "Not Low", y = MUTIscore), fill = "lightblue") +
  geom_boxplot(data = category1_data, aes(x = "Low", y = MUTIscore), fill = "lightgreen") +
  geom_point(data = NestTraitDat5, aes(x = ifelse(NestSite_Low == 0, "Not Low", "Low"), y = MUTIscore), color = "deepskyblue4", size = 2, shape = 21, fill = "deepskyblue4", alpha = 0.3, position = position_jitter(width = 0.2)) +
  theme_classic() +
  xlab("Nest Site") +
  ylab("MUTI score") +
  scale_x_discrete(labels = c("Not Low", "Low"))

print(correct_low_MUTI_plot)

saveRDS(correct_low_MUTI_plot, here("Results", "MUTI_low_nest.rds"))




##############################################################################
#arrange all nest site low plots in a grid

UAI_low_plot <- readRDS(here("Results", "UAI_low_nest.rds"))
MUTI_low_plot <- readRDS(here("Results", "MUTI_low_nest.rds"))

low_nesting_plots <- grid.arrange(UAI_low_plot, MUTI_low_plot, ncol=2, nrow = 1)


#

# now time for nest site high 

#


######################## UAI and %  Nest Site HIGH ##########################
# 0 = not high
# 1 = HIGH

# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UAIDataUT <- C_Nest_dat2 %>% filter(!is.na(aveUAI)) 
NestData7 <- UAIDataUT %>% filter(!is.na(NestSite_High)) 
length(NestData7$NestSite_High)
#796 species with UAI and NestSite_High

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
#796


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

#let's plot the UAI & Nest site High model that DOESN'T exclude ambiguous-nesting species for now. We can redo [Nest Site Low (ONLY)] plot later 

###############This is the correct code for boxplot 



if(!require(crayon)){
  install.packages("crayon")
  require(crayon)
}
library(crayon)
#plot it 



# Filter data for category 0 and category 1
category0_data <- subset(NestTraitDat7, NestSite_High == 0)
#0 = not high
category1_data <- subset(NestTraitDat7, NestSite_High == 1)
#1 = high

# Plotting boxplots side by side

NestTraitDat7$NestSite_High <- as.numeric(as.character(NestTraitDat7$NestSite_High))

# now, plot! 

correct_high_UAI_plot <- ggplot() +
  geom_boxplot(data = category0_data, aes(x = "Not High", y = aveUAI), fill = "lightblue") +
  geom_boxplot(data = category1_data, aes(x = "High", y = aveUAI), fill = "lightgreen") +
  geom_point(data = NestTraitDat7, aes(x = ifelse(NestSite_High == 0, "Not High", "High"), y = aveUAI), color = "deepskyblue4", size = 2, shape = 21, fill = "deepskyblue4", alpha = 0.3, position = position_jitter(width = 0.2)) +
  theme_classic() +
  xlab("Nest Site") +
  ylab("MUTI score") +
  scale_x_discrete(labels = c("Not High", "High"))

print(correct_high_UAI_plot)

saveRDS(correct_high_UAI_plot, here("Results", "UAI_high_nest.rds"))



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


MUTI_GLS_nest_high <- gls(MUTIscore~ NestSite_High + Mass_log, data = NestTraitDat8, 
                          correlation = corPagel(0.5, phy=Nestphy8,fixed=F, form = ~Species_Jetz), 
                          method = "ML") 
#check out the model
check_model(MUTI_GLS_nest_high) 
qqnorm(resid(MUTI_GLS_nest_high)) 
qqline(resid(MUTI_GLS_nest_high))
hist(resid(MUTI_GLS_nest_high))


summary(MUTI_GLS_nest_high) 
confint(MUTI_GLS_nest_high)


#let's plot the MUTI & Nest site High model that DOESN'T exclude ambiguous-nesting species for now. We can redo [Nest Site Low (ONLY)] plot later 

###############This is the correct code for boxplot 


if(!require(crayon)){
  install.packages("crayon")
  require(crayon)
}
library(crayon)
#plot it 



# Filter data for category 0 and category 1
category0_data <- subset(NestTraitDat8, NestSite_High == 0)
#0 = not high
category1_data <- subset(NestTraitDat8, NestSite_High == 1)
#1 = high

# Plotting boxplots side by side

NestTraitDat8$NestSite_High <- as.numeric(as.character(NestTraitDat8$NestSite_High))

# now, plot! 

correct_high_MUTI_plot <- ggplot() +
  geom_boxplot(data = category0_data, aes(x = "Not High", y = MUTIscore), fill = "lightblue") +
  geom_boxplot(data = category1_data, aes(x = "High", y = MUTIscore), fill = "lightgreen") +
  geom_point(data = NestTraitDat8, aes(x = ifelse(NestSite_High == 0, "Not High", "High"), y = MUTIscore), color = "deepskyblue4", size = 2, shape = 21, fill = "deepskyblue4", alpha = 0.3, position = position_jitter(width = 0.2)) +
  theme_classic() +
  xlab("Nest Site") +
  ylab("MUTI score") +
  scale_x_discrete(labels = c("Not High", "High"))

print(correct_high_MUTI_plot)

saveRDS(correct_high_MUTI_plot, here("Results", "MUTI_high_nest.rds"))




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
UN_M_nest_high <- phylolm(Urban~ NestSite_High + Mass_log, data=NestTraitDat9,
                            phy=Nestphy9, model="lambda") 

# time to check out the model 
qqnorm(UN_M_nest_high$residuals)
qqline(UN_M_nest_high$residuals) # what is happening? two separate lines bc of binomial... but is this the correct model check for binomial regression?
#the two lines do not have overlap... maybe this is good? 
hist(UN_M_nest_high$residuals, breaks = 20) 

#lets get those values for our results table 
summary(UN_M_nest_high)
confint(UN_M_nest_high)




#let's plot the UN & Nest site High model that DOESN'T exclude ambiguous-nesting species for now. We can redo [Nest Site Low (ONLY)] plot later 

###############This is the correct code for boxplot 


if(!require(crayon)){
  install.packages("crayon")
  require(crayon)
}
library(crayon)
#plot it 



# Filter data for category 0 and category 1
category0_data <- subset(NestTraitDat9, NestSite_High == 0)
#0 = not high
category1_data <- subset(NestTraitDat9, NestSite_High == 1)
#1 = high

# Plotting boxplots side by side

NestTraitDat9$NestSite_High <- as.numeric(as.character(NestTraitDat9$NestSite_High))

# now, plot! 

correct_high_UN_plot <- ggplot() +
  geom_boxplot(data = category0_data, aes(x = "Not High", y = Urban), fill = "lightblue") +
  geom_boxplot(data = category1_data, aes(x = "High", y = Urban), fill = "lightgreen") +
  geom_point(data = NestTraitDat9, aes(x = ifelse(NestSite_High == 0, "Not High", "High"), y = Urban), color = "deepskyblue4", size = 2, shape = 21, fill = "deepskyblue4", alpha = 0.3, position = position_jitter(width = 0.2)) +
  theme_classic() +
  xlab("Nest Site") +
  ylab("Urban") +
  scale_x_discrete(labels = c("Not High", "High"))

print(correct_high_UN_plot)

saveRDS(correct_high_UN_plot, here("Results", "UN_high_nest.rds"))

#

#

##############################################################################
#arrange all nest site low plots in a grid

UAI_high_plot <- readRDS(here("Results", "UAI_high_nest.rds"))
MUTI_high_plot <- readRDS(here("Results", "MUTI_high_nest.rds"))
UN_high_plot <- readRDS(here("Results", "UN_high_nest.rds"))

low_nesting_plots <- grid.arrange(UAI_high_plot, MUTI_high_plot, UN_high_plot, ncol=3, nrow = 1)






##################################################################
########################## SOCIAL TRAITS ########################## 



###################### Prep

#load in "Coastal_Species_Social.rds" - since this contains coastal species and all diet trait variables :) 

C_Social_dat <- readRDS(here("Outputs", "Coastal_Species_Social.rds"))
str(C_Social_dat)

C_Social_dat2 <- C_Social_dat %>%
  mutate(Species_Jetz  = str_replace(Species_Jetz, " ", "_"))
str(C_Social_dat2)

C_Social_dat2$Urban <- ifelse(C_Social_dat2$Urban == "U", 1, 0)
View(C_Social_dat2)
colnames(C_Social_dat2)



######################## UAI and territoriality ##########################

# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UAIDataUT <- C_Social_dat2 %>% filter(!is.na(aveUAI)) 
SocialData1 <- UAIDataUT %>% filter(!is.na(territoriality))
length(SocialData1$territoriality)
#766 species with UAI and territoriality

colnames(SocialData1)

###### add and pair tree

# add rownames to data
row.names(SocialData1) <- SocialData1$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Socialphydat1 <- treedata(tree_out,SocialData1, sort=T)

Socialphy1 <- Socialphydat1$phy
SocialTraitDat1 <- as.data.frame(Socialphydat1$data)

str(SocialTraitDat1)
length(SocialTraitDat1$territoriality)
#766


### convert traits of interest to numeric

SocialTraitDat1$aveUAI <- as.numeric(SocialTraitDat1$aveUAI)
SocialTraitDat1$Mass_log <- as.numeric(SocialTraitDat1$Mass_log)
SocialTraitDat1$territoriality <- as.numeric(SocialTraitDat1$territoriality)


#lets run the model!


UAI_GLS_territory <- gls(aveUAI~ territoriality + Mass_log, data = SocialTraitDat1, 
                         correlation = corPagel(0.5, phy=Socialphy1,fixed=F, form = ~Species_Jetz), 
                         method = "ML") 
#check out the model
check_model(UAI_GLS_territory) 
qqnorm(resid(UAI_GLS_territory)) 
qqline(resid(UAI_GLS_territory)) 
hist(resid(UAI_GLS_territory)) 


summary(UAI_GLS_territory) 
confint(UAI_GLS_territory)


#let's plot the UAI & Territoriality model

###############This is the correct code for boxplot 


if(!require(crayon)){
  install.packages("crayon")
  require(crayon)
}
library(crayon)
#plot it 



# Filter data for category 0 and category 1
category0_data <- subset(SocialTraitDat1, territoriality == 0)
#0 = not high
category1_data <- subset(SocialTraitDat1, territoriality == 1)
#1 = high
category2_data <- subset(SocialTraitDat1, territoriality == 2)

# Plotting boxplots side by side

SocialTraitDat1$territoriality <- as.numeric(as.character(SocialTraitDat1$territoriality))

# now, plot! 

correct_territorial_UAI_plot <- ggplot() +
  geom_boxplot(data = category0_data, aes(x = factor(0, labels = "Non-Territorial"), y = aveUAI), fill = "lightblue") +
  geom_boxplot(data = category1_data, aes(x = factor(1, labels = "Sesonally territorial"), y = aveUAI), fill = "lightgreen") +
  geom_boxplot(data = category2_data, aes(x = factor(2, labels = "Year-Round territorial"), y = aveUAI), fill = "orange") +
  geom_point(data = SocialTraitDat1, aes(x = factor(territoriality, levels = c(0, 1, 2), labels = c("Non-Territorial", "Sesonally territorial", "Year-Round territorial")), y = aveUAI), 
             color = "deepskyblue4", size = 2, shape = 21, fill = "deepskyblue4", alpha = 0.3, position = position_jitter(width = 0.2)) +
  theme_classic() +
  xlab("Territoriality") +
  ylab("aveUAI") +
  scale_x_discrete(labels = c("Non-Territorial", "Sesonally territorial", "Year-Round territorial"))

print(correct_territorial_UAI_plot)


saveRDS(correct_territorial_UAI_plot, here("Results", "UAI_territorial.rds"))







