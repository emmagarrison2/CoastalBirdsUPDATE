# The objective of this script is to run phylogenetic linear models JUST body mass

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

L5<-ggplot(data=MASS.UAI.DF,aes(x=Mass_log, y=fit)) +
  geom_line(color="indianred3",lwd=1.5) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="indianred3",alpha=.2, lwd=.1)+xlim(-1000,1000) +
  
  geom_point(data=MassTraitDat1,aes(x=jitter(Mass_log, 1), y =aveUAI),color="blue", bg="blue",alpha=.2, size=1,pch=21) +
  coord_cartesian(ylim = c(-0, 4.2), xlim =c(0,10.5)) + theme_classic() +
  theme(axis.text.x =  element_text(color="black", size = 13),axis.text.y =  element_text(color="black", size = 13), axis.title.x = element_text(size =14), axis.title.y = element_text(size =14)) +
  xlab("log(Body Mass)") + 
  ylab("Average UAI") 

L5




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

L5<-ggplot(data=MASS.MUTI.DF,aes(x=Mass_log, y=fit)) +
  geom_line(color="indianred3",lwd=1.5) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="indianred3",alpha=.2, lwd=.1)+xlim(-1000,1000) +
  
  geom_point(data= MassTraitDat2, aes(x=jitter(Mass_log, 1), y =MUTIscore),color="blue", bg="blue",alpha=.2, size=1,pch=21) +
  coord_cartesian(ylim = c(-4, 4), xlim =c(0,10.5)) + theme_classic() +
  theme(axis.text.x =  element_text(color="black", size = 13),axis.text.y =  element_text(color="black", size = 13), axis.title.x = element_text(size =14), axis.title.y = element_text(size =14)) +
  xlab("log(Body Mass)") + 
  ylab("MUTI") 

L5



######################## UN and body mass ##########################

# lets first simplify a NEW database by removing records where we don't have an UN / brood_value
UNDataUT <- C_mass_dat2 %>% filter(!is.na(Urban)) 
MassData3 <- UNDataUT %>% filter(!is.na(Mass_log)) 
length(MassData3$Mass_log)
#129 species with UN and Mass_log

colnames(MassData3)

###### add and pair tree

# add rownames to data
row.names(MassData3) <- MassData3$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Massphydat3 <- treedata(tree_out,MassData3, sort=T)

Massphy3 <- Massphydat3$phy
MassTraitDat3 <- as.data.frame(Massphydat3$data)

str(MassTraitDat3)
length(MassTraitDat3$Mass_log)
#129


### convert traits of interest to numeric

MassTraitDat3$Urban <- as.numeric(MassTraitDat3$Urban)
MassTraitDat3$Mass_log <- as.numeric(MassTraitDat3$Mass_log)



#lets run the model using Phylolm!  

#(have to use lambda, until we figure out a way to make GLS work with binomial linear regression)
UN_M_mass <- phylolm(Urban~ Mass_log, data=MassTraitDat3,
                            phy=Massphy3, model="lambda") 

# time to check out the model 
qqnorm(UN_M_mass$residuals)
qqline(UN_M_mass$residuals)
hist(UN_M_mass$residuals, breaks = 12) 

#lets get those values for our results table 
summary(UN_M_mass)
confint(UN_M_mass)
