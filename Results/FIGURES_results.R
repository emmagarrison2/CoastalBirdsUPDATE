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

if(!require(stringr)){
  install.packages("stringr")
  require(stringr)
}
library(stringr)


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


plot(Mass.UAIdb, ask =FALSE, xlab = "Body Mass", ylab = "UAI", main ="", lines=list(multiline=TRUE, col=c("#FFD67B")), confint=list(style="auto"), auto.key=FALSE, ylim=c(-20,10))

MASS.UAI.DF<-data.frame(Mass.UAIdb)

BodyMass_UAI_plot<-ggplot(data=MASS.UAI.DF,aes(x=Mass_log, y=fit)) +
  geom_line(color="#FFD67B",lwd=1.5) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="#FFD67B",alpha=.3, lwd=.1)+xlim(-1000,1000) +
  
  geom_point(data=MassTraitDat1,aes(x=jitter(Mass_log, 1), y =aveUAI),color="#E48816", bg="#E48816",alpha=.6, size=1,pch=21) +
  coord_cartesian(ylim = c(-0, 4.2), xlim =c(0,10.5)) + theme_classic() +
  theme(axis.text.x =  element_text(color="black", size = 10),
        axis.text.y =  element_text(color="black", size = 10), 
        axis.title.x = element_text(size =14), 
        axis.title.y = element_text(size =14)) +
  xlab("log(Body Mass)") + 
  ylab("UAI") 

BodyMass_UAI_plot

saveRDS(BodyMass_UAI_plot, here("Results", "FIG_BodyMass_UAI.rds"))





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


######################## UAI and % brood value ##########################


# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UAIDataUT <- C_LifeHist_dat2 %>% filter(!is.na(aveUAI)) 
LifehistData1 <- UAIDataUT %>% filter(!is.na(brood_value)) 
length(LifehistData1$brood_value)
#480 species with UAI and brood_value

###### add and pair tree

# add rownames to data
row.names(LifehistData1) <- LifehistData1$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Lifehistphydat1 <- treedata(tree_out,LifehistData1, sort=T)

Lifehistphy1 <- Lifehistphydat1$phy
LifehistTraitDat1 <- as.data.frame(Lifehistphydat1$data)

str(LifehistTraitDat1)
length(LifehistTraitDat1$brood_value)
#480

### convert traits of interest to numeric

LifehistTraitDat1$aveUAI <- as.numeric(LifehistTraitDat1$aveUAI)
LifehistTraitDat1$Mass_log <- as.numeric(LifehistTraitDat1$Mass_log)
LifehistTraitDat1$brood_value <- as.numeric(LifehistTraitDat1$brood_value)


#lets run the model!


UAI_GLS_bv <- gls(aveUAI~ brood_value + Mass_log, data = LifehistTraitDat1, 
                  correlation = corPagel(0.5, phy=Lifehistphy1,fixed=F, form = ~Species_Jetz), 
                  method = "ML") 
#check out the model
check_model(UAI_GLS_bv) ## low collinearity - which is good! - normality of residuals line does not fall on line, but is in a straight line 
qqnorm(resid(UAI_GLS_bv)) 
qqline(resid(UAI_GLS_bv)) #most points fall on the line 
hist(resid(UAI_GLS_bv)) #roughly normal dist of residuals


summary(UAI_GLS_bv) # strong phylogenetic relationship between traits and response, but no remaining influence of any predictor trait and response
confint(UAI_GLS_bv)



############
#plot it 

library(effects)
bv.UAIdb <- predictorEffect("brood_value" , UAI_GLS_bv)


plot(bv.UAIdb, ask =FALSE, xlab = "Brood Value", ylab = "UAI", main ="", lines=list(multiline=TRUE, col=c("#FFD67B")), confint=list(style="auto"), auto.key=FALSE, ylim=c(-20,10))

BV.UAI.DF<-data.frame(bv.UAIdb)

BV_UAI_plot<-ggplot(data=BV.UAI.DF,aes(x=brood_value, y=fit)) +
  geom_line(color="#FFD67B",lwd=1.5) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="#FFD67B",alpha=.3, lwd=.1)+xlim(-1000,1000) +
  
  geom_point(data=LifehistTraitDat1,aes(x=jitter(brood_value, 1), y =aveUAI),color="#E48816", bg="#E48816",alpha=.6, size=2,pch=21) +
  coord_cartesian(ylim = c(-0, 3.5), xlim =c(-6,-1.7)) + theme_classic() +
  theme(axis.text.x =  element_text(color="black", size = 10),
        axis.text.y =  element_text(color="black", size = 10), 
        axis.title.x = element_text(size =14), 
        axis.title.y = element_text(size =14)) +
  xlab("Brood Value") + 
  ylab("UAI") 

BV_UAI_plot

saveRDS(BV_UAI_plot, here("Results", "FIG_BroodValue_UAI_outlier.rds"))


#################  WITHOUT OUTLIER #################

################ so, brood value has this outlier... let's try again with UAI and Brood Value to see if this relationship is still significantly negative 

LifehistDataFiltered <- LifehistData1 %>% filter(brood_value >= -5)

# Add rownames to filtered data
row.names(LifehistDataFiltered) <- LifehistDataFiltered$Species_Jetz

# Pair with the tree
LifehistphydatFiltered <- treedata(tree_out, LifehistDataFiltered, sort=T)

LifehistphyFiltered <- LifehistphydatFiltered$phy
LifehistTraitDatFiltered <- as.data.frame(LifehistphydatFiltered$data)

# Convert traits of interest to numeric
LifehistTraitDatFiltered$aveUAI <- as.numeric(LifehistTraitDatFiltered$aveUAI)
LifehistTraitDatFiltered$Mass_log <- as.numeric(LifehistTraitDatFiltered$Mass_log)
LifehistTraitDatFiltered$brood_value <- as.numeric(LifehistTraitDatFiltered$brood_value)

# Run the GLS model on the filtered dataset
UAI_GLS_bv_filtered <- gls(aveUAI ~ brood_value + Mass_log, data = LifehistTraitDatFiltered, 
                           correlation = corPagel(0.5, phy=LifehistphyFiltered, fixed=F, form = ~Species_Jetz), 
                           method = "ML")

# Check the model
check_model(UAI_GLS_bv_filtered)
qqnorm(resid(UAI_GLS_bv_filtered)) 
qqline(resid(UAI_GLS_bv_filtered))
hist(resid(UAI_GLS_bv_filtered))

# Model summary and confidence intervals
summary(UAI_GLS_bv_filtered) # still significant (p = 0.0288)
confint(UAI_GLS_bv_filtered) # still significant (C.I. = (-0.349, -0.0195)

############


#plot it 

library(effects)
bv.filtered.UAIdb <- predictorEffect("brood_value" , UAI_GLS_bv_filtered)


plot(bv.filtered.UAIdb, ask =FALSE, xlab = "Brood Value", ylab = "UAI", main ="", lines=list(multiline=TRUE, col=c("#FFD67B")), confint=list(style="auto"), auto.key=FALSE, ylim=c(-20,10))

BV.filtered.UAI.DF<-data.frame(bv.filtered.UAIdb)

BV_no_outlier_UAI_plot<-ggplot(data=BV.filtered.UAI.DF,aes(x=brood_value, y=fit)) +
  geom_line(color="#FFD67B",lwd=1.5) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="#FFD67B",alpha=.3, lwd=.1)+xlim(-1000,1000) +
  
  geom_point(data=LifehistTraitDatFiltered,aes(x=jitter(brood_value, 1), y =aveUAI),color="#E48816", bg="#E48816",alpha=.6, size=2,pch=21) +
  coord_cartesian(ylim = c(-0, 3.5), xlim =c(-4.3,-1.8)) + theme_classic() +
  theme(axis.text.x =  element_text(color="black", size = 10),
        axis.text.y =  element_text(color="black", size = 10), 
        axis.title.x = element_text(size =14, margin = margin(t = 8)), 
        axis.title.y = element_text(size =14)) +
  xlab("Brood Value") + 
  ylab("UAI") 

BV_no_outlier_UAI_plot

saveRDS(BV_no_outlier_UAI_plot, here("Results", "FIG_BroodValue_UAI_no_outlier.rds"))





######################## UAI and % clutch size ##########################


# lets first simplify a NEW database by removing records where we don't have an UAI / brood_value
UAIDataUT <- C_LifeHist_dat2 %>% filter(!is.na(aveUAI)) 
LifehistData4 <- UAIDataUT %>% filter(!is.na(clutch_size)) 
length(LifehistData4$clutch_size)
#738 species with UAI and clutch_size

###### add and pair tree

# add rownames to data
LifehistData4 <- as.data.frame(LifehistData4)
row.names(LifehistData4) <- LifehistData4$Species_Jetz

tree_out<- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

Lifehistphydat4 <- treedata(tree_out,LifehistData4, sort=T)

Lifehistphy4 <- Lifehistphydat4$phy
LifehistTraitDat4 <- as.data.frame(Lifehistphydat4$data)

str(LifehistTraitDat4)
length(LifehistTraitDat4$clutch_size)
#738

### convert traits of interest to numeric

LifehistTraitDat4$aveUAI <- as.numeric(LifehistTraitDat4$aveUAI)
LifehistTraitDat4$Mass_log <- as.numeric(LifehistTraitDat4$Mass_log)
LifehistTraitDat4$clutch_size <- as.numeric(LifehistTraitDat4$clutch_size)



#lets run the model!


UAI_GLS_clutch <- gls(aveUAI~ clutch_size + Mass_log, data = LifehistTraitDat4, 
                      correlation = corPagel(0.5, phy=Lifehistphy4,fixed=F, form = ~Species_Jetz), 
                      method = "ML") 
#check out the model
check_model(UAI_GLS_clutch) ## low collinearity - which is good! - normality of residuals line does not fall on line, but is in a straight line 
qqnorm(resid(UAI_GLS_clutch)) 
qqline(resid(UAI_GLS_clutch)) #most points fall on the line 
hist(resid(UAI_GLS_clutch)) #roughly normal dist of residuals


summary(UAI_GLS_clutch) # strong phylogenetic relationship between traits and response, but no remaining influence of any predictor trait and response
confint(UAI_GLS_clutch)


############


#plot it 

library(effects)
clutch.UAIdb <- predictorEffect("clutch_size" , UAI_GLS_clutch)


plot(clutch.UAIdb, ask =FALSE, xlab = "Clutch Size", ylab = "UAI", main ="", lines=list(multiline=TRUE, col=c("#FFD67B")), confint=list(style="auto"), auto.key=FALSE, ylim=c(-20,10))

clutch.UAI.DF<-data.frame(clutch.UAIdb)

clutch_UAI_plot<-ggplot(data=clutch.UAI.DF,aes(x=clutch_size, y=fit)) +
  geom_line(color="#FFD67B",lwd=1.5) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="#FFD67B",alpha=.3, lwd=.1)+xlim(-1000,1000) +
  
  geom_point(data=LifehistTraitDat4,aes(x=jitter(clutch_size, 1), y =aveUAI),color="#E48816", bg="#E48816",alpha=.6, size=2,pch=21) +
  coord_cartesian(ylim = c(-0, 4.1), xlim =c(0,15.5)) + theme_classic() +
  theme(axis.text.x =  element_text(color="black", size = 10),
        axis.text.y =  element_text(color="black", size = 10), 
        axis.title.x = element_text(size =14, margin = margin(t = 8)), 
        axis.title.y = element_text(size =14)) +
  xlab("Clutch Size") + 
  ylab("UAI") 

clutch_UAI_plot


saveRDS(clutch_UAI_plot, here("Results", "FIG_clutch_UAI.rds"))




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



# Plotting boxplots side by side

LifehistTraitDat11$developmental_mode <- as.numeric(as.character(LifehistTraitDat11$developmental_mode))

# now, plot! 


correct_develop_MUTI_plot <- ggplot() +
  geom_boxplot(data = category0_data, aes(x = "Precocial", y = MUTIscore), fill = "olivedrab1", alpha = 0.4) +
  geom_boxplot(data = category1_data, aes(x = "Altricial", y = MUTIscore), fill = "olivedrab", alpha = 0.4) +
  geom_point(data = LifehistTraitDat11, aes(x = ifelse(developmental_mode == 0, "Precocial", "Altricial"), y = MUTIscore), color = "olivedrab4", size = 2, shape = 21, fill = "olivedrab4", alpha = 0.6, position = position_jitter(width = 0.2)) +
  theme_classic() +
  xlab("Developmental Mode") +
  ylab("MUTI") +
  theme(axis.text.x =  element_text(color="black", size = 10),
        axis.text.y =  element_text(color="black", size = 10), 
        axis.title.x = element_text(size =14, margin = margin(t = 8)), 
        axis.title.y = element_text(size =14)) +
  scale_x_discrete(limits = c("Precocial", "Altricial"), labels = c("Precocial", "Altricial"))  # Specify the order here

print(correct_develop_MUTI_plot)

saveRDS(correct_develop_MUTI_plot, here("Results", "FIG_Develop_MUTI.rds"))




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



################# Bar Chart 
View(LifehistTraitDat12)


Develop_UN <- LifehistTraitDat12 %>% 
  select(Urban, developmental_mode)
colnames(Develop_UN)
nrow(Develop_UN) #129, looks good! 


Develop_UN$developmental_mode <- factor(Develop_UN$developmental_mode, levels = c(0, 1), labels = c("Precocial", "Altricial"))

str(Develop_UN)
View(Develop_UN)

#stacked bar chart ---> IT WORKS!! 

#first, calculate percentages 

Develop_UN_summary <- Develop_UN %>%
  group_by(developmental_mode, Urban) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(developmental_mode) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  mutate(Group_Urban = paste(developmental_mode, Urban, sep = " & "))  # Combine Group and Urban for legend

View(Develop_UN_summary)

library(dplyr)

Develop_UN_plot <- ggplot(Develop_UN_summary, aes(x = developmental_mode, y = percentage, fill = Group_Urban)) + 
  geom_bar(stat = "identity", position = "fill", color = "black") + 
  scale_y_continuous(labels = scales::percent) +  # Show percentages on the y-axis
  # Custom colors for each Group & Urban combination
  scale_fill_manual(
    values = c(
      "Precocial & 0" = "#B8D5E9",   # Color for "Precocial & Nonurban"
      "Precocial & 1" = "#6C9CCC",        # Color for "Precocial & Urban"
      "Altricial & 0" = "#B8D5E9",  # Color for "Altricial & Nonurban"
      "Altricial & 1" = "#6C9CCC"        # Color for "Altricial & Urban"
    )
  ) + 
  labs(
    x = "Developmental Mode", 
    y = "UN", 
    fill = ""  # Legend title
  ) + 
  theme_classic() + 
  theme(
    axis.text.x = element_text(color="black", size = 10),
    axis.text.y = element_text(color="black", size = 10),
    axis.title.x = element_text(size = 14, margin = margin(t = 8)), 
    axis.title.y = element_text(size = 12), 
    legend.position = "none"
  ) + 
  annotate(
    "text", 
    x = c(1, 1, 2, 2), 
    y = c(0.15, 0.63, 0.21, 0.7), 
    label = c("Urban", "Non-Urban", "Urban", "Non-Urban"))



print(Develop_UN_plot)

saveRDS(Develop_UN_plot, here("Results", "FIG_Develop_UN.rds"))




##############################################################################
#plot all the life history results together... 

BV_UAI_plot_no_outlier <- readRDS(here("Results", "FIG_BroodValue_UAI_no_outlier.rds"))
clutch_UAI_plot <- readRDS(here("Results", "FIG_clutch_UAI.rds"))
develop_MUTI_plot <- readRDS(here("Results", "FIG_Develop_MUTI.rds"))
develop_UN_plot <- readRDS(here("Results", "FIG_Develop_UN.rds"))


#using patchwork package 
if(!require(patchwork)){
  install.packages("patchwork")
  require(patchwork)
}
library(patchwork)


all_life_history_plots <- BV_UAI_plot_no_outlier + clutch_UAI_plot + develop_MUTI_plot + develop_UN_plot + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 14))

print(all_life_history_plots)

saveRDS(all_life_history_plots, here("Results", "FIG_all_life_history.rds"))


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
not_low_data <- subset(NestTraitDat4, NestSite_Low == 0)
#0 = not high
low_data <- subset(NestTraitDat4, NestSite_Low == 1)
#1 = high

# Plotting boxplots side by side

NestTraitDat4$NestSite_Low <- as.numeric(as.character(NestTraitDat4$NestSite_Low))

# now, plot! 

correct_low_UAI_plot <- ggplot() +
  geom_boxplot(data = not_low_data, aes(x = "Not Low", y = aveUAI), fill = "#FFD67B", alpha = 0.4) +
  geom_boxplot(data = low_data, aes(x = "Low", y = aveUAI), fill = "#E48816", alpha = 0.4) +
  geom_point(data = NestTraitDat4, aes(x = ifelse(NestSite_Low == 0, "Not Low", "Low"), y = aveUAI), color = "#E48816", size = 1, shape = 21, fill = "#E48816", alpha = 0.7, position = position_jitter(width = 0.2)) +
  theme_classic() +
  xlab("Nest Site") +
  ylab("UAI") +
  scale_x_discrete(limits = c("Not Low", "Low"), labels = c("Not Low", "Low"))  + # Specify the order here 
  theme(
  axis.title.x = element_text(color="black", size = 14),    # Adjust font size for x-axis label
  axis.title.y = element_text(color="black", size = 14),    # Adjust font size for y-axis label
  axis.text.x = element_text(size = 10, margin = margin(t = 8)),     # Adjust font size for x-axis tick labels
  axis.text.y = element_text(size = 10)      # Adjust font size for y-axis tick labels
)

print(correct_low_UAI_plot)

saveRDS(correct_low_UAI_plot, here("Results", "FIG_low_nest_UAI.rds"))





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

# now, plot! olivedrab1, olivedrab, olivedrab4

correct_low_MUTI_plot <- ggplot() +
  geom_boxplot(data = category0_data, aes(x = "Not Low", y = MUTIscore), fill = "olivedrab1", alpha = 0.4) +
  geom_boxplot(data = category1_data, aes(x = "Low", y = MUTIscore), fill = "olivedrab", alpha = 0.4) +
  geom_point(data = NestTraitDat5, aes(x = ifelse(NestSite_Low == 0, "Not Low", "Low"), y = MUTIscore), color = "olivedrab4", size = 1, shape = 21, fill = "olivedrab4", alpha = 0.8, position = position_jitter(width = 0.2)) +
  theme_classic() +
  xlab("Nest Site") +
  ylab("MUTI") +
  scale_x_discrete(limits = c("Not Low", "Low"), labels = c("Not Low", "Low")) + 
  theme(
    axis.title.x = element_text(color="black", size = 14),    # Adjust font size for x-axis label
    axis.title.y = element_text(color="black", size = 14),    # Adjust font size for y-axis label
    axis.text.x = element_text(size = 10, margin = margin(t = 8)),     # Adjust font size for x-axis tick labels
    axis.text.y = element_text(size = 10)      # Adjust font size for y-axis tick labels
  )

print(correct_low_MUTI_plot)

saveRDS(correct_low_MUTI_plot, here("Results", "FIG_low_nest_MUTI.rds"))






##############################################################################
#arrange all nest site low plots in a grid

UAI_low_plot <- readRDS(here("Results", "FIG_low_nest_UAI.rds"))
MUTI_low_plot <- readRDS(here("Results", "FIG_low_nest_MUTI.rds"))

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
not_high_data <- subset(NestTraitDat7, NestSite_High == 0)
#0 = not high
high_data <- subset(NestTraitDat7, NestSite_High == 1)
#1 = high

# Plotting boxplots side by side

NestTraitDat7$NestSite_High <- as.numeric(as.character(NestTraitDat7$NestSite_High))

# now, plot! 

correct_high_UAI_plot <- ggplot() +
  geom_boxplot(data = not_high_data, aes(x = "Not High", y = aveUAI), fill = "#FFD67B", alpha = 0.4) +
  geom_boxplot(data = high_data, aes(x = "High", y = aveUAI), fill = "#E48816", alpha = 0.4) +
  geom_point(data = NestTraitDat7, aes(x = ifelse(NestSite_High == 0, "Not High", "High"), y = aveUAI), color = "#E48816", size = 1, shape = 21, fill = "#E48816", alpha = 0.7, position = position_jitter(width = 0.2)) +
  theme_classic() +
  xlab("Nest Site") +
  ylab("UAI") +
  scale_x_discrete(limits = c("Not High", "High"), labels = c("Not High", "High"))   + # Specify the order here
  theme(
    axis.title.x = element_text(size = 14),    # Adjust font size for x-axis label
    axis.title.y = element_text(size = 14),    # Adjust font size for y-axis label
    axis.text.x = element_text(size = 10, margin = margin(t = 8)),     # Adjust font size for x-axis tick labels
    axis.text.y = element_text(size = 10)      # Adjust font size for y-axis tick labels
  )

print(correct_high_UAI_plot)

saveRDS(correct_high_UAI_plot, here("Results", "FIG_high_nest_UAI.rds"))



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


#let's plot the MUTI & Nest site High model that DOESN'T exclude ambiguous-nesting species for now. We can redo [Nest Site High (ONLY)] plot later 

###############This is the correct code for boxplot 


if(!require(crayon)){
  install.packages("crayon")
  require(crayon)
}
library(crayon)
#plot it 



# Filter data for category 0 and category 1
not_high_data <- subset(NestTraitDat8, NestSite_High == 0)
#0 = not high
mean(not_high_data$MUTIscore)
# -0.261

high_data <- subset(NestTraitDat8, NestSite_High == 1)
#1 = high
mean(high_data$MUTIscore)
# 0.170 



# now, plot! olivedrab1, olivedrab, olivedrab4

correct_high_MUTI_plot <- ggplot() +
  geom_boxplot(data = not_high_data, aes(x = "Not High", y = MUTIscore), fill = "olivedrab1", alpha = 0.4) +
  geom_boxplot(data = high_data, aes(x = "High", y = MUTIscore), fill = "olivedrab", alpha = 0.4) +
  geom_point(data = NestTraitDat8, aes(x = ifelse(NestSite_High == 0, "Not High", "High"), y = MUTIscore), color = "olivedrab4", size = 1, shape = 21, fill = "olivedrab4", alpha = 0.8, position = position_jitter(width = 0.2)) +
  theme_classic() +
  xlab("Nest Site") +
  ylab("MUTI") +
  scale_x_discrete(limits = c("Not High", "High"), labels = c("Not High", "High"))   +
  theme(axis.title.x = element_text(size = 14),    # Adjust font size for x-axis label
        axis.title.y = element_text(size = 14),    # Adjust font size for y-axis label
        axis.text.x = element_text(size = 10, margin = margin(t = 8)),     # Adjust font size for x-axis tick labels
        axis.text.y = element_text(size = 10))
    
print(correct_high_MUTI_plot)


saveRDS(correct_high_MUTI_plot, here("Results", "FIG_high_nest_MUTI.rds"))




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


View(NestTraitDat9)
UN_Nest_High <- NestTraitDat9 %>% select(Urban, NestSite_High)
View(UN_Nest_High)


#PLOT -- let's use the one with all species (so, not exlcuding ambiguous nesters (ONLY))

#stacked bar chart ---> IT WORKS!! 


#relabel 0 = not high, and 1 = high 
UN_Nest_High$NestSite_High <- factor(UN_Nest_High$NestSite_High, levels = c(0, 1), labels = c("Not High", "High"))


#first, calculate percentages 
UN_Nest_High_Summary <- UN_Nest_High %>%
  group_by(NestSite_High, Urban) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(NestSite_High) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  mutate(Group_Urban = paste(NestSite_High, Urban, sep = " & "))  # Combine Group and Urban for legend


View(UN_Nest_High_Summary)

# light - B8D5E9, dark - 6C9CCC

library(dplyr)

Nest_High_UN_plot <- ggplot(UN_Nest_High_Summary, aes(x = NestSite_High, y = percentage, fill = Group_Urban)) + 
  geom_bar(stat = "identity", position = "fill", color = "black") + 
  scale_y_continuous(labels = scales::percent) +  # Show percentages on the y-axis
  # Custom colors for each Group & Urban combination
  scale_fill_manual(
    values = c(
      "Not High & 0" = "#B8D5E9",   # Color for "Precocial & Nonurban"
      "Not High & 1" = "#6C9CCC",        # Color for "Precocial & Urban"
      "High & 0" = "#B8D5E9",  # Color for "Altricial & Nonurban"
      "High & 1" = "#6C9CCC"        # Color for "Altricial & Urban"
    )
  ) + 
  labs(
    x = "Nest Site", 
    y = "UN", 
    fill = ""  # Legend title
  ) + 
  theme_classic() + 
  theme(
    axis.title.x = element_text(size = 14, margin = margin(t = 8)), 
    axis.title.y = element_text(size = 12), 
    axis.text.x = element_text(size = 10), 
    axis.text.y = element_text(size = 10), 
    legend.position = "none"
  ) + 
  annotate(
    "text", 
    x = c(1, 1, 2, 2), 
    y = c(0.12, 0.63, 0.21, 0.7), 
    label = c("Urban", "Non-Urban", "Urban", "Non-Urban"))

print(Nest_High_UN_plot)
 
saveRDS(Nest_High_UN_plot, here("Results", "FIG_high_nest_UN.rds"))



##############################################################################
#arrange all nest site high plots in a grid

UAI_high_plot <- readRDS(here("Results", "FIG_high_nest_UAI.rds"))
MUTI_high_plot <- readRDS(here("Results", "FIG_high_nest_MUTI.rds"))
UN_high_plot <- readRDS(here("Results", "FIG_high_nest_UN.rds"))

low_nesting_plots <- grid.arrange(UAI_high_plot, MUTI_high_plot, UN_high_plot, ncol=3, nrow = 1)

#USE THE FOLLOWING in RESULTS SECTION 

#using patchwork package 
if(!require(patchwork)){
  install.packages("patchwork")
  require(patchwork)
}
library(patchwork)


FIG_All_nest_high_plots <- UAI_high_plot + MUTI_high_plot + UN_high_plot + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 14))

 print(FIG_All_nest_high_plots)

 saveRDS(FIG_All_nest_high_plots, here("Results", "FIG_all_nest_high"))

 
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


#time to plot!

library(effects)
# light - B8D5E9, dark - 6C9CCC
Safety.UNdb <- predictorEffect("nest.safety" , UN_M_nest_safety)


plot(Safety.UNdb, ask =FALSE, xlab = "Nest Safety", ylab = "UN", main ="", lines=list(multiline=TRUE, col=c("#B8D5E9")), confint=list(style="auto"), auto.key=FALSE, ylim=c(-20,10))

SAFETY.UN.DF<-data.frame(Safety.UNdb)

SAFETY_UN_plot<-ggplot(data=SAFETY.UN.DF,aes(x=nest.safety, y=fit)) +
  geom_line(color="#6C9CCC",lwd=1.5) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="#B8D5E9",alpha=.7, lwd=.1)+xlim(-1000,1000) +
  
  geom_point(data=NestTraitDat12,aes(x=jitter(nest.safety, 1), y =Urban),color="#6C9CCC", bg="#6C9CCC",alpha=.7, size=1,pch=21) +
  coord_cartesian(ylim = c(0, 1), xlim =c(0.9,4.1)) + theme_classic() +
  theme(axis.text.x =  element_text(color="black", size = 10),
        axis.text.y =  element_text(color="black", size = 10), 
        axis.title.x = element_text(size =14, margin = margin(t = 8)), 
        axis.title.y = element_text(size =14)) +
  xlab("Nest Safety") + 
  ylab("UN") 

SAFETY_UN_plot


#save, even though this must be fixed to plot NEW phyloglm accurately 

saveRDS(SAFETY_UN_plot, here("Results", "need_fix_FIG_nest_safety_UN.rds"))





####################################################################################
#okay, now let's arrange all nesting traits! 


UAI_high_plot <- readRDS(here("Results", "FIG_high_nest_UAI.rds"))
MUTI_high_plot <- readRDS(here("Results", "FIG_high_nest_MUTI.rds"))
UN_high_plot <- readRDS(here("Results", "FIG_high_nest_UN.rds"))


UAI_low_plot <- readRDS(here("Results", "FIG_low_nest_UAI.rds"))
MUTI_low_plot <- readRDS(here("Results", "FIG_low_nest_MUTI.rds"))


UN_safety_plot <- readRDS(here("Results", "need_fix_FIG_nest_safety_UN.rds"))


low_nesting_plots <- grid.arrange(UAI_low_plot, MUTI_low_plot, ncol=2, nrow = 1)



all_nesting_plots <- grid.arrange(UAI_low_plot, MUTI_low_plot, UAI_high_plot, MUTI_high_plot, UN_high_plot, UN_safety_plot, ncol=3, nrow = 2)

#specialty plot 

layout <- rbind(c(1, 2, NA),  # First row: UAI_low_plot, MUTI_low_plot
                c(3, 4, 5),   # Second row: UAI_high_plot, MUTI_high_plot, UN_high_plot
                c(6, NA, NA)) # Third row: UN_safety_plot

all_nesting_plots <- grid.arrange(UAI_low_plot, MUTI_low_plot, UAI_high_plot, 
                                  MUTI_high_plot, UN_high_plot, UN_safety_plot, 
                                  layout_matrix = layout)







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



#test 



