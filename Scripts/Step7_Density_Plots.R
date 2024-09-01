##### The objective of this script is to make a density plot to show the distribution 
# of Urban Tolerance scores across each of the three indexes (MUTI, UAI, and UN)

#Install/load required packages 


if(!require(ggplot2)){
  install.packages("ggplot2")
  require(ggplot2)
}
library(ggplot2)

if(!require(tidyquant)){
  install.packages("tidyquant")
  require(tidyquant)
}
library(tidyquant)

if(!require(tidyverse)){
  install.packages("tidyverse")
  require(tidyverse)
}
library(tidyverse)

if(!require(prediction)){
  install.packages("prediction")
  require(prediction)
}
library(prediction)


if(!require(here)){
  install.packages("here")
  require(here)
}
library(here)


if(!require(ggdist)){
  install.packages("ggdist")
  require(ggdist)
}
library(ggdist)



if(!require(ggthemes)){
  install.packages("ggthemes")
  require(ggthemes)
}
library(ggthemes)



if(!require(gridExtra)){
  install.packages("gridExtra")
  require(gridExtra)
}
library(gridExtra)



#pull in refined dataset of coastal species and index scores

AllIndexesCoastal <- readRDS(here("Data", "Coastal_Birds_List.rds"))
str(AllIndexesCoastal)
AllIndexesCoastal$Urban <- ifelse(AllIndexesCoastal$Urban == "U", 1, 0)
colnames(AllIndexesCoastal)

summary(AllIndexesCoastal$aveUAI)
summary(AllIndexesCoastal$Urban)
summary(AllIndexesCoastal$MUTIscore)

# let's see how many species have ALL three index scores! 

species_with_all_indexes <- AllIndexesCoastal %>% 
  filter(!is.na(aveUAI)) %>% 
  filter(!is.na(Urban)) %>% 
  filter(!is.na(MUTIscore)) 

nrow(species_with_all_indexes)
#49 species that are representatives of all three indexes! 
unique(species_with_all_indexes$Species_Jetz)

#######################################################################
#####density plot for aveUAI 

UAICoastal <- AllIndexesCoastal %>% 
  filter(!is.na(aveUAI) & is.finite(aveUAI))
nrow(UAICoastal)#798

UAI_density <- UAICoastal %>%
  ggplot(aes(x=aveUAI, fill=aveUAI)) + 
  stat_halfeye(
    adjust=0.5,
    justification=-0.2,
    .width= 0,
    point_colour = NA, 
    fill ='#FFC107', 
    alpha = 0.8
  ) + 
  geom_boxplot(
    width = 0.12, 
    color = 'black',
    fill= '#FFC107',
    outlier.colour = '#FFC107',
    alpha = 0.5
  ) +  
  labs (
    #title = "Coastal Score Distribution", 
    x = "Average UAI", 
    y = ""
  ) + 
  theme_classic() + 
  theme(
    #plot.title = element_text(hjust = 0.5), # Center the plot title
    axis.title = element_text(size = 12)    # Adjust axis title size
    # axis.text.y=element_blank() # Remove y-axis text labels
    # axis.ticks.y=element_blank()  # Remove y-axis ticks
  )


print(UAI_density)




#Raincloud density plot for MUTI 

MUTICoastal <- AllIndexesCoastal %>% 
  filter(!is.na(MUTIscore) & is.finite(MUTIscore))
nrow(MUTICoastal)#130

MUTI_density <- MUTICoastal %>%
  ggplot(aes(x=MUTIscore, fill=MUTIscore)) + 
  stat_halfeye(
    adjust=0.5,
    justification= -0.2,
    .width= 0,
    point_colour = NA, 
    fill ='#004D40', 
    alpha = 0.8
  ) + 
  geom_boxplot(
    width = 0.12, 
    color = 'black',
    fill= '#004D40',
    outlier.colour = '#004D40',
    alpha = 0.5
  ) +  
  labs (
    #title = "Coastal Score Distribution", 
    x = "Multivariate Urban Tolerance Index", 
    y = ""
  ) + 
  theme_classic() + 
  theme(
    #plot.title = element_text(hjust = 0.5), # Center the plot title
    axis.title = element_text(size = 12)   
    #axis.text.y=element_blank(), # Remove y-axis text labels
    #axis.ticks.y=element_blank()  # Remove y-axis ticks
  )    # Adjust axis title size


print(MUTI_density)



#histogram to display density of UN values (Urban or Non-Urban)  

UNCoastal <- AllIndexesCoastal %>% 
  filter(!is.na(Urban) & is.finite(Urban))
nrow(UNCoastal)

UN_density <- UNCoastal %>%
  ggplot(aes(x=Urban, fill=Urban)) + 
  geom_histogram(
    binwidth = 0.25,
    color = 'black', 
    fill ='#1E88E5', 
    alpha = 0.5
  )  +
  scale_y_continuous(
    breaks = c(0, 20, 40, 60, 80, 100),
    labels = c("0", "20", "40", "60", "80", "100")
  ) + 
  scale_x_continuous(
    breaks = c(0, 1),              # Set x-axis breaks
    labels = c("Non Urban", "Urban")           # Set x-axis labels 
  ) + 
  labs (
    #title = "Coastal Score Distribution", 
    x = "Urban Classification", 
    y = ""
  ) + 
  theme_classic() + 
  theme(
    #plot.title = element_text(hjust = 0.5), # Center the plot title
    axis.title = element_text(size = 12),   # Adjust axis title size
    axis.title.y = element_text(margin = margin(r = 10)) 
    #axis.text.y=element_blank(), # Remove y-axis text labels
    #axis.ticks.y=element_blank()  # Remove y-axis ticks
  )

print(UN_density)



#arrange them all side by side 
all_density <- grid.arrange(UAI_density, MUTI_density, UN_density, ncol=3, top= "Coastal Score Distributions", left = "Number of Species")


#?ggplot
