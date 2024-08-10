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

