### Objective: find list of birds that have ALL THREE UT indexes 

# load packages
library(ape)
library(geiger)
library(nlme)
library(effects)
library(ggplot2)
library(ggeffects)
library(sjPlot)
library(dplyr)
library(here)
library(tidyverse)

coastal_birds <- readRDS(here("Outputs", "Coastal_Species_w_Mass.rds"))
View(coastal_birds)
nrow(coastal_birds)
# 807 species 

# filter the dataset so that we only keep species with ALL THREE urban tolerance
# scores (UAI, MUTI, and UN)
coastal_birds_filtered <- coastal_birds %>%
  filter(!is.na(aveUAI) & !is.na(MUTIscore) & !is.na(Urban))
nrow(coastal_birds_filtered)
# there are 49 species with all three scores 

# save this list 

saveRDS(coastal_birds_filtered, here("Notes", "birds_w_all_three_scores.rds"))


