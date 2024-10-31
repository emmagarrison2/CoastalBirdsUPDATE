##########################################################
######### Coastal Birds Urban Tolerance Project ##########
##########################################################
# Combine model outputs into tables
# Authors: Sarah L. Jennings, Emma M. Garrison
##########################################################

# load packages
library(here)
library(tidyverse)
library(broom)
library(broom.mixed)
library(openxlsx)
library(nlme)
library(phylolm)

#################### Sensory models #######################

#### Load models ####
# C.T.
UAI_GLS_C.T <- readRDS(here("Models/UAI", "UAI_GLS_C.T.rds"))
MUTI_GLS_C.T <- readRDS(here("Models/MUTI", "MUTI_GLS_C.T.rds"))
phyglm_UN_CT <- readRDS(here("Models/UN", "phyglm_UN_CT_scale.rds"))

# Peak freq                                   
UAI_GLS_pf <- readRDS(here("Models/UAI", "UAI_GLS_pf.rds"))
MUTI_GLS_pf <- readRDS(here("Models/MUTI", "MUTI_GLS_pf.rds"))
phyglm_UN_pf <- readRDS(here("Models/UN", "phyglm_UN_pf_scale.rds"))

#### Create tidy model outputs #####
# C.T models
UAI_CT_tidy <- tidy(UAI_GLS_C.T, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "UAI",
         model_type = "gls",
         trait_group = "sensory",
         predictor_trait = "C.T ratio",
         sample_size = UAI_GLS_C.T$dims$N,
         lambda = 0, # if lambda was fixed in the model, need to manually enter it here
         alpha = as.numeric(NA)) %>%
  select(-p.value)

MUTI_CT_tidy <-tidy(MUTI_GLS_C.T, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "MUTI",
         model_type = "gls",
         trait_group = "sensory",
         predictor_trait = "C.T ratio",
         sample_size = MUTI_GLS_C.T$dims$N,
         lambda = 0.1, # if lambda was fixed in the model, need to manually enter it here
         alpha = as.numeric(NA)) %>%
  select(-p.value)

UN_CT_tidy <- as.data.frame(summary(phyglm_UN_CT)$coefficients) %>% 
  rownames_to_column(., var ="term") %>%
  select(-p.value) %>%
  rename(estimate = Estimate, std.error = StdErr, # rename columns to match UAI and MUTI outputs
         statistic = z.value, conf.low = lowerbootCI, 
         conf.high = upperbootCI) %>%
  mutate(index = "UN",
         model_type = "phyloglm",
         trait_group = "sensory",
         predictor_trait = "C.T ratio",
         sample_size = phyglm_UN_CT$n,
         lambda = as.numeric(NA),
         alpha = phyglm_UN_CT$alpha)

# peak frequency models
UAI_pf_tidy <- tidy(UAI_GLS_pf, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "UAI",
         model_type = "gls",
         trait_group = "sensory",
         predictor_trait = "peak frequency",
         sample_size = UAI_GLS_pf$dims$N,
         lambda = as.numeric(UAI_GLS_pf$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model)
         alpha = as.numeric(NA)) %>%
  select(-p.value)

MUTI_pf_tidy <-tidy(MUTI_GLS_pf, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "MUTI",
         model_type = "gls",
         trait_group = "sensory",
         predictor_trait = "peak frequency",
         sample_size = MUTI_GLS_pf$dims$N,
         lambda = as.numeric(MUTI_GLS_pf$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model), 
         alpha = as.numeric(NA)) %>%
  select(-p.value)

UN_pf_tidy <- as.data.frame(summary(phyglm_UN_pf)$coefficients) %>% 
  rownames_to_column(., var ="term") %>%
  select(-p.value) %>%
  rename(estimate = Estimate, std.error = StdErr, # rename columns to match UAI and MUTI outputs
         statistic = z.value, conf.low = lowerbootCI, 
         conf.high = upperbootCI) %>%
  mutate(index = "UN",
         model_type = "phyloglm",
         trait_group = "sensory",
         predictor_trait = "peak frequency",
         sample_size = phyglm_UN_pf$n,
         lambda = as.numeric(NA),
         alpha = phyglm_UN_pf$alpha)

#### Combine outputs #####
# combine, re-order columns, round to 4 decimal places
sensory_tidy <- bind_rows(UAI_CT_tidy, MUTI_CT_tidy, UN_CT_tidy,
                          UAI_pf_tidy, MUTI_pf_tidy, UN_pf_tidy) %>%
  mutate_if(is.numeric, round, 3) %>% # round numeric columns to 3 decimal places
  select(index:sample_size, term:alpha)

#################### Diet models #######################

#### Load models ####
# Diet invert
UAI_GLS_invert <- readRDS(here("Models/UAI", "UAI_GLS_invert.rds"))
MUTI_GLS_invert <- readRDS(here("Models/MUTI", "MUTI_GLS_invert.rds"))
phyglm_UN_invert <- readRDS(here("Models/UN", "phyglm_UN_Invert_fix.rds"))

# Diet vert                                 
UAI_GLS_vert <- readRDS(here("Models/UAI", "UAI_GLS_vert.rds"))
MUTI_GLS_vert <- readRDS(here("Models/MUTI", "MUTI_GLS_vert.rds"))
phyglm_UN_vert <- readRDS(here("Models/UN", "phyglm_UN_Vert_fix.rds"))

# Diet plant/seed
UAI_GLS_PS <- readRDS(here("Models/UAI", "UAI_GLS_PS.rds"))
MUTI_GLS_PS <- readRDS(here("Models/MUTI", "MUTI_GLS_PS.rds"))
phyglm_UN_PS <- readRDS(here("Models/UN", "phyglm_UN_PS_scale.rds"))

# Diet fruit/nut
UAI_GLS_FN <- readRDS(here("Models/UAI", "UAI_GLS_FN.rds"))
MUTI_GLS_FN <- readRDS(here("Models/MUTI", "MUTI_GLS_FN.rds"))
phyglm_UN_FN <- readRDS(here("Models/UN", "phyglm_UN_FN_scale.rds"))


#### Create tidy model outputs #####
# Invert models
UAI_invert_tidy <- tidy(UAI_GLS_invert, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "UAI",
         model_type = "gls",
         trait_group = "diet",
         predictor_trait = "% invertebrates",
         sample_size = UAI_GLS_invert$dims$N,
         lambda = as.numeric(UAI_GLS_invert$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model)
         alpha = as.numeric(NA)) %>%
  select(-p.value)

MUTI_invert_tidy <-tidy(MUTI_GLS_invert, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "MUTI",
         model_type = "gls",
         trait_group = "diet",
         predictor_trait = "% invertebrates",
         sample_size = MUTI_GLS_invert$dims$N,
         lambda = as.numeric(MUTI_GLS_invert$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model), 
         alpha = as.numeric(NA)) %>%
  select(-p.value)

UN_invert_tidy <- as.data.frame(summary(phyglm_UN_invert)$coefficients) %>% 
  rownames_to_column(., var ="term") %>%
  select(-p.value) %>%
  rename(estimate = Estimate, std.error = StdErr, # rename columns to match UAI and MUTI outputs
         statistic = z.value, conf.low = lowerbootCI, 
         conf.high = upperbootCI) %>%
  mutate(index = "UN",
         model_type = "phyloglm",
         trait_group = "diet",
         predictor_trait = "% invertebrates",
         sample_size = phyglm_UN_invert$n,
         lambda = as.numeric(NA),
         alpha = phyglm_UN_invert$alpha)

# Vertebrate models
UAI_vert_tidy <- tidy(UAI_GLS_vert, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "UAI",
         model_type = "gls",
         trait_group = "diet",
         predictor_trait = "% vertebrates",
         sample_size = UAI_GLS_vert$dims$N,
         lambda = as.numeric(UAI_GLS_vert$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model)
         alpha = as.numeric(NA)) %>%
  select(-p.value)

MUTI_vert_tidy <-tidy(MUTI_GLS_vert, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "MUTI",
         model_type = "gls",
         trait_group = "diet",
         predictor_trait = "% vertebrates",
         sample_size = MUTI_GLS_vert$dims$N,
         lambda = as.numeric(MUTI_GLS_vert$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model), 
         alpha = as.numeric(NA)) %>%
  select(-p.value)

UN_vert_tidy <- as.data.frame(summary(phyglm_UN_vert)$coefficients) %>% 
  rownames_to_column(., var ="term") %>%
  select(-p.value) %>%
  rename(estimate = Estimate, std.error = StdErr, # rename columns to match UAI and MUTI outputs
         statistic = z.value, conf.low = lowerbootCI, 
         conf.high = upperbootCI) %>%
  mutate(index = "UN",
         model_type = "phyloglm",
         trait_group = "diet",
         predictor_trait = "% vertebrates",
         sample_size = phyglm_UN_vert$n,
         lambda = as.numeric(NA),
         alpha = phyglm_UN_vert$alpha)

# Plant/seed models
UAI_PS_tidy <- tidy(UAI_GLS_PS, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "UAI",
         model_type = "gls",
         trait_group = "diet",
         predictor_trait = "% plant/seed",
         sample_size = UAI_GLS_PS$dims$N,
         lambda = as.numeric(UAI_GLS_PS$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model)
         alpha = as.numeric(NA)) %>%
  select(-p.value)

MUTI_PS_tidy <-tidy(MUTI_GLS_PS, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "MUTI",
         model_type = "gls",
         trait_group = "diet",
         predictor_trait = "% plant/seed",
         sample_size = MUTI_GLS_PS$dims$N,
         lambda = as.numeric(MUTI_GLS_PS$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model), 
         alpha = as.numeric(NA)) %>%
  select(-p.value)

UN_PS_tidy <- as.data.frame(summary(phyglm_UN_PS)$coefficients) %>% 
  rownames_to_column(., var ="term") %>%
  select(-p.value) %>%
  rename(estimate = Estimate, std.error = StdErr, # rename columns to match UAI and MUTI outputs
         statistic = z.value, conf.low = lowerbootCI, 
         conf.high = upperbootCI) %>%
  mutate(index = "UN",
         model_type = "phyloglm",
         trait_group = "diet",
         predictor_trait = "% plant/seed",
         sample_size = phyglm_UN_PS$n,
         lambda = as.numeric(NA),
         alpha = phyglm_UN_PS$alpha)

# Fruit/nectar models
UAI_FN_tidy <- tidy(UAI_GLS_FN, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "UAI",
         model_type = "gls",
         trait_group = "diet",
         predictor_trait = "% fruit/nectar",
         sample_size = UAI_GLS_FN$dims$N,
         lambda = as.numeric(UAI_GLS_FN$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model)
         alpha = as.numeric(NA)) %>%
  select(-p.value)

MUTI_FN_tidy <-tidy(MUTI_GLS_FN, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "MUTI",
         model_type = "gls",
         trait_group = "diet",
         predictor_trait = "% fruit/nectar",
         sample_size = MUTI_GLS_FN$dims$N,
         lambda = as.numeric(MUTI_GLS_FN$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model), 
         alpha = as.numeric(NA)) %>%
  select(-p.value)

UN_FN_tidy <- as.data.frame(summary(phyglm_UN_FN)$coefficients) %>% 
  rownames_to_column(., var ="term") %>%
  select(-p.value) %>%
  rename(estimate = Estimate, std.error = StdErr, # rename columns to match UAI and MUTI outputs
         statistic = z.value, conf.low = lowerbootCI, 
         conf.high = upperbootCI) %>%
  mutate(index = "UN",
         model_type = "phyloglm",
         trait_group = "diet",
         predictor_trait = "% fruit/nectar",
         sample_size = phyglm_UN_FN$n,
         lambda = as.numeric(NA),
         alpha = phyglm_UN_FN$alpha)

#### Combine outputs #####
# combine, re-order columns, round to 4 decimal places
diet_tidy <- bind_rows(UAI_invert_tidy, MUTI_invert_tidy, UN_invert_tidy,
                       UAI_vert_tidy, MUTI_vert_tidy, UN_vert_tidy,
                       UAI_PS_tidy, MUTI_PS_tidy, UN_PS_tidy,
                       UAI_FN_tidy, MUTI_FN_tidy, UN_FN_tidy) %>%
  mutate_if(is.numeric, round, 3) %>% # round numeric columns to 3 decimal places
  select(index:sample_size, term:alpha)

#################### Life history models #######################

#### Load models ####
# brood value
UAI_GLS_brood <- readRDS(here("Models/UAI", "UAI_GLS_bv.rds"))
MUTI_GLS_brood <- readRDS(here("Models/MUTI", "MUTI_GLS_bv.rds"))
phyglm_UN_brood <- readRDS(here("Models/UN", "phyglm_UN_bv_scale.rds"))

# clutch size                               
UAI_GLS_clutch <- readRDS(here("Models/UAI", "UAI_GLS_clutch.rds"))
MUTI_GLS_clutch <- readRDS(here("Models/MUTI", "MUTI_GLS_clutch.rds"))
phyglm_UN_clutch <- readRDS(here("Models/UN", "phyglm_UN_clutch_fix.rds"))

# longevity
UAI_GLS_long <- readRDS(here("Models/UAI", "UAI_GLS_long.rds"))
MUTI_GLS_long <- readRDS(here("Models/MUTI", "MUTI_GLS_long.rds"))
phyglm_UN_long <- readRDS(here("Models/UN", "phyglm_UN_long_fix.rds"))

# developmental mode
UAI_GLS_develop <- readRDS(here("Models/UAI", "UAI_GLS_develop.rds"))
MUTI_GLS_develop <- readRDS(here("Models/MUTI", "MUTI_GLS_develop.rds"))
phyglm_UN_develop <- readRDS(here("Models/UN", "phyglm_UN_develop_fix.rds"))


#### Create tidy model outputs #####

# Brood value models
UAI_brood_tidy <- tidy(UAI_GLS_brood, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "UAI",
         model_type = "gls",
         trait_group = "life history",
         predictor_trait = "brood value",
         sample_size = UAI_GLS_brood$dims$N,
         lambda = as.numeric(UAI_GLS_brood$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model)
         alpha = as.numeric(NA)) %>%
  select(-p.value)

MUTI_brood_tidy <-tidy(MUTI_GLS_brood, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "MUTI",
         model_type = "gls",
         trait_group = "life history",
         predictor_trait = "brood value",
         sample_size = MUTI_GLS_brood$dims$N,
         lambda = as.numeric(MUTI_GLS_brood$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model), 
         alpha = as.numeric(NA)) %>%
  select(-p.value)

UN_brood_tidy <- as.data.frame(summary(phyglm_UN_brood)$coefficients) %>% 
  rownames_to_column(., var ="term") %>%
  select(-p.value) %>%
  rename(estimate = Estimate, std.error = StdErr, # rename columns to match UAI and MUTI outputs
         statistic = z.value, conf.low = lowerbootCI, 
         conf.high = upperbootCI) %>%
  mutate(index = "UN",
         model_type = "phyloglm",
         trait_group = "life history",
         predictor_trait = "brood value",
         sample_size = phyglm_UN_brood$n,
         lambda = as.numeric(NA),
         alpha = phyglm_UN_brood$alpha)

# Clutch size models
UAI_clutch_tidy <- tidy(UAI_GLS_clutch, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "UAI",
         model_type = "gls",
         trait_group = "life history",
         predictor_trait = "clutch size",
         sample_size = UAI_GLS_clutch$dims$N,
         lambda = as.numeric(UAI_GLS_clutch$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model)
         alpha = as.numeric(NA)) %>%
  select(-p.value)

MUTI_clutch_tidy <-tidy(MUTI_GLS_clutch, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "MUTI",
         model_type = "gls",
         trait_group = "life history",
         predictor_trait = "clutch size",
         sample_size = MUTI_GLS_clutch$dims$N,
         lambda = as.numeric(MUTI_GLS_clutch$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model), 
         alpha = as.numeric(NA)) %>%
  select(-p.value)

UN_clutch_tidy <- as.data.frame(summary(phyglm_UN_clutch)$coefficients) %>% 
  rownames_to_column(., var ="term") %>%
  select(-p.value) %>%
  rename(estimate = Estimate, std.error = StdErr, # rename columns to match UAI and MUTI outputs
         statistic = z.value, conf.low = lowerbootCI, 
         conf.high = upperbootCI) %>%
  mutate(index = "UN",
         model_type = "phyloglm",
         trait_group = "life history",
         predictor_trait = "clutch size",
         sample_size = phyglm_UN_clutch$n,
         lambda = as.numeric(NA),
         alpha = phyglm_UN_clutch$alpha)

# Longevity models
UAI_long_tidy <- tidy(UAI_GLS_long, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "UAI",
         model_type = "gls",
         trait_group = "life history",
         predictor_trait = "longevity",
         sample_size = UAI_GLS_long$dims$N,
         lambda = as.numeric(UAI_GLS_long$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model)
         alpha = as.numeric(NA)) %>%
  select(-p.value)

MUTI_long_tidy <-tidy(MUTI_GLS_long, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "MUTI",
         model_type = "gls",
         trait_group = "life history",
         predictor_trait = "longevity",
         sample_size = MUTI_GLS_long$dims$N,
         lambda = as.numeric(MUTI_GLS_long$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model), 
         alpha = as.numeric(NA)) %>%
  select(-p.value)

UN_long_tidy <- as.data.frame(summary(phyglm_UN_long)$coefficients) %>% 
  rownames_to_column(., var ="term") %>%
  select(-p.value) %>%
  rename(estimate = Estimate, std.error = StdErr, # rename columns to match UAI and MUTI outputs
         statistic = z.value, conf.low = lowerbootCI, 
         conf.high = upperbootCI) %>%
  mutate(index = "UN",
         model_type = "phyloglm",
         trait_group = "life history",
         predictor_trait = "longevity",
         sample_size = phyglm_UN_long$n,
         lambda = as.numeric(NA),
         alpha = phyglm_UN_long$alpha)

# Developmental mode models
UAI_develop_tidy <- tidy(UAI_GLS_develop, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "UAI",
         model_type = "gls",
         trait_group = "life history",
         predictor_trait = "developmental mode",
         sample_size = UAI_GLS_develop$dims$N,
         lambda = as.numeric(UAI_GLS_develop$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model)
         alpha = as.numeric(NA)) %>%
  select(-p.value)

MUTI_develop_tidy <-tidy(MUTI_GLS_develop, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "MUTI",
         model_type = "gls",
         trait_group = "life history",
         predictor_trait = "developmental mode",
         sample_size = MUTI_GLS_develop$dims$N,
         lambda = as.numeric(MUTI_GLS_develop$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model), 
         alpha = as.numeric(NA)) %>%
  select(-p.value)

UN_develop_tidy <- as.data.frame(summary(phyglm_UN_develop)$coefficients) %>% 
  rownames_to_column(., var ="term") %>%
  select(-p.value) %>%
  rename(estimate = Estimate, std.error = StdErr, # rename columns to match UAI and MUTI outputs
         statistic = z.value, conf.low = lowerbootCI, 
         conf.high = upperbootCI) %>%
  mutate(index = "UN",
         model_type = "phyloglm",
         trait_group = "life history",
         predictor_trait = "developmental mode",
         sample_size = phyglm_UN_develop$n,
         lambda = as.numeric(NA),
         alpha = phyglm_UN_develop$alpha)

#### Combine outputs #####
# combine, re-order columns, round to 4 decimal places
lifehistory_tidy <- bind_rows(UAI_brood_tidy, MUTI_brood_tidy, UN_brood_tidy,
                              UAI_clutch_tidy, MUTI_clutch_tidy, UN_clutch_tidy,
                              UAI_long_tidy, MUTI_long_tidy, UN_long_tidy,
                              UAI_develop_tidy, MUTI_develop_tidy, UN_develop_tidy) %>%
  mutate_if(is.numeric, round, 3) %>% # round numeric columns to 3 decimal places
  select(index:sample_size, term:alpha)

#################### Sexual selection models #######################

#### Load models ####
# brightness
UAI_GLS_bright <- readRDS(here("Models/UAI", "UAI_GLS_bright.rds"))
MUTI_GLS_bright <- readRDS(here("Models/MUTI", "MUTI_GLS_bright.rds"))
phyglm_UN_bright <- readRDS(here("Models/UN", "phyglm_UN_brightness_scale.rds"))

# hue                              
UAI_GLS_hue <- readRDS(here("Models/UAI", "UAI_GLS_hue.rds"))
MUTI_GLS_hue <- readRDS(here("Models/MUTI", "MUTI_GLS_hue.rds"))
phyglm_UN_hue <- readRDS(here("Models/UN", "phyglm_UN_hue_scale.rds"))

# SS intensity males
UAI_GLS_ssm <- readRDS(here("Models/UAI", "UAI_GLS_ssm.rds"))
MUTI_GLS_ssm <- readRDS(here("Models/MUTI", "MUTI_GLS_ssm.rds"))
phyglm_UN_ssm <- readRDS(here("Models/UN", "phyglm_UN_ssm_fix.rds"))

# SS intensity females
UAI_GLS_ssf <- readRDS(here("Models/UAI", "UAI_GLS_ssf.rds"))
MUTI_GLS_ssf <- readRDS(here("Models/MUTI", "MUTI_GLS_ssf.rds"))
phyglm_UN_ssf <- readRDS(here("Models/UN", "phyglm_UN_ssf_fix.rds"))


#### Create tidy model outputs #####

# brightness models
UAI_bright_tidy <- tidy(UAI_GLS_bright, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "UAI",
         model_type = "gls",
         trait_group = "sexual selection",
         predictor_trait = "dichromatism: brightness",
         sample_size = UAI_GLS_bright$dims$N,
         lambda = as.numeric(UAI_GLS_bright$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model)
         alpha = as.numeric(NA)) %>%
  select(-p.value)

MUTI_bright_tidy <-tidy(MUTI_GLS_bright, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "MUTI",
         model_type = "gls",
         trait_group = "sexual selection",
         predictor_trait = "dichromatism: brightness",
         sample_size = MUTI_GLS_bright$dims$N,
         lambda = 0.4, # manually input lambda as it was fixed in the model, 
         alpha = as.numeric(NA)) %>%
  select(-p.value)

UN_bright_tidy <- as.data.frame(summary(phyglm_UN_bright)$coefficients) %>% 
  rownames_to_column(., var ="term") %>%
  select(-p.value) %>%
  rename(estimate = Estimate, std.error = StdErr, # rename columns to match UAI and MUTI outputs
         statistic = z.value, conf.low = lowerbootCI, 
         conf.high = upperbootCI) %>%
  mutate(index = "UN",
         model_type = "phyloglm",
         trait_group = "sexual selection",
         predictor_trait = "dichromatism: brightness",
         sample_size = phyglm_UN_bright$n,
         lambda = as.numeric(NA),
         alpha = phyglm_UN_bright$alpha)

# hue models
UAI_hue_tidy <- tidy(UAI_GLS_hue, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "UAI",
         model_type = "gls",
         trait_group = "sexual selection",
         predictor_trait = "dichromatism: hue",
         sample_size = UAI_GLS_hue$dims$N,
         lambda = as.numeric(UAI_GLS_hue$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model)
         alpha = as.numeric(NA)) %>%
  select(-p.value)

MUTI_hue_tidy <-tidy(MUTI_GLS_hue, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "MUTI",
         model_type = "gls",
         trait_group = "sexual selection",
         predictor_trait = "dichromatism: hue",
         sample_size = MUTI_GLS_hue$dims$N,
         lambda = as.numeric(MUTI_GLS_hue$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model), 
         alpha = as.numeric(NA)) %>%
  select(-p.value)

UN_hue_tidy <- as.data.frame(summary(phyglm_UN_hue)$coefficients) %>% 
  rownames_to_column(., var ="term") %>%
  select(-p.value) %>%
  rename(estimate = Estimate, std.error = StdErr, # rename columns to match UAI and MUTI outputs
         statistic = z.value, conf.low = lowerbootCI, 
         conf.high = upperbootCI) %>%
  mutate(index = "UN",
         model_type = "phyloglm",
         trait_group = "sexual selection",
         predictor_trait = "dichromatism: hue",
         sample_size = phyglm_UN_hue$n,
         lambda = as.numeric(NA),
         alpha = phyglm_UN_hue$alpha)

# ssm models
UAI_ssm_tidy <- tidy(UAI_GLS_ssm, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "UAI",
         model_type = "gls",
         trait_group = "sexual selection",
         predictor_trait = "sexual selection intensity: males",
         sample_size = UAI_GLS_ssm$dims$N,
         lambda = as.numeric(UAI_GLS_ssm$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model)
         alpha = as.numeric(NA)) %>%
  select(-p.value)

MUTI_ssm_tidy <-tidy(MUTI_GLS_ssm, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "MUTI",
         model_type = "gls",
         trait_group = "sexual selection",
         predictor_trait = "sexual selection intensity: males",
         sample_size = MUTI_GLS_ssm$dims$N,
         lambda = as.numeric(MUTI_GLS_ssm$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model), 
         alpha = as.numeric(NA)) %>%
  select(-p.value)

UN_ssm_tidy <- as.data.frame(summary(phyglm_UN_ssm)$coefficients) %>% 
  rownames_to_column(., var ="term") %>%
  select(-p.value) %>%
  rename(estimate = Estimate, std.error = StdErr, # rename columns to match UAI and MUTI outputs
         statistic = z.value, conf.low = lowerbootCI, 
         conf.high = upperbootCI) %>%
  mutate(index = "UN",
         model_type = "phyloglm",
         trait_group = "sexual selection",
         predictor_trait = "sexual selection intensity: males",
         sample_size = phyglm_UN_ssm$n,
         lambda = as.numeric(NA),
         alpha = phyglm_UN_ssm$alpha)

# ssf mode models
UAI_ssf_tidy <- tidy(UAI_GLS_ssf, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "UAI",
         model_type = "gls",
         trait_group = "sexual selection",
         predictor_trait = "sexual selection intensity: females",
         sample_size = UAI_GLS_ssf$dims$N,
         lambda = as.numeric(UAI_GLS_ssf$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model)
         alpha = as.numeric(NA)) %>%
  select(-p.value)

MUTI_ssf_tidy <-tidy(MUTI_GLS_ssf, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "MUTI",
         model_type = "gls",
         trait_group = "sexual selection",
         predictor_trait = "sexual selection intensity: females",
         sample_size = MUTI_GLS_ssf$dims$N,
         lambda = as.numeric(MUTI_GLS_ssf$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model), 
         alpha = as.numeric(NA)) %>%
  select(-p.value)

UN_ssf_tidy <- as.data.frame(summary(phyglm_UN_ssf)$coefficients) %>% 
  rownames_to_column(., var ="term") %>%
  select(-p.value) %>%
  rename(estimate = Estimate, std.error = StdErr, # rename columns to match UAI and MUTI outputs
         statistic = z.value, conf.low = lowerbootCI, 
         conf.high = upperbootCI) %>%
  mutate(index = "UN",
         model_type = "phyloglm",
         trait_group = "sexual selection",
         predictor_trait = "sexual selection intensity: females",
         sample_size = phyglm_UN_ssf$n,
         lambda = as.numeric(NA),
         alpha = phyglm_UN_ssf$alpha)

#### Combine outputs #####
# combine, re-order columns, round to 4 decimal places
sexualselection_tidy <- bind_rows(UAI_bright_tidy, MUTI_bright_tidy, UN_bright_tidy,
                                  UAI_hue_tidy, MUTI_hue_tidy, UN_hue_tidy,
                                  UAI_ssm_tidy, MUTI_ssm_tidy, UN_ssm_tidy,
                                  UAI_ssf_tidy, MUTI_ssf_tidy, UN_ssf_tidy) %>%
  mutate_if(is.numeric, round, 3) %>% # round numeric columns to 3 decimal places
  select(index:sample_size, term:alpha)

#################### Social models #######################

#### Load models ####
# territoriality
UAI_GLS_territory <- readRDS(here("Models/UAI", "UAI_GLS_territory.rds"))
MUTI_GLS_territory <- readRDS(here("Models/MUTI", "MUTI_GLS_territory.rds"))
phyglm_UN_territory <- readRDS(here("Models/UN", "phyglm_UN_territorial_scale.rds"))

# cooperative breeding - only UAI                             
UAI_GLS_cooperative <- readRDS(here("Models/UAI", "UAI_GLS_cooperative.rds"))

#### Create tidy model outputs #####

# territoriality models
UAI_territory_tidy <- tidy(UAI_GLS_territory, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "UAI",
         model_type = "gls",
         trait_group = "social",
         predictor_trait = " territoriality",
         sample_size = UAI_GLS_territory$dims$N,
         lambda = as.numeric(UAI_GLS_territory$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model)
         alpha = as.numeric(NA)) %>%
  select(-p.value)

MUTI_territory_tidy <-tidy(MUTI_GLS_territory, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "MUTI",
         model_type = "gls",
         trait_group = "social",
         predictor_trait = " territoriality",
         sample_size = MUTI_GLS_territory$dims$N,
         lambda = 0.4, # manually input lambda as it was fixed in the model, 
         alpha = as.numeric(NA)) %>%
  select(-p.value)

UN_territory_tidy <- as.data.frame(summary(phyglm_UN_territory)$coefficients) %>% 
  rownames_to_column(., var ="term") %>%
  select(-p.value) %>%
  rename(estimate = Estimate, std.error = StdErr, # rename columns to match UAI and MUTI outputs
         statistic = z.value, conf.low = lowerbootCI, 
         conf.high = upperbootCI) %>%
  mutate(index = "UN",
         model_type = "phyloglm",
         trait_group = "social",
         predictor_trait = " territoriality",
         sample_size = phyglm_UN_territory$n,
         lambda = as.numeric(NA),
         alpha = phyglm_UN_territory$alpha)

# cooperative breeding models
UAI_cooperative_tidy <- tidy(UAI_GLS_cooperative, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "UAI",
         model_type = "gls",
         trait_group = "social",
         predictor_trait = " cooperative breeding",
         sample_size = UAI_GLS_cooperative$dims$N,
         lambda = as.numeric(UAI_GLS_cooperative$modelStruct$corStruct), # extract lambda from model when not-fixed (selected by model)
         alpha = as.numeric(NA)) %>%
  select(-p.value)

#### Combine outputs #####
# combine, re-order columns, round to 4 decimal places
social_tidy <- bind_rows(UAI_territory_tidy, MUTI_territory_tidy, UN_territory_tidy,
                         UAI_cooperative_tidy) %>%
  mutate_if(is.numeric, round, 3) %>% # round numeric columns to 3 decimal places
  select(index:sample_size, term:alpha)

################## Nest Models ############################

#### Load models ####
# Site: low
UAI_GLS_low <- readRDS(here("Models/UAI", "UAI_GLS_nest_low.rds"))
MUTI_GLS_low <- readRDS(here("Models/MUTI", "MUTI_GLS_nest_low.rds"))
phyglm_UN_low <- readRDS(here("Models/UN", "phyglm_UN_nest_low_fix.rds"))

# Site: low only                            
UAI_GLS_low_only <- readRDS(here("Models/UAI", "UAI_GLS_nest_low_only.rds"))
MUTI_GLS_low_only <- readRDS(here("Models/MUTI", "MUTI_GLS_nest_low_only.rds"))
phyglm_UN_low_only <- readRDS(here("Models/UN", "phyglm_UN_nest_low_only_scale.rds"))

# Site: high
UAI_GLS_high <- readRDS(here("Models/UAI", "UAI_GLS_nest_high.rds"))
MUTI_GLS_high <- readRDS(here("Models/MUTI", "MUTI_GLS_nest_high.rds"))
phyglm_UN_high <- readRDS(here("Models/UN", "phyglm_UN_nest_high_fix.rds"))

# Strategy (open/enclosed) - not run for UN
UAI_GLS_neststr <- readRDS(here("Models/UAI", "UAI_GLS_neststr.rds"))
MUTI_GLS_neststr <- readRDS(here("Models/MUTI", "MUTI_GLS_neststr.rds"))

# Nest safety
UAI_GLS_nest_safety <- readRDS(here("Models/UAI", "UAI_GLS_nest_safety.rds"))
MUTI_GLS_nest_safety <- readRDS(here("Models/MUTI", "MUTI_GLS_nest_safety.rds"))
phyglm_UN_nest_safety <- readRDS(here("Models/UN", "phyglm_UN_nest_safety_scale.rds"))

#### Create tidy model outputs #####

# nest site low models
UAI_low_tidy <- tidy(UAI_GLS_low, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "UAI",
         model_type = "gls",
         trait_group = "nest",
         predictor_trait = "nest site low",
         sample_size = UAI_GLS_low$dims$N,
         lambda = as.numeric(UAI_GLS_low$modelStruct$corStruct), 
         alpha = as.numeric(NA)) %>%
  select(-p.value)

MUTI_low_tidy <-tidy(MUTI_GLS_low, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "MUTI",
         model_type = "gls",
         trait_group = "nest",
         predictor_trait = "nest site low",
         sample_size = MUTI_GLS_low$dims$N,
         lambda = as.numeric(MUTI_GLS_low$modelStruct$corStruct), 
         alpha = as.numeric(NA)) %>%
  select(-p.value)

UN_low_tidy <- as.data.frame(summary(phyglm_UN_low)$coefficients) %>% 
  rownames_to_column(., var ="term") %>%
  select(-p.value) %>%
  rename(estimate = Estimate, std.error = StdErr, # rename columns to match UAI and MUTI outputs
         statistic = z.value, conf.low = lowerbootCI, 
         conf.high = upperbootCI) %>%
  mutate(index = "UN",
         model_type = "phyloglm",
         trait_group = "nest",
         predictor_trait = "nest site low",
         sample_size = phyglm_UN_low$n,
         lambda = as.numeric(NA),
         alpha = phyglm_UN_low$alpha)

# nest site low ONLY models
UAI_low_only_tidy <- tidy(UAI_GLS_low_only, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "UAI",
         model_type = "gls",
         trait_group = "nest",
         predictor_trait = "nest site low (only)",
         sample_size = UAI_GLS_low_only$dims$N,
         lambda = as.numeric(UAI_GLS_low_only$modelStruct$corStruct), 
         alpha = as.numeric(NA)) %>%
  select(-p.value)

MUTI_low_only_tidy <-tidy(MUTI_GLS_low_only, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "MUTI",
         model_type = "gls",
         trait_group = "nest",
         predictor_trait = "nest site low (only)",
         sample_size = MUTI_GLS_low_only$dims$N,
         lambda = 0.3, 
         alpha = as.numeric(NA)) %>%
  select(-p.value)

UN_low_only_tidy <- as.data.frame(summary(phyglm_UN_low_only)$coefficients) %>% 
  rownames_to_column(., var ="term") %>%
  select(-p.value) %>%
  rename(estimate = Estimate, std.error = StdErr, # rename columns to match UAI and MUTI outputs
         statistic = z.value, conf.low = lowerbootCI, 
         conf.high = upperbootCI) %>%
  mutate(index = "UN",
         model_type = "phyloglm",
         trait_group = "nest",
         predictor_trait = "nest site low (only)",
         sample_size = phyglm_UN_low_only$n,
         lambda = as.numeric(NA),
         alpha = phyglm_UN_low_only$alpha)

# nest site high models
UAI_high_tidy <- tidy(UAI_GLS_high, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "UAI",
         model_type = "gls",
         trait_group = "nest",
         predictor_trait = "nest site high",
         sample_size = UAI_GLS_high$dims$N,
         lambda = as.numeric(UAI_GLS_high$modelStruct$corStruct), 
         alpha = as.numeric(NA)) %>%
  select(-p.value)

MUTI_high_tidy <-tidy(MUTI_GLS_high, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "MUTI",
         model_type = "gls",
         trait_group = "nest",
         predictor_trait = "nest site high",
         sample_size = MUTI_GLS_high$dims$N,
         lambda = as.numeric(MUTI_GLS_high$modelStruct$corStruct), 
         alpha = as.numeric(NA)) %>%
  select(-p.value)

UN_high_tidy <- as.data.frame(summary(phyglm_UN_high)$coefficients) %>% 
  rownames_to_column(., var ="term") %>%
  select(-p.value) %>%
  rename(estimate = Estimate, std.error = StdErr, # rename columns to match UAI and MUTI outputs
         statistic = z.value, conf.low = lowerbootCI, 
         conf.high = upperbootCI) %>%
  mutate(index = "UN",
         model_type = "phyloglm",
         trait_group = "nest",
         predictor_trait = "nest site high",
         sample_size = phyglm_UN_high$n,
         lambda = as.numeric(NA),
         alpha = phyglm_UN_high$alpha)

# nest safety models
UAI_nest_safety_tidy <- tidy(UAI_GLS_nest_safety, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "UAI",
         model_type = "gls",
         trait_group = "nest",
         predictor_trait = "nest safety",
         sample_size = UAI_GLS_nest_safety$dims$N,
         lambda = as.numeric(UAI_GLS_nest_safety$modelStruct$corStruct), 
         alpha = as.numeric(NA)) %>%
  select(-p.value)

MUTI_nest_safety_tidy <-tidy(MUTI_GLS_nest_safety, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "MUTI",
         model_type = "gls",
         trait_group = "nest",
         predictor_trait = "nest safety",
         sample_size = MUTI_GLS_nest_safety$dims$N,
         lambda = as.numeric(MUTI_GLS_nest_safety$modelStruct$corStruct), 
         alpha = as.numeric(NA)) %>%
  select(-p.value)

UN_nest_safety_tidy <- as.data.frame(summary(phyglm_UN_nest_safety)$coefficients) %>% 
  rownames_to_column(., var ="term") %>%
  select(-p.value) %>%
  rename(estimate = Estimate, std.error = StdErr, # rename columns to match UAI and MUTI outputs
         statistic = z.value, conf.low = lowerbootCI, 
         conf.high = upperbootCI) %>%
  mutate(index = "UN",
         model_type = "phyloglm",
         trait_group = "nest",
         predictor_trait = "nest safety",
         sample_size = phyglm_UN_nest_safety$n,
         lambda = as.numeric(NA),
         alpha = phyglm_UN_nest_safety$alpha)

# nest structure (open/enclosed) models
UAI_neststr_tidy <- tidy(UAI_GLS_neststr, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "UAI",
         model_type = "gls",
         trait_group = "nest",
         predictor_trait = "nest structure (open/enclosed)",
         sample_size = UAI_GLS_neststr$dims$N,
         lambda = as.numeric(UAI_GLS_neststr$modelStruct$corStruct), 
         alpha = as.numeric(NA)) %>%
  select(-p.value)

MUTI_neststr_tidy <-tidy(MUTI_GLS_neststr, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "MUTI",
         model_type = "gls",
         trait_group = "nest",
         predictor_trait = "nest structure (open/enclosed)",
         sample_size = MUTI_GLS_neststr$dims$N,
         lambda = as.numeric(MUTI_GLS_neststr$modelStruct$corStruct), 
         alpha = as.numeric(NA)) %>%
  select(-p.value)

#### Combine outputs #####
# combine, re-order columns, round to 4 decimal places
nest_tidy <- bind_rows(UAI_low_tidy, MUTI_low_tidy, UN_low_tidy,
                       UAI_low_only_tidy, MUTI_low_only_tidy, UN_low_only_tidy,
                       UAI_high_tidy, MUTI_high_tidy, UN_high_tidy,
                       UAI_nest_safety_tidy, MUTI_nest_safety_tidy, UN_nest_safety_tidy,
                       UAI_neststr_tidy, MUTI_neststr_tidy) %>%
  mutate_if(is.numeric, round, 3) %>% # round numeric columns to 3 decimal places
  select(index:sample_size, term:alpha)

################## Body Mass Models ############################

#### Load models ####
# body mass
UAI_GLS_mass <- readRDS(here("Models/UAI", "UAI_GLS_mass.rds"))
MUTI_GLS_mass <- readRDS(here("Models/MUTI", "MUTI_GLS_mass.rds"))
phyglm_UN_mass <- readRDS(here("Models/UN", "phyglm_UN_Mass_fix.rds"))

#### Create tidy model outputs #####

# body mass models
UAI_mass_tidy <- tidy(UAI_GLS_mass, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "UAI",
         model_type = "gls",
         trait_group = "body mass",
         predictor_trait = "log(body mass)",
         sample_size = UAI_GLS_mass$dims$N,
         lambda = as.numeric(UAI_GLS_mass$modelStruct$corStruct), 
         alpha = as.numeric(NA)) %>%
  select(-p.value)

MUTI_mass_tidy <-tidy(MUTI_GLS_mass, conf.int = T, conf.level = 0.95) %>%
  mutate(index = "MUTI",
         model_type = "gls",
         trait_group = "body mass",
         predictor_trait = "log(body mass)",
         sample_size = MUTI_GLS_mass$dims$N,
         lambda = as.numeric(MUTI_GLS_mass$modelStruct$corStruct), 
         alpha = as.numeric(NA)) %>%
  select(-p.value)

UN_mass_tidy <- as.data.frame(summary(phyglm_UN_mass)$coefficients) %>% 
  rownames_to_column(., var ="term") %>%
  select(-p.value) %>%
  rename(estimate = Estimate, std.error = StdErr, # rename columns to match UAI and MUTI outputs
         statistic = z.value, conf.low = lowerbootCI, 
         conf.high = upperbootCI) %>%
  mutate(index = "UN",
         model_type = "phyloglm",
         trait_group = "body mass",
         predictor_trait = "log(body mass)",
         sample_size = phyglm_UN_mass$n,
         lambda = as.numeric(NA),
         alpha = phyglm_UN_mass$alpha)

#### Combine outputs #####
# combine, re-order columns, round to 4 decimal places
bodymass_tidy <- bind_rows(UAI_mass_tidy, MUTI_mass_tidy, UN_mass_tidy) %>%
  mutate_if(is.numeric, round, 3) %>% # round numeric columns to 3 decimal places
  select(index:sample_size, term:alpha)

###########################################################################
###########################################################################
######## Put all results into a single spreadsheet and export #############

# make each data frame of model outputs into a separate sheet in the excel file
names <- list(
  "Body Mass" = bodymass_tidy, 
  "Sensory" = sensory_tidy,
  "Diet" = diet_tidy,
  "Life History" = lifehistory_tidy,
  "Social" = social_tidy,
  "Sexual Selection" = sexualselection_tidy,
  "Nest" = nest_tidy)

write.xlsx(names, file = here("Results", "RESULTS_TraitModelTables.xlsx"))
