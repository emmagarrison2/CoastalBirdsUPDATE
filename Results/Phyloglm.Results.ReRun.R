############### Objective of this script is to import Sarah's phyloglm rds 
# files and obtain model statistics! 

# this script can be deleted later, once results are documented.

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

#######################################################
###################### BODY MASS #################
#######################################################


######################## UN AND BODY MASS ######################## 

# load model
phyglm_UN_Mass_fix <- readRDS(here("Coastal Bird Phyloglm Models", "phyglm_UN_Mass_fix.rds"))

# as alpha is at upper bounds, we can also look at a regular logistic model
# use logistf package for this
# this runs as logistic regression with Firth's correction (a bias reduction method)
# Ives and Garland 2010 recommend log.alpha.bound = 4 as the limits (and this is default setting in phyloglm)
# they specify 4 because when the model reaches the upper bounds of this limit,
# the model estimates should be very similar to a logistic model using Firth's correction


# get values for results table for the fixed alpha model with scale(Mass_log)
summary(phyglm_UN_Mass_fix)
# use the bootstrapped CIs in the summary

# get alpha, t, and half life for the model
(phyglm_UN_Mass_fix$mean.tip.height) # this is t (mean tip height) of the tree
(alpha_mass <- phyglm_UN_Mass_fix$alpha) # this is alpha
(hl_mass <- log(2)/alpha_mass) # this is the half-life for the model
# compare the value for half life with the mean tip height of the tree
# compared to t, the half life is small -> therefore we conclude there is low phylogenetic signal


#######################################################
###################### SENSORY TRAITS #################
#######################################################


######################## UN AND C.T. (Dim light vision) ######################## 

# load model
phyglm_UN_CT_scale <- readRDS(here("Coastal Bird Phyloglm Models", "phyglm_UN_CT_scale.rds"))

# get values for results table for the fixed alpha model with scale(Mass_log)
summary(phyglm_UN_CT_scale)
# use the bootstrapped CIs in the summary

# get alpha, t, and half life for the model
(alpha_CT <- phyglm_UN_CT_scale$alpha) # alpha
(phyglm_UN_CT_scale$mean.tip.height) # t (aka mean tip height)
(hl_CT<- log(2)/alpha_CT) # half-life
# compared to t, this is a small Half-Life -> low phylogenetic signal

######################## UN AND Peak Frequency ######################## 

# load model
phyglm_UN_pf_scale <- readRDS(here("Coastal Bird Phyloglm Models", "phyglm_UN_pf_scale.rds"))

summary(phyglm_UN_pf_scale)

# get alpha, t, and half life for the model
(phyglm_UN_pf_scale$mean.tip.height) # t
(alpha_pfreq <- phyglm_UN_pf_scale$alpha) # alpha
(hl_pfreq<- log(2)/alpha_pfreq) # half life
# low half life relative to mean tip height -> low phylogenetic signal




#######################################################
###################### DIET TRAITS #################
#######################################################


######################## UN AND % DIET INVERTEBRATES ######################## 

phyglm_UN_Invert_fix_4.05 <- readRDS(here("Coastal Bird Phyloglm Models", "phyglm_UN_Invert_fix.rds"))

summary(phyglm_UN_Invert_fix_4.05)

# get alpha, t, and half life for the model
# using phyglm_UN_Invert_fix_4.05 as final model
(phyglm_UN_Invert_fix_4.05$mean.tip.height) # t
(alpha_invert <- phyglm_UN_Invert_fix_4.05$alpha) # alpha
(hl_invert <- log(2)/alpha_invert) # half life
# small half life -> low phylogenetic signal


######################## UN AND % DIET VERTEBRATES ######################## 

# load model
phyglm_UN_Vert_fix <- readRDS(here("Coastal Bird Phyloglm Models", "phyglm_UN_Vert_fix.rds"))

summary(phyglm_UN_Vert_fix)

# get alpha, t, and half life for the model
(phyglm_UN_Vert_fix$mean.tip.height) # t
(alpha_vert <- phyglm_UN_Vert_fix$alpha) # alpha
(hl_vert <- log(2)/alpha_vert) # half life
# compared to t, this is a small half life


######################## UN AND % DIET PLANT/SEED ######################## 

# load model
phyglm_UN_PS_scale <- readRDS(here("Coastal Bird Phyloglm Models", "phyglm_UN_PS_scale.rds"))

summary(phyglm_UN_PS_scale)

# get alpha, t, and half life for the model
(phyglm_UN_PS_scale$mean.tip.height) # t
(alpha_PS <- phyglm_UN_PS_scale$alpha) # alpha
(hl_PS <- log(2)/alpha_PS) # half life
# small compared to t -> low phylogenetic signal


######################## UN AND % DIET FRUIT/NECTAR ######################## 

# load model
phyglm_UN_FN_scale <- readRDS(here("Coastal Bird Phyloglm Models", "phyglm_UN_FN_scale.rds"))

summary(phyglm_UN_FN_scale)

# get alpha, t, and half life for the model
(phyglm_UN_FN_scale$mean.tip.height) # t
(alpha_FN <- phyglm_UN_FN_scale$alpha) # alpha
(hl_FN <- log(2)/alpha_FN) # half life
# compared to t, this is a small half life



#######################################################
###################### LIFE HISTORY TRAITS ############
#######################################################


######################## UN and % BROOD VALUE ##########################

# load model
phyglm_UN_bv_scale <- readRDS(here("Coastal Bird Phyloglm Models", "phyglm_UN_bv_scale.rds"))

summary(phyglm_UN_bv_scale)

# get alpha, t, and half life for the model
(phyglm_UN_bv_scale$mean.tip.height) # t
(alpha_bv <- phyglm_UN_bv_scale$alpha) # alpha
(hl_bv <- log(2)/alpha_bv) # half-life
# small half life relative to t -> low phylogenetic signal


######################## UN and % CLUTCH SIZE ##########################

# load model
phyglm_UN_clutch_fix <- readRDS(here("Coastal Bird Phyloglm Models", "phyglm_UN_clutch_fix.rds"))
           
summary(phyglm_UN_clutch_fix)

# get alpha, t, and half life for the model
(phyglm_UN_clutch_fix$mean.tip.height) # t
(alpha_clutch <- phyglm_UN_clutch_fix$alpha) # alpha
(hl_clutch <- log(2)/alpha_clutch) # half life
# small half life relative to t -> low phylogenetic signal



######################## UN and % LONGEVITY ##########################

# load model
phyglm_UN_long_fix <- readRDS(here("Coastal Bird Phyloglm Models", "phyglm_UN_long_fix.rds"))

summary(phyglm_UN_long_fix)

# get alpha, t, and half life for the model
(phyglm_UN_long_fix$mean.tip.height) # t
(alpha_clutch <- phyglm_UN_long_fix$alpha) # alpha
(hl_clutch <- log(2)/alpha_clutch) # half life
# small half life relative to t -> low phylogenetic signal



######################## UN and % DEVELOPMENTAL MODE ##########################

# load model
phyglm_UN_develop_fix <- readRDS(here("Coastal Bird Phyloglm Models", "phyglm_UN_develop_fix.rds"))

summary(phyglm_UN_develop_fix)

# get alpha, t, and half life for the model
(phyglm_UN_develop_fix$mean.tip.height) # t
(alpha_clutch <- phyglm_UN_develop_fix$alpha) # alpha
(hl_clutch <- log(2)/alpha_clutch) # half life
# small half life relative to t -> low phylogenetic signal




#######################################################
###################### NESTING TRAITS ############
#######################################################

  ######################## UN and NEST SITE LOW ##########################

# load model
phyglm_UN_nest_low_fix <- readRDS(here("Coastal Bird Phyloglm Models", "phyglm_UN_nest_low_fix.rds"))

summary(phyglm_UN_nest_low_fix)

# get alpha, t, and half life for the model
(phyglm_UN_nest_low_fix$mean.tip.height) # t
(alpha_clutch <- phyglm_UN_nest_low_fix$alpha) # alpha
(hl_clutch <- log(2)/alpha_clutch) # half life
# small half life relative to t -> low phylogenetic signal


######################## UN and NEST SITE LOW # # # ONLY # # # ##########################

# load model
phyglm_UN_nest_low_only_scale <- readRDS(here("Coastal Bird Phyloglm Models", "phyglm_UN_nest_low_only_scale.rds"))

summary(phyglm_UN_nest_low_only_scale)

# get alpha, t, and half life for the model
(phyglm_UN_nest_low_only_scale$mean.tip.height) # t
(alpha_clutch <- phyglm_UN_nest_low_only_scale$alpha) # alpha
(hl_clutch <- log(2)/alpha_clutch) # half life
# small half life relative to t -> low phylogenetic signal


######################## UN and NEST SITE HIGH ##########################

# load model
phyglm_UN_nest_high_fix <- readRDS(here("Coastal Bird Phyloglm Models", "phyglm_UN_nest_high_fix.rds"))

summary(phyglm_UN_nest_high_fix)

# get alpha, t, and half life for the model
(phyglm_UN_nest_high_fix$mean.tip.height) # t
(alpha_clutch <- phyglm_UN_nest_high_fix$alpha) # alpha
(hl_clutch <- log(2)/alpha_clutch) # half life
# small half life relative to t -> low phylogenetic signal


######################## UN and NEST SITE HIGH # # # ONLY # # # ##########################

# load model
phyglm_UN_nest_high_only_scale <- readRDS(here("Coastal Bird Phyloglm Models", "phyglm_UN_nest_high_only_scale.rds"))

summary(phyglm_UN_nest_high_only_scale)

# get alpha, t, and half life for the model
(phyglm_UN_nest_high_only_scale$mean.tip.height) # t
(alpha_clutch <- phyglm_UN_nest_high_only_scale$alpha) # alpha
(hl_clutch <- log(2)/alpha_clutch) # half life
# small half life relative to t -> low phylogenetic signal


######################## UN and NEST STRUCTURE (ENCLOSED/OPEN) ##########################



# too small sample size, do not run. 



######################## UN and % NEST SAFETY ##########################

# load model
phyglm_UN_nest_safety_scale <- readRDS(here("Coastal Bird Phyloglm Models", "phyglm_UN_nest_safety_scale.rds"))

summary(phyglm_UN_nest_safety_scale)

# get alpha, t, and half life for the model
(phyglm_UN_nest_safety_scale$mean.tip.height) # t
(alpha_clutch <- phyglm_UN_nest_safety_scale$alpha) # alpha
(hl_clutch <- log(2)/alpha_clutch) # half life
# small half life relative to t -> low phylogenetic signal



#######################################################
###################### SEXUAL SELECTION TRAITS ############
#######################################################

######################## UN and DICHROMATISM (BRIGHTNESS) ##########################

# load model
phyglm_UN_brightness_scale <- readRDS(here("Coastal Bird Phyloglm Models", "phyglm_UN_brightness_scale.rds"))

summary(phyglm_UN_brightness_scale)

# get alpha, t, and half life for the model
(phyglm_UN_brightness_scale$mean.tip.height) # t
(alpha_clutch <- phyglm_UN_brightness_scale$alpha) # alpha
(hl_clutch <- log(2)/alpha_clutch) # half life
# small half life relative to t -> low phylogenetic signal


######################## UN and DICHROMATISM (HUE) ##########################

# load model
phyglm_UN_hue_scale <- readRDS(here("Coastal Bird Phyloglm Models", "phyglm_UN_hue_scale.rds"))

summary(phyglm_UN_hue_scale)

# get alpha, t, and half life for the model
(phyglm_UN_hue_scale$mean.tip.height) # t
(alpha_clutch <- phyglm_UN_hue_scale$alpha) # alpha
(hl_clutch <- log(2)/alpha_clutch) # half life
# small half life relative to t -> low phylogenetic signal


######################## UN and SEXUAL SELECTION INTENSITY ON MALES ##########################

# load model
phyglm_UN_ssm_fix <- readRDS(here("Coastal Bird Phyloglm Models", "phyglm_UN_ssm_fix.rds"))

summary(phyglm_UN_ssm_fix)

# get alpha, t, and half life for the model
(phyglm_UN_ssm_fix$mean.tip.height) # t
(alpha_clutch <- phyglm_UN_ssm_fix$alpha) # alpha
(hl_clutch <- log(2)/alpha_clutch) # half life
# small half life relative to t -> low phylogenetic signal


######################## UN and SEXUAL SELECTION INTENSITY ON FEMALES ##########################

# load model
phyglm_UN_ssf_fix <- readRDS(here("Coastal Bird Phyloglm Models", "phyglm_UN_ssf_fix.rds"))

summary(phyglm_UN_ssf_fix)

# get alpha, t, and half life for the model
(phyglm_UN_ssf_fix$mean.tip.height) # t
(alpha_clutch <- phyglm_UN_ssf_fix$alpha) # alpha
(hl_clutch <- log(2)/alpha_clutch) # half life
# small half life relative to t -> low phylogenetic signal



#######################################################
###################### SOCIAL TRAITS ############
#######################################################

######################## UN and TERRITORIALITY ##########################

# load model
phyglm_UN_territorial_scale <- readRDS(here("Coastal Bird Phyloglm Models", "phyglm_UN_territorial_scale.rds"))

summary(phyglm_UN_territorial_scale)

# get alpha, t, and half life for the model
(phyglm_UN_territorial_scale$mean.tip.height) # t
(alpha_clutch <- phyglm_UN_territorial_scale$alpha) # alpha
(hl_clutch <- log(2)/alpha_clutch) # half life
# small half life relative to t -> low phylogenetic signal

######################## UN and COOPERATIVE BREEDING ##########################


# too small sample size, do not run. 


