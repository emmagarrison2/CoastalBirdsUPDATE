#### OBJECTIVE - Calculate Phylogenetic Signal (lambda) of Discrete Traits (0,1)

# This will require importing joined rds files of Coastal sp, index values, and predictor trait groups.

# This value is important for... ?? ensuring that it is appropriate to use phylo linear model, 
# instead of regular linear model 

#load/install packages 

library(here)
library(tidyverse)
library(ape)
library(geiger)
library(picante)
library(phytools)

#list of discrete traits - 
#### developmental mode 
#### nest site low 
#### nest site high 
#### nest strategy (open or enclosed)
#### territoriality
#### cooperative breeding 




######################## Phylogenetic Signal (lambda) for NESTING TRAITS ######################## 

#import Coastal_Species_Sensory.rds 

Nest_Traits <- readRDS(here("Outputs", "Coastal_Species_Nest.rds"))
head(Nest_Traits)
colnames(Nest_Traits)

#should I be going off of... Species_Jetz, Species_BirdLife, or Species_eBird???? prob Jetz because I join it to Jetz tree.... 

#do some reformatting and re-naming... 
Nest_Traits <- Nest_Traits %>%
  rename(Species = Species_Jetz)
colnames(Nest_Traits)
#nice 

#now, reformat Species column so that it is in Aaaa_aaaa format! 
Nest_Traits <- Nest_Traits %>%
  mutate(Species = str_replace(Species, " ", "_"))
View(Nest_Traits)

### we need to trim Jetz tree that we imported above to get a tree that only contains coastal species

# create data frame with traits and species as the rownames
Nestspp <- Nest_Traits %>%
  column_to_rownames(var="Species")
str(Nestspp)


### load phylogenetic tree
jetztree <- read.tree(here("Data", "Jetz_ConsensusPhy.tre")) 
jetztree$tip.label # look at species name formatting for tree tips to confirm it matches formatting we are using - yes!

# check that the tree is ultrametric (do all the tree tips line up?)
is.ultrametric(jetztree)
#true! 

# check that all the species in our data are also present in the tree
check_nest <- name.check(jetztree, Nestspp) # there should be many more species in the tree than our data
check_nest

# trim the tree to match the data by dropping all the extra species identified in the previous step
jetztree_nest <-drop.tip(jetztree, check_nest$tree_not_data)

# check again whether tree has same species as data. Should say "OK"
name.check(jetztree_nest, Nestspp)
#OK 
jetztree_nest #807 which is good 

######
          # nest structure (open/enclosed) DISCRETE 

Nestspp$NestStr <- as.factor(Nestspp$NestStr)

# get vector of nest safety values with species names as rownames
nest_str_vect <- setNames(Nestspp$NestStr, rownames(Nestspp))
View(nest_str_vect)

# drop species with NA values
nest_str_vect[!is.na(nest_str_vect)]

# get measure of lambda for nest safety
fit <- fitDiscrete(phy = jetztree_nest, dat = nest_str_vect, model = "ER", transform = "lambda")
print(fit)
#saying that fitted 'lambda' at model parameter is 0.000 


#####
      # nest site low - DISCRETE 

Nestspp$NestSite_Low <- as.factor(Nestspp$NestSite_Low)

# get vector of nest safety values with species names as rownames
nest_low_vect <- setNames(Nestspp$NestSite_Low, rownames(Nestspp))

# drop species with NA values
nest_low_vect[!is.na(nest_low_vect)]

# get measure of lambda for nest safety
fit_low <- fitDiscrete(phy = jetztree_nest, dat = nest_low_vect, model = "ER", transform = "lambda")
print(fit_low)
#saying that fitted 'lambda' at model parameter is 0.000 




#####
# nest site high - DISCRETE 

Nestspp$NestSite_High <- as.factor(Nestspp$NestSite_High)

# get vector of nest safety values with species names as rownames
nest_high_vect <- setNames(Nestspp$NestSite_High, rownames(Nestspp))
is.atomic(nest_high_vect)
#TRUE


# drop species with NA values
nest_high_vect[!is.na(nest_high_vect)]

# get measure of lambda for nest safety
lambda_high <- fitDiscrete(phy = jetztree_nest, dat = nest_high_vect, model = "ER", transform = "lambda")
print(lambda_high)
#saying that fitted 'lambda' at model parameter is 0.000 




