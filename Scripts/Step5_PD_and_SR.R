##### the objective of this script is to find the phylogenetic distance (PD) and species richness (SR) 
# of each urban tolerance index (UAI, MUTI, and UN)

#install/load appropriate packages 

library(here)
library(tidyverse)
library(ape)
library(geiger)
library(picante)
library(phytools)
library(phyr) 
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(colorspace)
library(treeio)

######################## Phylogenetic distance ########################
# in this section, we will use the pd function to get Phylogenetic distance (PD) and
# species richness (SR) for each Urban Tolerance Index.
# Additionally, we will count the number of distinct families in each urban tolerance index


# import coastal species list (these are joined with scores from all 3 indexes)
# trait data not needed, as PD is across the entire index 

Coastal <- read.csv(here("Data", "Coastal_Birds_List.csv"), header=T)
head(Coastal)
colnames(Coastal)
nrow(Coastal) #807

# Reformat Species_Jetz column so that it is in Aaaa_aaaa format
# Put into a new column called Species
Coastal <- Coastal %>%
  mutate(Species = str_replace(Species_Jetz, " ", "_"))
head(Coastal)

# get coastal UAI species
CoastalUAI <- Coastal %>% filter(!is.na(aveUAI))
nrow(CoastalUAI) # 798 species

# get coastal UN species
CoastalUN <- Coastal %>% filter(!is.na(Urban))
nrow(CoastalUN) # 129 species

# get coastal MUTI species
CoastalMUTI <- Coastal %>% filter(!is.na(MUTIscore))
nrow(CoastalMUTI) # 130 species

# import Jetz phylogeny

### load phylogenetic tree
jetztree <- read.tree(here("Data", "Jetz_ConsensusPhy.tre")) 
jetztree$tip.label # look at species name formatting for tree tips to confirm it matches formatting we are using - yes!

# check that the tree is ultrametric (do all the tree tips line up?)
is.ultrametric(jetztree)
#True! 

######
# reorganize the data
# the picante package which allows us to calculate phylogenetic distance (PD)
# they have the data organized a specific way and the following steps put our data in the same format
ForPD <- Coastal %>% dplyr::select(Species, aveUAI, Urban, MUTIscore) %>%
  rename(UAI = aveUAI, UN = Urban, MUTI = MUTIscore) %>%
  mutate(
    UAI = if_else(!is.na(UAI), 1, 0), # convert species with index to 1 and species without index to 0
    UN  = if_else(!is.na(UN), 1, 0),
    MUTI = if_else(!is.na(MUTI), 1, 0)) %>% 
  pivot_longer(!Species, names_to="UrbanIndex", values_to="Score") %>%
  pivot_wider(., names_from = Species, values_from = Score) %>%
  column_to_rownames(., var="UrbanIndex")

View(ForPD)

# prune the jetz phylogeny to get only coastal bird species
coastal_jetz <- prune.sample(ForPD, jetztree)
coastal_jetz #807 tips 

# make sure the list of species in ForPD and the species in the pruned tree are in the same order
ForPD <- ForPD[, coastal_jetz$tip.label]

# get Faith's phylogenetic distance for each group of species (one for UAI, UTI, and UN)
PD <- pd(ForPD, coastal_jetz) # pd function from picante package
PD
# The pd function returns two values for each community, the PD and the species richness (SR)
# Faith’s PD is defined as the total branch length spanned by the tree for the species in the group

#interesting that UN is less phylogenetically diverse than MUTI, although they have almost = # of sp 

# Examine several other phylogenetic species diversity metrics

# measures described in Helmus et al. 2007 Am Nat
# link to paper: https://www.journals.uchicago.edu/doi/10.1086/511334

# We can obtain psv and psr for our data (the authors also present an evenness measure that we can't calculate).
# These measures are explained in the paper as follows:

# The three metrics we present here are derived statistically by considering the 
# value of some unspecified neutral trait shared by all species in a community. 
# As this neutral trait evolves up a phylogenetic tree, 
# speciation occurs, and from this point forward, evolution proceeds independently along each phylogenetic lineage.
# Our metric of phylogenetic species variability (PSV) quantifies how phylogenetic relatedness 
# decreases the variance of this hypothetical unselected trait shared by all species in the community. 
# To calculate PSV, only information about the phylogenetic relatedness of species in a community is needed, 
# not information about any particular trait. 
# Nonetheless, framing this measure in the context of a hypothetical neutral trait gives a metric 
# that has not only an intuitive interpretation but also appealing statistical properties. 
# The second metric quantifies phylogenetic species richness (PSR) as the number of species in a community 
# multiplied by the community’s PSV. 
# This metric is directly comparable to the traditional metric of species richness but includes phylogenetic relatedness. 

# PSV = phylogenetic species variability
# Bound between zero and one
# approaches zero as relatedness of the species increases 
# therefore, higher values therefore reflect greater diversity
PSV <- psv(ForPD, coastal_jetz, compute.var=F)
PSV
#UAI - 0.8044951 
#UN - 0.7606829 
#MUTI - 0.7970667 

# PSR = phylogenetic species richness
# Can take on any value, but is directly comparable to actual species richness
# if PSR is much lower than the actual richness, it would indicate a community where species are closely related (i.e. congenerics)
PSR <- psr(ForPD, coastal_jetz, compute.var=F)
PSR
# PSR column gives the phylogenetic species richness and SR gives the actual species richness
#UAI - 641.98712 
#UN - 98.12809 
#MUTI - 103.61867 

#### how many bird families are there for each urban index?

# organize the data again to enable this
colnames(Coastal)

families <- Coastal %>% dplyr::select(Species, Family_Sci, aveUAI, Urban, MUTIscore) %>%
  rename(UAI = aveUAI, UN = Urban, MUTI = MUTIscore) %>%
  mutate(
    UAI = if_else(!is.na(UAI), 1, 0), # convert species with index to 1 and species without index to 0
    MUTI  = if_else(!is.na(MUTI), 1, 0),
    UN = if_else(!is.na(UN), 1, 0)) 

# now find the families associated with each urban tolerance index

#UAI
familiesUAI <- families %>% filter(UAI==1) %>% 
  distinct(Family_Sci) # this will print a list of families in UAI
nrow(familiesUAI) # gives number of families in UAI = 81

#UN
familiesUN <- families %>% filter(UN==1) %>% 
  distinct(Family_Sci) # this will print a list of families in UN
nrow(familiesUN) # gives number of families in UN = 24

#MUTI
familiesMUTI <- families %>% filter(MUTI==1) %>% 
  distinct(Family_Sci) # this will print a list of families in MUTI
nrow(familiesMUTI) # gives number of families in MUTI = 33


# combine the family counts with PD, PSV, PSR, and Species Richness (SR)
PD$Family <- as.vector(c(nrow(familiesUAI), nrow(familiesUN), nrow(familiesMUTI)))

phy_measures <- PD %>% rownames_to_column(., var="index") %>%
  left_join(., PSV) %>% left_join(., PSR) %>%
  select(Family, SR, PD, PSVs, PSR) %>%
  mutate(across(c("PSVs", "PSR")), (round(., 3)))
  
phy_measures


###############################################
###### Make phylogenetic tree figure #########

# coastal bird phylogeny created above
coastal_jetz

# import eBird taxonomy to add Order info for each species
ebird_tax <- read.csv(here("Data","ebird_taxonomy_v2022.csv"), header=T) 
colnames(ebird_tax)

# add order to Coastal bird species list
Coastal_order <- ebird_tax %>% 
  rename(Species_eBird = SCI_NAME) %>% 
  select(Species_eBird, ORDER1) %>%
  left_join(Coastal, .)

# are any species missing an order?
Coastal_order %>% filter(is.na(ORDER1)) # 1 species
# this species is in Passeriformes
# manually make this edit
Coastal_order$ORDER1[Coastal_order$Species_eBird == "Sericornis citreogularis"] <- "Passeriformes"
Coastal_order %>% filter(is.na(ORDER1)) %>% nrow() # this should now be equal to zero

# prep list of species and their orders to be joined to phylogeny
orders <- Coastal_order %>% select(Species, ORDER1) %>%
  rename(label=Species, order = ORDER1)
length(unique(orders$order)) # 23 avian orders are represented

# how many species per Order?
spp_order <- orders %>% group_by(order) %>% count()
print(spp_order, n=Inf)

# what percentage are Passeriformes?
orders %>% filter(order=="Passeriformes") %>% nrow()
98/807 # % 12% are passerines
# this means that 88% of the coastal species list are non-passerines

# extract phylogeny and convert into tibble
phy_tibble <- as_tibble(coastal_jetz)

# add order information to the phylogeny
phy_join <- left_join(phy_tibble, orders)

# find nodes for each Order
ordernode <- phy_join %>% group_by(order) %>% summarize(min=min(node), max=max(node)) %>% left_join(spp_order)
ordernode

# convert back into class phylo
new_tree <- as.phylo(phy_join) 

# create a list of order names that are associated with each tree tip label (species name) 
order_info <- split(phy_join$label, phy_join$order)
order_info

# update the tree with the order as the grouping info using list created in previous step
order_coastal_tree <- groupOTU(new_tree, order_info) 
# this will allow us to color branches of tree based on Order

# create a custom color palette for the tree
# specifically, I want each Order to alternate colors (between two shades of gray) as you move around the tree
# to do this, I had to figure out the organization of the avian Orders around the circular tree, so they are listed in order below
# 
scale_color_custom <- function(...){
  ggplot2:::manual_scale(
    'color', 
    values = setNames(rep(c("gray25", "gray50"), 12), # two shades of gray repeated 12 times each for 23 Orders
                      c("Charadriiformes","Passeriformes", "Psittaciformes",
                        "Falconiformes", "Accipitriformes","Cathartiformes",
                        "Coraciiformes", "Strigiformes", "Pelecaniformes", 
                        "Suliformes","Ciconiiformes", "Procellariiformes",
                        "Sphenisciformes","Gaviiformes","Gruiformes", 
                        "Cuculiformes","Podicipediformes", "Phoenicopteriformes", 
                        "Columbiformes","Phaethontiformes", "Caprimulgiformes", 
                        "Anseriformes", "Galliformes")), 
    ...
  )
}



circ_order <- ggtree(order_coastal_tree, layout='circular', aes(color=group)) + 
  guides(color="none") +
  scale_color_custom()
circ_order

# figuring out where some of the most species-rich orders are located on the tree
# using a node in the middle of the range of nodes for each order as location for label -> there is probably a more sophisticated way to do this
# look at object ordernode to see nodes for each Order
ggtree(order_coastal_tree , layout='circular', aes(color=group)) + 
  scale_color_custom() +
  guides(color="none") +
  geom_cladelab(node = 765, label = "Gruiformes", fontsize = 3, vjust =-1) + # 60 species
  geom_cladelab(node = 275, label = "Charadriiformes",   fontsize = 3) + # 228 species
  geom_cladelab(node = 62, label = "Anseriformes", fontsize = 3) + # 126 species
  geom_cladelab(node = 660, label = "Pelecaniformes", fontsize = 3, hjust=1, vjust=-1) + # 77 species
  geom_cladelab(node = 415, label = "Passeriformes", fontsize = 3, vjust = 1, hjust =1) + # 98 species
  geom_cladelab(node = 515, label = "Accipitriformes", fontsize = 3, hjust=1, vjust=2) + # 50 species
  geom_cladelab(node = 720, label = "Suliformes", fontsize = 3, hjust=0.8, vjust=-2) + # 38 species
  geom_cladelab(node = 480, label = "Coraciiformes", fontsize = 3, hjust=1) # 29 species
 
# get species values for each urban tolerance index
un_dat <- Coastal %>% column_to_rownames(., var="Species") %>% select(Urban)
uai_dat <- Coastal %>% column_to_rownames(., var="Species") %>% select(aveUAI)
muti_dat <- Coastal %>% column_to_rownames(., var="Species") %>% select(MUTIscore)

# begin to build plot
# start with UN
hcl_palettes(palette="Blues", n = 9, plot = TRUE)
sequential_hcl(9, "Blues")

# using two shades from palette above that are high contrast
UN <- gheatmap(circ_order, un_dat, offset=-1, width=.1, colnames =F) +
  scale_fill_manual(values=c("#7FABD3", "#273871"), name = "UN", na.translate = F, labels=c("Non-Urban", "Urban")) # use na.translate = F to not plot species with NAs
UN  


# add MUTI
# try the following palettes: Green-Yellow (not good), YlGn (okay), ag_GrnYl (good!), Greens 3 (okay)
p1 <- UN + new_scale_fill()
UN_MUTI <-gheatmap(p1, muti_dat, offset=12, width=.1, colnames = F) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl", name="MUTI", na.value="white") 
UN_MUTI

# add UAI
# could try the following palettes: Oranges (good), OrRd (too red), YlOrBr (ugly), YlOrRd (too much red), Peach (nope), Heat (too many colors) and/or OrYel (ok)
p2 <- UN_MUTI + new_scale_fill()
UN_MUTI_UAI <-gheatmap(p2, uai_dat, offset=25, width=.1, colnames = F) +
  scale_fill_continuous_sequential(palette = "Oranges", name="UAI", na.value="white") 
UN_MUTI_UAI

ggsave(here("Results", "UN_MUTI_UAI.png"), width=27, height = 16, units="cm", dpi=300)


