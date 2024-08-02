

UAI_eBirdJetzBirdLife_final

################# JOIN BODY MASS #################
# body mass from avonet (Tobias et al. 2022 Ecol Letters)


######################### JOIN SENSORY TRAITS ####################### 
# sensory traits: dim light vision and dominant vocal frequency
## Why ARE WE USING HU & CARDOSO 2009 FOR VOCAL FREQUENCY?

####### Dim light vision (eye CT ratio) ######
# uses BirdTree/Jetz taxonomy

####### Dominant vocal frequency ######
# Hu and Cardoso 2009


########## JOIN DIET TRAITS ##########
# diet traits: % invertebrates, % vertebrates, % plant/seed, % fruit/nectar
# from Wilman et al. 2014 (elton traits)

# import elton traits
# uses BirdTree/Jetz taxonomy






########## JOIN LIFE HISTORY TRAITS ##########
# life history traits include: clutch size, longevity, brood value, developmental mode

####### Clutch Size #######
# from Myhrvold
# may need to use both BirdLife and BirdTree/Jetz

####### Longevity #######
# from Bird et al. 2020
# uses BirdLife

####### Brood Value #######
# calculate most species manually using traits in Bird and Myhrvold
# fix one sketchy species
# fill in additional from Sayol et al. 2020. I think this works with Jetz but need to check


####### Developmental Mode #######
# from Wang and Kimball 2016
# which taxonomy?


########## JOIN NEST TRAITS ##########
# nest traits include: nest site (low, high), nest structure (open or enclosed), nest safety

# import nest site and structure info from Chia et al 2023
# uses the BirdLife taxonomy
nests_chia <- read.csv("~/Desktop/MAPS Project/Trait Files/nest.csv", header=T)
head(nests_chia)

nests_names <- nests_chia %>%
  rename(Species_BirdLife = Scientific_name)

UAI_nests <- left_join(UAI_eBirdJetzBirdLife_final, nests_names)


## simplify the nest traits




########## JOIN SEXUAL SELECTION TRAITS ##########
# sexual selection traits include: plumage brightness and hue, intensity of sexual selection on Males and Females