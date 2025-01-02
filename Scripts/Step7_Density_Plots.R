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


if(!require(dplyr)){
  install.packages("dplyr")
  require(dplyr)
}
library(dplyr)

if(!require(patchwork)){
  install.packages("patchwork")
  require(patchwork)
}
library(patchwork)



#pull in refined dataset of coastal species and index scores

AllIndexesCoastal <- readRDS(here("Data", "Coastal_Birds_List.rds"))
str(AllIndexesCoastal)
AllIndexesCoastal$Urban <- ifelse(AllIndexesCoastal$Urban == "U", 1, 0)
colnames(AllIndexesCoastal)
View(AllIndexesCoastal)
nrow(AllIndexesCoastal) #807

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


#all-birds .rds for comparison against Coastal Birds Index Densities 

AllIndexes_AllBirds <- readRDS(here("Outputs", "UAI_MUTI_UN_final.rds"))
str(AllIndexes_AllBirds)
AllIndexes_AllBirds$Urban <- ifelse(AllIndexes_AllBirds$Urban == "U", 1, 0)
colnames(AllIndexes_AllBirds)
nrow(AllIndexes_AllBirds)#4433




#######################################################################
#####density plot for aveUAI 

# refine all birds to just be UAI score, species, and group (N/A)
# change all "NA" to = all
# refine coastal birds just be UAI score, species, and group (coastal = yes)
# change all "yes" to = coastal 
# bind rows, which just stacks the rows from coastal underneath the rows for all. 
# then, do density plot 

#new df with reduced columns of AllIndexes_AllBirds

nrow(AllIndexes_AllBirds) #4433
colnames(AllIndexes_AllBirds)
AllIndexes_AllBirds$Group <- "All"
colnames(AllIndexes_AllBirds)
AllBirds_UAI <- AllIndexes_AllBirds %>% 
  filter(!is.na(aveUAI) & is.finite(aveUAI))
nrow(AllBirds_UAI) #4347

#average all-birds UAI 
mean(AllBirds_UAI$aveUAI) #1.182

AllBirds_UAI_2 <- AllBirds_UAI %>% 
  select(aveUAI, Group)
colnames(AllBirds_UAI)
nrow(AllBirds_UAI) #4347, looks good! 

#new df with reduced columnds of AllIndexesCoastal

nrow(AllIndexesCoastal) #807
colnames(AllIndexesCoastal)
AllIndexesCoastal$Group <- "Coastal"
colnames(AllIndexesCoastal)
View(AllIndexesCoastal)
CoastalBirds_UAI <- AllIndexesCoastal %>% 
  filter(!is.na(aveUAI) & is.finite(aveUAI))
nrow(CoastalBirds_UAI) #798

#average coastal-birds UAI 
mean(CoastalBirds_UAI$aveUAI) #1.428

CoastalBirds_UAI_2 <- CoastalBirds_UAI %>% 
  select(aveUAI, Group)
colnames(CoastalBirds_UAI)
nrow(CoastalBirds_UAI_2) #798, looks good! 

###time to bind rows for UAI 

combined_UAI_for_density <- bind_rows(AllBirds_UAI_2, CoastalBirds_UAI_2)
nrow(combined_UAI_for_density) #5145
4347 + 798 # = 5145
#correct number of rows 


# Create the plot
UAI_density <- ggplot() + 
  # Density plots according to groups 
  geom_density(
    data = combined_UAI_for_density, 
    aes(x = aveUAI, fill = Group), 
    alpha = 0.5, 
    adjust = 0.5, 
  ) + 
  # Boxplot for both groups i believe? 
  geom_boxplot(
    data = combined_UAI_for_density,
    aes(x = aveUAI, fill = Group), 
    width = 0.12, 
    color = 'black',
    alpha = 0.5, 
    position = position_nudge(y = -0.1), 
    show.legend = F
  ) +  
  scale_fill_manual(
    values = c("All" = "#FFCDA1", "Coastal" = "#ED6F00")
  ) +
  labs(
    x = "UAI", 
    y = "Density"
  ) + 
  theme_classic() + 
  theme(
    axis.title.x = element_text(size = 14, margin = margin(t = 8)),  # Adjust axis title size
    axis.title.y = element_text(margin = margin(r = 10), size = 11), 
    legend.position = c(.95, .85),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),            
    legend.title = element_blank(),         # Remove the legend title
    legend.text = element_text(size = 10)
  )

print(UAI_density)

saveRDS(UAI_density, here("Results", "UAI_Density_Comparison.rds"))





#######################################################################
#####density plot for MUTI 

#new df with reduced columns of AllIndexes_AllBirds

AllBirds_MUTI <- AllIndexes_AllBirds %>% 
  filter(!is.na(MUTIscore) & is.finite(MUTIscore))
nrow(AllBirds_MUTI) #431

View(AllBirds_MUTI)

AllBirds_MUTI_2 <- AllBirds_MUTI %>% 
  select(MUTIscore, Group)
colnames(AllBirds_MUTI_2)
nrow(AllBirds_MUTI_2) #431, looks good! 

#average all-birds MUTI 
mean(AllBirds_MUTI_2$MUTIscore) #-0.00511

#new df with reduced columnds of AllIndexesCoastal

CoastalBirds_MUTI <- AllIndexesCoastal %>% 
  filter(!is.na(MUTIscore) & is.finite(MUTIscore))
nrow(CoastalBirds_MUTI) #130

CoastalBirds_MUTI_2 <- CoastalBirds_MUTI %>% 
  select(MUTIscore, Group)
colnames(CoastalBirds_MUTI_2)
nrow(CoastalBirds_MUTI_2) #130, looks good! 

#average coastal-birds MUTI 
mean(CoastalBirds_MUTI_2$MUTIscore) #-0.01237143

###time to bind rows for UAI 

combined_MUTI_for_density <- bind_rows(AllBirds_MUTI_2, CoastalBirds_MUTI_2)
nrow(combined_MUTI_for_density) #561
431 + 130 # = 561
#correct number of rows 

# Create the plot
MUTI_density <- ggplot() + 
  # Density plots according to groups 
  geom_density(
    data = combined_MUTI_for_density, 
    aes(x = MUTIscore, fill = Group), 
    alpha = 0.5, 
    adjust = 0.5, 
  ) + 
  # Boxplot for both groups i believe? 
  geom_boxplot(
    data = combined_MUTI_for_density,
    aes(x = MUTIscore, fill = Group), 
    width = 0.054, 
    color = 'black',
    alpha = 0.5, 
    position = position_nudge(y = -0.045), 
    show.legend = F
  ) +  
  scale_fill_manual(
    values = c("All" = "#B5DD60", "Coastal" = "#17A77E")
  ) +
  labs(
    x = "MUTI", 
    y = "Density"
  ) + 
  theme_classic() + 
  theme(
    axis.title.x = element_text(size = 14, margin = margin(t = 8)),  # Adjust axis title size
    axis.title.y = element_text(margin = margin(r = 13), size = 11), 
    legend.position = c(.95, .85),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),            
    legend.title = element_blank(), # Remove the legend title
    legend.text = element_text(size = 10)
  )

print(MUTI_density)


saveRDS(MUTI_density, here("Results", "MUTI_Density_Comparison.rds"))



###################################### UN DENSITY PLOT ###################################### 

########################## percentages, stacked barplot for UN 


#new df with reduced columns of AllIndexes_AllBirds

AllBirds_UN <- AllIndexes_AllBirds %>% 
  filter(!is.na(Urban) & is.finite(Urban))
nrow(AllBirds_UN) #533
View(AllBirds_UN)
colnames(AllBirds_UN)

AllBirds_UN_2 <- AllBirds_UN %>% 
  select(Urban, Group)
colnames(AllBirds_UN_2)
nrow(AllBirds_UN_2) #533, looks good! 

#new df with reduced columnds of AllIndexesCoastal

CoastalBirds_UN <- AllIndexesCoastal %>% 
  filter(!is.na(Urban) & is.finite(Urban))
nrow(CoastalBirds_UN) #129

CoastalBirds_UN_2 <- CoastalBirds_UN %>% 
  select(Urban, Group)
colnames(CoastalBirds_UN_2)
nrow(CoastalBirds_UN_2) #129, looks good! 

View(CoastalBirds_UN_2)


#how many birds are urban? 
sum(CoastalBirds_UN_2$Urban == 1, na.rm = TRUE) #41 
41/129 # = 0.318 --> 31.8% of coastal species urban 


###time to bind rows for UAI 

combined_UN_for_density <- bind_rows(AllBirds_UN_2, CoastalBirds_UN_2)
nrow(combined_UN_for_density) #662
533 + 129 # = 662
#correct number of rows 



combined_UN_for_density_summary <- combined_UN_for_density %>%
  group_by(Group, Urban) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(Group) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  mutate(Group_Urban = paste(Group, Urban, sep = " & "))  # Combine Group and Urban for legend

View(combined_UN_for_density_summary)

library(dplyr)

UN_stacked_barplot <- ggplot(combined_UN_for_density_summary, aes(x = Group, y = percentage, fill = Group_Urban)) + 
  geom_bar(stat = "identity", position = "fill", color = "black") + 
  scale_y_continuous(labels = scales::percent) +  # Show percentages on the y-axis
  # Custom colors for each Group & Urban combination
  scale_fill_manual(
    values = c(
      "All & 0" = "#DEEEF7",   # Color for "All & Nonurban"
      "All & 1" = "#A1C4E0",        # Color for "All & Urban"
      "Coastal & 0" = "#7FABD3",  # Color for "Coastal & Nonurban"
      "Coastal & 1" = "#5C90C6"        # Color for "Coastal & Urban"
    )
  ) + 
  labs(
    x = "UN", 
    y = "Percentage", 
    fill = ""  # Legend title
  ) + 
  theme_classic() + 
  theme(
    axis.title.x = element_text(size = 14, margin = margin(t = 8)), 
    axis.title.y = element_text(size = 11), 
    legend.position = "none"
  ) + 
  annotate(
    "text", 
    x = c(1, 1, 2, 2), 
    y = c(0.2, 0.7, 0.16, 0.65), 
    label = c("Urban", "Non-Urban", "Urban", "Non-Urban"),
    size = 4)

print(UN_stacked_barplot)

saveRDS(UN_stacked_barplot, here("Results", "UN_Density_Comparison.rds"))


######################################################################
######################################################################
######################################################################

#arrange all three plots side by side 

# hcl.colors(8, palette = "Blues 2")
#C6C8D7
#AEB2CD
#9299C2
#727EB5

# hcl.colors(8, palette = "Blues 3")
#BAD5FA
#91BAEB
#5E9BD8
#007BC0


# hcl.colors(8, palette = "Blues")
#B8D5E9
#93BADB
#6C9CCC
#417CBD

#let's try using the patchwork package... since grid.arrange may mess up the 
#formatting (UN doesn't have legend)

all_density_3 <- UAI_density + MUTI_density + UN_stacked_barplot + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 14))

print(all_density_3)

saveRDS (all_density_3, here("Results", "All_Density_Plots_edited.rds"))


