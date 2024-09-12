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
    position = position_nudge(y = -0.1)
  ) +  
  scale_fill_manual(
    values = c("All" = "#FFD67B", "Coastal" = "#E48816")
  ) +
  labs(
    x = "UAI", 
    y = ""
  ) + 
  theme_classic() + 
  theme(
    axis.title = element_text(size = 12),  # Adjust axis title size
    legend.position = "top",             # Position the legend on the right
    legend.title = element_blank()         # Remove the legend title
  )

print(UAI_density)

saveRDS(UAI_density, here("Outputs", "UAI_Density_Comparison.rds"))






#######################################################################
#####density plot for MUTI 

#new df with reduced columns of AllIndexes_AllBirds

AllBirds_MUTI <- AllIndexes_AllBirds %>% 
  filter(!is.na(MUTIscore) & is.finite(MUTIscore))
nrow(AllBirds_MUTI) #431

AllBirds_MUTI_2 <- AllBirds_MUTI %>% 
  select(MUTIscore, Group)
colnames(AllBirds_MUTI_2)
nrow(AllBirds_MUTI_2) #431, looks good! 

#new df with reduced columnds of AllIndexesCoastal

CoastalBirds_MUTI <- AllIndexesCoastal %>% 
  filter(!is.na(MUTIscore) & is.finite(MUTIscore))
nrow(CoastalBirds_MUTI) #130

CoastalBirds_MUTI_2 <- CoastalBirds_MUTI %>% 
  select(MUTIscore, Group)
colnames(CoastalBirds_MUTI_2)
nrow(CoastalBirds_MUTI_2) #130, looks good! 

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
    position = position_nudge(y = -0.045)
  ) +  
  scale_fill_manual(
    values = c("All" = "olivedrab1", "Coastal" = "olivedrab")
  ) +
  labs(
    x = "MUTI", 
    y = ""
  ) + 
  theme_classic() + 
  theme(
    axis.title = element_text(size = 12),  # Adjust axis title size
    legend.position = "top",             # Position the legend on the right
    legend.title = element_blank()         # Remove the legend title
  )

print(MUTI_density)


saveRDS(MUTI_density, here("Outputs", "MUTI_Density_Comparison.rds"))


#######################################################################
#####density plot for UN 

#new df with reduced columns of AllIndexes_AllBirds

AllBirds_UN <- AllIndexes_AllBirds %>% 
  filter(!is.na(Urban) & is.finite(Urban))
nrow(AllBirds_UN) #533

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

###time to bind rows for UAI 

combined_UN_for_density <- bind_rows(AllBirds_UN_2, CoastalBirds_UN_2)
nrow(combined_UN_for_density) #662
533 + 129 # = 662
#correct number of rows 



#histogram to display density of UN values (Urban or Non-Urban)  
#comparing all_birds distribution and coastal_birds distribution 

UN_density <- ggplot() + 
  # Density plots according to groups 
  geom_histogram(
    data = combined_UN_for_density, 
    aes(x = Urban, fill = Group), 
    alpha = 0.7, 
    binwidth = 0.3
  ) +  
  scale_fill_manual(
    values = c("All" = "cadetblue1", "Coastal" = "cadetblue")
  )  + 
  scale_x_continuous(
    breaks = c(0, 1),         # Specify the breaks (original values on the x-axis)
    labels = c("Non-Urban", "Urban")  # Labels to display for each break
  ) +
  labs (
    #title = "Coastal Score Distribution", 
    x = "UN", 
    y = ""
  ) + 
  theme_classic() + 
  theme(
    axis.title = element_text(size = 12),   # Adjust axis title size
    axis.title.y = element_text(margin = margin(r = 10)) 
  )

print(UN_density)

saveRDS(UN_density, here("Outputs", "UN_Density_Comparison.rds"))




###################### ARRANGE ###################### 
###################### SIDE BY SIDE ###################### 

#arrange plots 
all_density <- grid.arrange(UAI_density, MUTI_density, UN_density, ncol=3, top= "Coastal Score Distributions", left = "Number of Species")

saveRDS(all_density, here("Outputs", "All_Density_Plots.rds"))

#?ggplot



###################################### ALTERNATIVE OPTION FOR UN ###################################### 

########################## percentages, stacked barplot for UN 


library(dplyr)


UN_stacked_barplot <- ggplot(combined_UN_for_density_summary, aes(x = Group, y = percentage, fill = Group_Urban)) + 
  geom_bar(stat = "identity", position = "fill", color = "black") + 
  scale_y_continuous(labels = scales::percent) +  # Show percentages on the y-axis
  # Custom colors for each Group & Urban combination
  scale_fill_manual(
    values = c(
      "All & 0" = "darkslategray3",   # Color for "All & Nonurban"
      "All & 1" = "darkslategray4",        # Color for "All & Urban"
      "Coastal & 0" = "deepskyblue1",  # Color for "Coastal & Nonurban"
      "Coastal & 1" = "deepskyblue4"        # Color for "Coastal & Urban"
    ),
    labels = c(
      "All & 0" = "Non-urban",    # Custom labels
      "All & 1" = "Urban",
      "Coastal & 0" = "Non-urban",
      "Coastal & 1" = "Urban"
    )
  ) + 
  labs(
    x = "", 
    y = "Percentage", 
    fill = ""  # Legend title
  ) + 
  guides(
    fill = guide_legend(nrow = 2)  # Split the legend into 2 rows
  ) + 
  theme_classic() + 
  theme(
    axis.title = element_text(size = 12),   # Adjust axis title size
    legend.position = "top",                # Position the legend at the top
    legend.title = element_text(size = 12), # Adjust legend title size
    legend.text = element_text(size = 8)   # Adjust legend text size
  )
print(UN_stacked_barplot)




#new arrangement of all the plots (using stacked, percentage barplot)

all_density_2 <- grid.arrange(UAI_density, MUTI_density, UN_stacked_barplot, ncol=3, top= "Coastal Score Distributions", left = "Number of Species")

saveRDS(all_density_2, here("Outputs", "All_Density_Plots_option2.rds"))


