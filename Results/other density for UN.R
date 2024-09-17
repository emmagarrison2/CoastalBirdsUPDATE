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

View(CoastalBirds_UN_2)

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
    breaks = c(0, 0.9),         # Specify the breaks (original values on the x-axis)
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

UN_density + annotate("text", x = 0.9, y = 24, label = "Some text")
UN_density + annotate("text", x = 0.9, y = 150, label = "Some text")

UN_density + annotate("text", x = c(0.9, 0.9), y = c(24, 150), label = c("my label", "label 2"))

saveRDS(UN_density, here("Outputs", "UN_Density_Comparison.rds"))




###################### ARRANGE ###################### 
###################### SIDE BY SIDE ###################### 

#arrange plots 
all_density <- grid.arrange(UAI_density, MUTI_density, UN_density, ncol=3, top= "Coastal Score Distributions", left = "Number of Species")

saveRDS(all_density, here("Outputs", "All_Density_Plots.rds"))

#?ggplot


