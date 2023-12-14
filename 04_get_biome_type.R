# rm(list=ls())
# setwd("biome_shifts")

#--------------------------------------
library(data.table)
library(maptools)
source("00_utility_functions.R")
data("wrld_simpl")

# Load distribution data (thinned points)
thinned_points <- as.data.frame(fread("cleaned_points.csv"))

#--------------------------------------
# Add information on the main biome type for each point and species
biomes_for_points <- localityToBiome(thinned_points, lat="lat",lon="lon")
biomes_for_species <- getBiomes(biomes_for_points, species="species") # we want summaries for each species
# save(biomes_for_species, file="1_occurrence_data/biomes_for_species.Rsave")
biomes_for_species <- as.data.frame(biomes_for_species)
#write.csv(biomes_for_species, file="1_occurrence_data/biomes_for_species.csv")

# Binarizing biomes
closed_canopy <- c("Tropical & Subtropical Moist Broadleaf Forests", "Tropical & Subtropical Dry Broadleaf Forests", "Tropical & Subtropical Coniferous Forests", "Temperate Broadleaf & Mixed Forests", "Temperate Conifer Forests", "Boreal Forests/Taiga")
open_canopy <- c("Tropical & Subtropical Grasslands, Savannas & Shrubland", "Temperate Grasslands, Savannas & Shrublands", "Flooded Grasslands & Savannas", "Montane Grasslands & Shrublands", "Tundra","Deserts & Xeric Shrublands", "Mediterranean Forests, Woodlands & Scrub",  "Mangroves")

summary_biome <- data.frame(genera=row.names(biomes_for_species),closed_canopy=NA,open_canopy=NA,main_habitat=NA)
for(genus_index in 1:nrow(biomes_for_species)) {
  one_genus <- rownames(biomes_for_species)[genus_index]
  subset_biomes <- biomes_for_species[genus_index,]
  closed_canopy_points <- sum(subset_biomes[which(names(subset_biomes)%in%closed_canopy)])
  open_canopy_points <- sum(subset_biomes[which(names(subset_biomes)%in%open_canopy)])
  n_total <- closed_canopy_points+open_canopy_points
  if(closed_canopy_points!=0){
    summary_biome[genus_index,2] <- round(closed_canopy_points/n_total,2)
  } else {
    summary_biome[genus_index,2] <- round(0,2)
  }
  if(open_canopy_points!=0){
    summary_biome[genus_index,3] <- round(open_canopy_points/n_total,2)
  } else {
    summary_biome[genus_index,3] <- round(0,2)
  }
  summary_biome[genus_index,4] <- ifelse(summary_biome[genus_index,3]>0.5,"open","closed")
}
# write.csv(summary_biome, file="1_occurrence_data/summarized_biome.csv", row.names=F)

#--------------------------------------
# Adding tables for sensitivity analyses:
habitats <- read.csv("1_occurrence_data/summarized_biome.csv")
habitats$cutoff_40 <- NA
habitats$cutoff_60 <- NA

for(i in 1:nrow(habitats)) {
  habitats$cutoff_40[i] <- ifelse(habitats$open_canopy[i]>0.4,"open","closed")
  habitats$cutoff_60[i] <- ifelse(habitats$open_canopy[i]>0.6,"open","closed")
}

#--------------------------------------
# Adding columns for species that share biomes (25% distribution)
habitats$area_open <- NA
habitats$area_closed <- NA

for(i in 1:nrow(habitats)) {
  habitats$area_open[i] <- ifelse(habitats$open_canopy[i]>0.25,1,0)
  habitats$area_closed[i] <- ifelse(habitats$closed_canopy[i]>0.25,1,0)
}
# write.csv(habitats, file="1_occurrence_data/summarized_biome.csv", row.names=F)

#--------------------------------------
summarized_habitat <- read.csv("1_occurrence_data/summarized_biome.csv")
summarized_habitat$area_final <- NA

for (i in 1:nrow(summarized_habitat)){
  if (summarized_habitat[i, "area_open"] == 1 & summarized_habitat[i, "area_closed"] == 1){
    summarized_habitat[i, "area_final"] = "widespread"
  }
  if (summarized_habitat[i, "area_open"] == 0 & summarized_habitat[i, "area_closed"] == 1){
    summarized_habitat[i, "area_final"] = "closed"
  }
  if (summarized_habitat[i, "area_open"] == 1 & summarized_habitat[i, "area_closed"] == 0){
    summarized_habitat[i, "area_final"] = "open"
  }
}

table(summarized_habitat$area_final)
#    closed       open widespread 
# 10636       4818       3921 
