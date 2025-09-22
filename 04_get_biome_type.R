# rm(list=ls())
# setwd("~/Desktop/biome_shifts")

#--------------------------------------
library(data.table)
source("00_utility_functions.R")

# Load distribution data (thinned points)
all_thinned_points <- list.files("gbif_data", full.name=T)[grep("thinned_",list.files("gbif_data/"))]
all_points <- list()
for(i in 1:length(all_thinned_points)) {
  one_point_file <- all_thinned_points[i] 
  all_points[[i]] <- as.data.frame(fread(one_point_file))
}
thinned_points <- do.call(rbind,all_points)

all_trees <- load.trees("2_trees/")
treebank_info <- read.csv("treebank_info_simplified.csv")
all_trees <- subset(all_trees, names(all_trees) %in% treebank_info$label[treebank_info$large_and_well_sampled=="y"])

all_tips <- unname(unlist(lapply(all_trees,"[[", 4)))

#--------------------------------------
# Add information on the main biome type for each point and species
biomes_for_points <- localityToBiome(points=thinned_points, lat="lat",lon="lon")
biomes_for_points <- subset(biomes_for_points, biomes_for_points$species%in%all_tips)

biomes_for_species <- getBiomes(points=biomes_for_points, species="species") # we want summaries for each species
# save(biomes_for_species, file="1_occurrence_data/biomes_for_species.Rsave")
# biomes_for_species <- as.data.frame(biomes_for_species)
# write.csv(biomes_for_species, file="1_occurrence_data/biomes_for_species.csv")

# Binarizing biomes
closed_canopy <- c("Tropical & Subtropical Moist Broadleaf Forests", "Tropical & Subtropical Dry Broadleaf Forests", "Tropical & Subtropical Coniferous Forests", "Temperate Broadleaf & Mixed Forests", "Temperate Conifer Forests", "Boreal Forests/Taiga")
open_canopy <- c("Tropical & Subtropical Grasslands, Savannas & Shrubland", "Temperate Grasslands, Savannas & Shrublands", "Flooded Grasslands & Savannas", "Montane Grasslands & Shrublands", "Tundra","Deserts & Xeric Shrublands", "Mediterranean Forests, Woodlands & Scrub",  "Mangroves")

summary_biome <- data.frame(species=biomes_for_species$species,closed_canopy=NA,open_canopy=NA,main_habitat=NA)
for(species_index in 1:nrow(biomes_for_species)) {
  one_genus <- biomes_for_species$species[species_index]
  subset_biomes <- biomes_for_species[species_index,]
  closed_canopy_points <- sum(subset_biomes[which(names(subset_biomes)%in%closed_canopy)])
  open_canopy_points <- sum(subset_biomes[which(names(subset_biomes)%in%open_canopy)])
  n_total <- closed_canopy_points+open_canopy_points
  if(closed_canopy_points!=0){
    summary_biome[species_index,2] <- round(closed_canopy_points/n_total,2)
  } else {
    summary_biome[species_index,2] <- round(0,2)
  }
  if(open_canopy_points!=0){
    summary_biome[species_index,3] <- round(open_canopy_points/n_total,2)
  } else {
    summary_biome[species_index,3] <- round(0,2)
  }
  summary_biome[species_index,4] <- ifelse(summary_biome[species_index,3]>0.5,"open","closed")
}
write.csv(summary_biome, file="1_occurrence_data/summarized_biome.csv", row.names=F)

plot(sort(summary_biome$open_canopy))

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
write.csv(habitats, file="1_occurrence_data/summarized_biome.csv", row.names=F)

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
# closed       open widespread 
# 5372       3507       2072 

#--------------------------------------
#--------------------------------------
# Organizing datasets for MuHiSSE analyses

habitats <- read.csv("1_occurrence_data/summarized_biome.csv")

all_trees <- load.trees("2_trees/")
treebank_info <- read.csv("treebank_info_simplified.csv")
all_trees <- subset(all_trees, names(all_trees) %in% treebank_info$label[treebank_info$large_and_well_sampled=="y"])

colnames(habitats)[1] <- "species"
ntips_total <- data.frame(tree=names(all_trees),ntips_start=NA,ntips_end=NA)
for(i in 1:length(all_trees)) {
  habitats_subset <- subset(habitats, habitats$species %in% all_trees[[i]]$tip.label)
  if(nrow(habitats_subset)>0) {
    ntips_total[i,2] <- Ntip(all_trees[[i]])
    pruned_tree <- keep.tip(all_trees[[i]], habitats_subset$species)
    habitats_subset$species <- gsub(" ","_",habitats_subset$species)
    
    splitted <- strsplit(habitats_subset$species, "_")
    flag_genus_only <- c()
    for(j in 1:length(splitted)) {
      if(length(is.na(lapply(splitted[[j]], "[", 2)))==1) {
        flag_genus_only <- c(flag_genus_only, j)
      }
      if(grepl("^[[:upper:]]", splitted[[j]][2])) {
        flag_genus_only <- c(flag_genus_only, j)
      }
    }
    if(length(flag_genus_only)>0) {
      habitats_subset <- habitats_subset[-flag_genus_only,] 
    }
    
    save(habitats_subset, file=paste0("3_organized_datasets_muhisse/",names(all_trees)[i],"_area_score.Rsave"), row.names=F)
    save(pruned_tree, file=paste0("4_organized_trees_muhisse/",names(all_trees)[i],".Rsave"))
    
    ntips_total[i,3] <- Ntip(pruned_tree)    
  }
}

