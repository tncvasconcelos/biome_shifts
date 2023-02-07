# rm(list=ls())
setwd("~/Desktop/biome_shifts/")

library(ape)

# Load trees
all_trees_files <- list.files("4_organized_trees_geohisse/")
all_trees <- lapply(paste0("4_organized_trees_geohisse/",all_trees_files), read.tree)
names(all_trees) <- gsub(".tre","",all_trees_files)

# Load datasets
all_area_files <- list.files("3_organized_datasets_geohisse/")
all_areas <- lapply(paste0("3_organized_datasets_geohisse/",all_area_files), read.csv)
names(all_areas) <- gsub("_area_score.csv","",all_area_files)

# I still have to double check a lot of these because it looks like we lost a lot of data
# when filtering the points. Let's start with the following trees for a test:

pilot <- c("Acacia-Renner_et_al-2019","Lamiales-Fonseca-2021","Poaceae-Spriggs_et_al-2014",
"Carex-MartinBravo_et_al-2019" ,"Onagraceae-Freyman_&_Hohna-2018","Quercus-Hipp_et_al-2017",
"Viburnum-Landis_et_al-2021")

pilot_areas <- all_areas[pilot]
pilot_trees <- all_trees[pilot]


