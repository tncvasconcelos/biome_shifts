# rm(list=ls())
library(data.table)
library(raster)
library(sp)
library(sf)
library(rworldmap)
#install.packages("taxize")
source("00_utility_functions.R")

#-----------------------------
# Load WCVP dataset (needed for filtering)
dist_sample <- read.table("wcvp/wcvp_distribution.csv", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("wcvp/wcvp_names.csv", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

#-----------------------------
# Merge them in one big table
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")

# GBIF issues to remove
issues_to_remove <- read.csv("1_occurrence_data/gbif_issues_to_remove.csv")
#-----------------------------
# Looking at the POWO table and TDWG to clean occurrence points
path="wgsrpd-master/level3/level3.shp"
twgd_data_sf <- st_read(path)
twgd_data <- as(twgd_data_sf, "Spatial")

#-----------------------------
# Reading GBIF files (obs: gbif_data folder is in gitignore because the files are massive)
all_gbif_files <- list.files("gbif_data")
all_gbif_files <- subset(all_gbif_files, !grepl("_thinned_cleaned_points", all_gbif_files))
labels <- gsub(".csv","", all_gbif_files)
for(i in 1:length(labels)) {
  gbif_data <- fread(paste0("gbif_data/",all_gbif_files[i])) # load the table you downloaded from GBIF
  
  cleaned_points <- subset(gbif_data, gbif_data$scientificName!="")
  cleaned_points <- FilterWCVP_genus(points=cleaned_points, all_vars, twgd_data)
  
  for(issue_index in 1:nrow(issues_to_remove)) {
    cleaned_points <- subset(cleaned_points, !grepl(issues_to_remove$issues_to_remove[issue_index], cleaned_points$issue))
  }
  
  # Cleaning common problems:
  cleaned_points <- RemoveNoDecimal(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  cleaned_points <- RemoveDuplicates(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  # Thinning to eave only 5 points per 1x1 grid cell
  cleaned_points <- Thinning(gbif_data, lat = "decimalLatitude", lon="decimalLongitude", n = 5)
  
  write.csv(cleaned_points, file=paste0("gbif_data/",labels[i],"_thinned_cleaned_points.csv"), row.names=F)
}


