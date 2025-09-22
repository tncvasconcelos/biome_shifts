# rm(list=ls())
library(data.table)
library(raster)
library(sp)
library(sf)
library(rworldmap)
#install.packages("taxize")
source("00_utility_functions.R")

#-----------------------------
# Reading gbif files
# gbif_data1 <- fread("gbif_data/0021035-250811113504898.csv") # load the table you downloaded from GBIF

#-----------------------------
# Load WCVP dataset
dist_sample <- read.table("wcvp/wcvp_distribution.csv", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("wcvp/wcvp_names.csv", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

#-----------------------------
# Merge them in one big table
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")

# reference_table <- list.files("taxized_reference_tables", full.names = T)
# reference_table <- reference_table[grep("/reference_table",reference_table)]
# reference_table <- do.call(rbind, lapply(reference_table, read.csv))

issues_to_remove <- read.csv("1_occurrence_data/gbif_issues_to_remove.csv")
#-----------------------------
# Looking at the POWO table and TDWG to clean occurrence points
path="wgsrpd-master/level3/level3.shp"

# Read shapefile using sf
twgd_data_sf <- st_read(path)

twgd_data <- as(twgd_data_sf, "Spatial")

#-----------------------------
all_gbif_files <- list.files("gbif_data")
all_gbif_files <- subset(all_gbif_files, !grepl("_thinned_cleaned_points", all_gbif_files))
labels <- gsub(".csv","", all_gbif_files)
for(i in 1:length(labels)) {
  gbif_data <- fread(paste0("gbif_data/",all_gbif_files[i])) # load the table you downloaded from GBIF
  
  cleaned_points <- subset(gbif_data, gbif_data$scientificName!="")
  cleaned_points <- FilterWCVP_genus(points=cleaned_points, all_vars, twgd_data)
  
  #cleaned_points <- subset(cleaned_points, cleaned_points$basisOfRecord == "PRESERVED_SPECIMEN")
  for(issue_index in 1:nrow(issues_to_remove)) {
    cleaned_points <- subset(cleaned_points, !grepl(issues_to_remove$issues_to_remove[issue_index], cleaned_points$issue))
  }
  
  # subset_reference_table <- subset(reference_table, reference_table$gbif_name %in% unique(cleaned_points$scientificName))
  # if(nrow(subset_reference_table)>0){
  #   cleaned_points <- FilterWCVP(points=cleaned_points, all_vars, subset_reference_table, twgd_data) # This will filter the GBIF points according to WCVP for species
  # }
  
  # Cleaning common problems:
  cleaned_points <- RemoveNoDecimal(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  #cleaned_points <- RemoveCentroids(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  cleaned_points <- RemoveDuplicates(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  #cleaned_points <- RemoveOutliers(cleaned_points, species="scientificName", lon="decimalLongitude", lat="decimalLatitude")
  #cleaned_points <- RemoveZeros(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  #cleaned_points <- RemoveSeaPoints(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
  cleaned_points <- Thinning(gbif_data, lat = "decimalLatitude", lon="decimalLongitude", n = 5)
  
  write.csv(cleaned_points, file=paste0("gbif_data/",labels[i],"_thinned_cleaned_points.csv"), row.names=F)
  
}

#all_vars$whole_name <- paste(all_vars$taxon_name, all_vars$taxon_authors) 

# names_to_solve <- gbif_data$scientificName
# sources <- taxize::gna_data_sources()
# wcvp_name <- taxize::gna_verifier(names_to_solve,data_sources=sources$id[sources$title == "The World Checklist of Vascular Plants"],best_match_only=TRUE)$currentName
# gbif_data$scientificName <- wcvp_name
#


#------------------------
cleaned_points <- fread("cleaned_points.csv")

length(unique(cleaned_points$species))

# all_species <- c()
# all_genera <- c()
# for(i in 1:length(all_cleaned_points)) {
#   all_species <- c(all_species, length(unique(all_cleaned_points[[i]]$scientificName)))
#   all_genera <- c(all_genera, length(unique(all_cleaned_points[[i]]$genus)))
# }
# sum(all_species)
# sum(all_genera)


# Plotting to inspect distributions
# {; for(family_index in 1:length(all_cleaned_points)) {
#   points_cleaned <- all_cleaned_points[[family_index]]
#   genera <- unique(points_cleaned$genus)
#   genera <- subset(genera, genera!="")
#   pdf(paste0("plots/", names(all_cleaned_points)[family_index], "_points.pdf"))
#   for(genus_index in 1:length(genera)){
#     tmp_subset <- as.data.frame(points_cleaned[points_cleaned$genus==genera[genus_index],])
#     coord <- tmp_subset[,c("decimalLongitude","decimalLatitude")]
#     coordinates(coord) <- ~ decimalLongitude + decimalLatitude
#     plot(wrld_simpl)
#     plot(coord, col="red", add=T)
#     title(genera[genus_index])
#     print(genus_index)
#   }
#   dev.off()
# }
#  beepr::beep("fanfare"); } 


