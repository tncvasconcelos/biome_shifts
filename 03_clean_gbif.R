# rm(list=ls())
library(data.table)
library(maptools)
library(raster)
library(sp)
library(rgeos)
library(rworldmap)
data("wrld_simpl")
source("00_utility_functions.R")

#-----------------------------
# Load WCVP dataset
dist_sample <- read.table("wcvp/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("wcvp/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

#-----------------------------
# Merge them in one big table
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")

reference_table <- list.files("taxized_reference_tables", full.names = T)
reference_table <- reference_table[grep("/reference_table",reference_table)]
reference_table <- do.call(rbind, lapply(reference_table, read.csv))
#reference_table$gbif_name <- fix.names.taxize(reference_table$gbif_name)

# Reading gbif file
gbif_data <- fread("gbif_data/0214267-230224095556074.csv") # load the table you downloaded from GBIF

# Looking at the WCVP table and TDWG to clean GBIF points
#-----------------------------
# If local
path="wcvp/wgsrpd-master/level3/level3.shp"
#-----------------------------

twgd_data <- suppressWarnings(maptools::readShapeSpatial(path))
issues_to_remove <- read.csv("gbif_issues_to_remove.csv")

cleaned_points <- subset(gbif_data, gbif_data$scientificName!="")
cleaned_points <- FilterWCVP_genus(cleaned_points, all_vars, twgd_data)
# cleaned_points <- readRDS("cleaned_points1.Rdata")
# saveRDS(cleaned_points, file="cleaned_points1.Rdata")
cleaned_points <- subset(cleaned_points, cleaned_points$basisOfRecord == "PRESERVED_SPECIMEN")
for(issue_index in 1:nrow(issues_to_remove)) {
  cleaned_points <- subset(cleaned_points, !grepl(issues_to_remove$issues_to_remove[issue_index], cleaned_points$issue))
}

subset_reference_table <- subset(reference_table, reference_table$gbif_name %in% unique(cleaned_points$scientificName))
if(nrow(subset_reference_table)>0){
  cleaned_points <- FilterWCVP(cleaned_points, all_vars, subset_reference_table, twgd_data) # This will filter the GBIF points acording to WCVP for species
}
# Cleaning common problems:
cleaned_points <- RemoveNoDecimal(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
cleaned_points <- RemoveCentroids(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
cleaned_points <- RemoveDuplicates(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
#cleaned_points <- RemoveOutliers(cleaned_points, species="scientificName", lon="decimalLongitude", lat="decimalLatitude")
#cleaned_points <- RemoveZeros(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
#cleaned_points <- RemoveSeaPoints(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
#cleaned_points <- Thinning(cleaned_points, species="scientificName", lat = "decimalLatitude", lon="decimalLongitude", n = 3)

write.csv(cleaned_points, file="cleaned_points.csv", row.names=F)
 
#------------------------
#cleaned_points <- fread("cleaned_points.csv")

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


