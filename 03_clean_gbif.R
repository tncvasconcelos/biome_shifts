# rm(list=ls())
library(data.table)
library(maptools)
library(raster)
library(sp)
library(rgeos)
library(rworldmap)
data("wrld_simpl")


FilterWCVP_genus <- function(points, all_vars, twgd_data, lon="decimalLongitude", lat="decimalLatitude") {
  npoints_start <- nrow(points)
  tmp_points = as.data.frame(points)
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  tmp_points = subset(tmp_points, !is.na(tmp_points$x))
  tmp_points = subset(tmp_points, !is.na(tmp_points$y))
  # Load shape files and make sure they have the same name as the WCVP column with the TDWG areas
  #twgd_data <- suppressWarnings(maptools::readShapeSpatial(path))
  dubiousGBIF_ids <- c()
  all_genera <- unique(all_vars$genus)
  all_vars_genus_level <- all_vars[all_vars$taxon_rank=="Genus",]
  for(genus_index in 1:length(all_genera)) {
    gbif_subset <- subset(tmp_points, tmp_points$genus == all_genera[genus_index])
    if(nrow(gbif_subset)!=0) {
      wcvp_subset <- subset(all_vars_genus_level, all_vars_genus_level$genus == all_genera[genus_index])
      wcvp_subset <- subset(wcvp_subset, wcvp_subset$introduced==0)
      wcvp_subset <- subset(wcvp_subset, wcvp_subset$extinct==0)
      wcvp_subset <- subset(wcvp_subset, wcvp_subset$location_doubtful==0)
      occ_areas <- wcvp_subset$area_code_l3
      area_plus_buffer <- twgd_data[which(as.character(twgd_data$LEVEL3_COD) %in% occ_areas),]
      if(nrow(area_plus_buffer)>0) {
        coords <- gbif_subset[,c("x","y")]
        sp::coordinates(coords) <- ~ x + y
        answer <- which(is.na(sp::over(coords, area_plus_buffer)[,3]))
        if(length(answer) != 0) {
          dubiousGBIF_ids <- c(dubiousGBIF_ids, as.character(gbif_subset$gbifID[answer]))
        }
      }
    }
    cat(genus_index, "\r")
  }
  cleaned_points <- subset(points, !as.character(points$gbifID) %in% dubiousGBIF_ids)
  npoints_end <- nrow(cleaned_points)
  print(paste0(npoints_start - npoints_end, " points removed."))
  return(cleaned_points)
}


FilterWCVP <- function(points, all_vars, reference_table, twgd_data, species= "scientificName", lon="decimalLongitude", lat="decimalLatitude") {
  npoints_start <- nrow(points)
  tmp_points = as.data.frame(points)
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  tmp_points = subset(tmp_points, !is.na(tmp_points$x))
  tmp_points = subset(tmp_points, !is.na(tmp_points$y))
  # Load shape files and make sure they have the same name as the WCVP column with the TDWG areas
  #twgd_data <- suppressWarnings(maptools::readShapeSpatial(path))
  dubiousGBIF_ids <- c()
  for(species_index in 1:nrow(reference_table)) {
    gbif_subset <- subset(tmp_points, tmp_points$scientificName == reference_table$gbif_name[species_index])
    if(nrow(gbif_subset)!=0) {
      wcvp_subset <- subset(all_vars, all_vars$taxon_name == reference_table$wcvp_name[species_index])
      wcvp_subset <- subset(wcvp_subset, wcvp_subset$introduced==0)
      wcvp_subset <- subset(wcvp_subset, wcvp_subset$extinct==0)
      wcvp_subset <- subset(wcvp_subset, wcvp_subset$location_doubtful==0)
      occ_areas <- wcvp_subset$area_code_l3
      area_plus_buffer <- twgd_data[which(as.character(twgd_data$LEVEL3_COD) %in% occ_areas),]
      if(nrow(area_plus_buffer)>0) {
        coords <- gbif_subset[,c("x","y")]
        sp::coordinates(coords) <- ~ x + y
        answer <- which(is.na(sp::over(coords, area_plus_buffer)[,3]))
        if(length(answer) != 0) {
          dubiousGBIF_ids <- c(dubiousGBIF_ids, as.character(gbif_subset$gbifID[answer]))
        }
      }
    }
    cat(species_index, "\r")
  }
  cleaned_points <- subset(points, !as.character(points$gbifID) %in% dubiousGBIF_ids)
  npoints_end <- nrow(cleaned_points)
  print(paste0(npoints_start - npoints_end, " points removed."))
  return(cleaned_points)
}

# Getting climate means per TWGD region
MeanRasterWCVP <- function(path_raster="3_Landscape_instability/bio_1_instability.tif", path_tdwg="wgsrpd-master/level3/level3.shp") {
  twgd_data <- suppressWarnings(maptools::readShapeSpatial(path_tdwg))
  raster_example <- raster(path_raster)
  template_map <- raster_example
  template_map[!is.na(template_map[])] <- 0
  template_map <- aggregate(template_map, fact=25)
  raster_list <- list()
  for(area_index in 1:length(twgd_data)){
    one_area <- twgd_data[area_index,]
    cropped_raster <- NULL
    try(cropped_raster <- mask(crop(raster_example, one_area), one_area))
    if(!is.null(cropped_raster)) {
      mean_value <- mean(subset(cropped_raster[], !is.na(cropped_raster[])))
      template_raster <- cropped_raster
      template_raster[!is.na(template_raster[])] <- mean_value
      template_raster <- aggregate(template_raster, fact=25)
      template_raster <- raster::resample(template_raster, template_map)
      raster_list[[area_index]] <- raster::mask(template_raster, template_map) 
      print(area_index)
    } else { next }
  }
  raster_list[which(unlist(lapply(raster_list, is.null)))] <- NULL
  mm <- do.call(merge, raster_list)
  return(mm)
}

############################
#' Removes points in the sea
#' @param points A data.frame of distribution points with at least three columns where one column represents species names and the other two decimal coordinates
#' @param lat The name of the column in the data.frame with latitudes
#' @param lon The name of the column in the data.frame with longitudes
#' @param buffer A number of degrees around continental areas where points are still kept after filtering
RemoveSeaPoints <- function(points, lon="decimalLongitude", lat="decimalLatitude", buffer=0) {
  npoints_start <- nrow(points)
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  wrld_map <- rworldmap::getMap(resolution="low") # leaving both maps in the code for now, should probably drop one of them later
  coords <- tmp_points[,c("x","y")]
  sp::coordinates(coords) <- ~ x + y
  sp::proj4string(coords) <- sp::proj4string(wrld_map)
  country_plus_buffer <- raster::buffer(wrld_map, buffer) # adding buffer around polygons
  answer <- which(is.na(sp::over(coords, country_plus_buffer)))
  if(length(answer) > 0) {
    points <- points[-answer,]
    npoints_end <- nrow(points)
    print(paste0(npoints_start - npoints_end, " points removed."))
    return(points)
  } else {
    print("no points removed")
    return(points) }
}

#' Removes points that have 0 for both latitude and longitude
#' @param points A data.frame of distribution points with at least three columns where one column represents species names and the other two decimal coordinates
#' @param lat The name of the column in the data.frame with latitudes
#' @param lon The name of the column in the data.frame with longitudes
RemoveZeros <- function(points, lon="decimalLongitude", lat="decimalLatitude") {
  if(!inherits(points, "data.frame")) {
    stop("Argument points is not a data.frame.")
  }
  npoints_start <- nrow(points)
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  if(any(tmp_points$x==0 & tmp_points$y==0)) {
    points <- points[-which(tmp_points$x==0 & tmp_points$y==0),]
  }
  npoints_end <- nrow(points)
  print(paste0(npoints_start - npoints_end, " points removed."))
  return(points)
}

#' Removes outliers
#' @param points A data.frame of distribution points with at least three columns where one column represents species names and the other two decimal coordinates
#' @param species The name of the column in the data.frame with the names of species
#' @param lat The name of the column in the data.frame with latitudes
#' @param lon The name of the column in the data.frame with longitudes
RemoveOutliers <- function(points, species="scientificName", lon="decimalLongitude", lat="decimalLatitude") {
  npoints_start <- nrow(points)
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  colnames(tmp_points)[colnames(tmp_points)==species] <- "tmp_species"
  spp <- unique(tmp_points$species)
  all_points <- list()
  for(species_index in 1:length(spp)){
    sp0 <- tmp_points[tmp_points$species==spp[species_index],]
    out_lat <- grDevices::boxplot.stats(sp0$y)$out
    out_lon <- grDevices::boxplot.stats(sp0$x)$out
    sp0 <- sp0[!sp0$y %in% out_lat, ]
    sp0 <- sp0[!sp0$x %in% out_lon, ]
    all_points[[species_index]] <- sp0
  }
  points <- do.call(rbind, all_points)
  colnames(points)[colnames(points)=="x"] <- lon
  colnames(points)[colnames(points)=="y"] <- lat
  colnames(points)[colnames(points)=="tmp_species"] <- species
  npoints_end <- nrow(points)
  print(paste0(npoints_start - npoints_end, " points removed."))
  return(points)
}

#' Removes points that are located in the wrong country according to their GBIF labels
#' That will remove points that are not located in the countries where their labels say they were collected
#' @param points A data.frame of distribution points with at least five columns where one column represents species names and other two decimal coordinates.
#' @param lat The name of the column in the data.frame with latitudes
#' @param lon The name of the column in the data.frame with longitudes
#' @param buffer A number of degrees around each country where points are still considered part of that country
#' @details The input data.frame must have a column named countryCode and one named gbifID, as the .csv files downloaded directly from GBIF.
RemoveWrongCountries <- function(points, lon="decimalLongitude", lat="decimalLatitude", buffer=5, wrld_simpl="") {
  npoints_start <- nrow(points)
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  countries <- as.character(wrld_simpl[2]$ISO2)
  dubiousGBIF_ids <- c()
  for(country_index in 1:length(countries)) {
    tmp_country <- countries[country_index]
    if(length(which(tmp_points$countryCode %in% tmp_country)) > 0) {
      tmp_subset <- tmp_points[tmp_points$countryCode==tmp_country,]
      coords <- stats::na.omit(tmp_subset[,c("x","y")])
      sp::coordinates(coords) <- ~ x + y
      sp::proj4string(coords) <- sp::proj4string(wrld_simpl) <- "+proj=longlat +ellps=WGS84 +no_defs" 
      country_plus_buffer <- raster::buffer(wrld_simpl[country_index,], buffer) # adding buffer around country
      answer <- which(is.na(sp::over(coords, country_plus_buffer)))
      dubiousGBIF_ids <- c(dubiousGBIF_ids, tmp_subset$gbifID[answer])
    }
  }
  if(!is.null(dubiousGBIF_ids)) {
    points_cleaned <- points[-which(points$gbifID %in% dubiousGBIF_ids),]
  }
  npoints_end <- nrow(points_cleaned)
  print(paste0(npoints_start - npoints_end, " points removed."))
  return(points_cleaned)
}

#' Removes points that are located in country centroids
#' @param points A data.frame of distribution points with at least three columns where one column represents species names and the other two decimal coordinates
#' @param lat The name of the column in the data.frame with latitudes
#' @param lon The name of the column in the data.frame with longitudes
#' @param buffer A number in meters around each country centroid for points to be removed
RemoveCentroids <- function(points, lon="decimalLongitude", lat="decimalLatitude", buffer=75000) {
  npoints_start <- nrow(points)
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  wrld_map <- rworldmap::getMap(resolution="low")
  # here the buffer is in meters
  # Note: probably should do something about small countries where the buffer of 75km may be too broad
  coords <- tmp_points[,c("x","y")]
  sp::coordinates(coords) <- ~ x + y
  sp::proj4string(coords) <- sp::proj4string(wrld_map)
  centroids <- rgeos::gCentroid(wrld_map, byid=TRUE)
  centroids_plus_buffer <- raster::buffer(centroids, buffer) # adding buffer around polygons
  answer <- which(!is.na(sp::over(coords, centroids_plus_buffer)))
  if(length(answer) > 0) {
    points <- points[-answer,]
    npoints_end <- nrow(points)
    print(paste0(npoints_start - npoints_end, " points removed."))
    return(points)
  } else {
    print("no points removed")
    return(points) }
}

#' Removes duplicated latitudes and longitudes for the same species
#' @param points A data.frame of distribution points 
#' @param species The name of the column in the data.frame with the names of species
#' @param lat The name of the column in the data.frame with latitudes
#' @param lon The name of the column in the data.frame with longitudes
RemoveDuplicates <- function(points, species="scientificName", lon="decimalLongitude", lat="decimalLatitude") {
  npoints_start <- nrow(points)
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  colnames(tmp_points)[colnames(tmp_points)==species] <- "tmp_species"
  spp <- unique(tmp_points$species)
  all_points <- list()
  for(species_index in 1:length(spp)){
    tmp_subset <- as.data.frame(tmp_points[tmp_points$species==spp[species_index],])
    all_points[[species_index]] <- tmp_subset[-which(duplicated(tmp_subset[,c("x","y")])),]
  }
  points <- do.call(rbind, all_points)
  colnames(points)[colnames(points)=="x"] <- lon
  colnames(points)[colnames(points)=="y"] <- lat
  colnames(points)[colnames(points)=="tmp_species"] <- species
  npoints_end <- nrow(points)
  print(paste0(npoints_start - npoints_end, " points removed."))
  return(points)
}

#' Removes points with coordinates without decimal cases (probably innacurate)
#' @param points A data.frame of distribution points 
#' @param lat The name of the column in the data.frame with latitudes
#' @param lon The name of the column in the data.frame with longitudes
RemoveNoDecimal <- function(points, lon="decimalLongitude", lat="decimalLatitude") {
  npoints_start <- nrow(points)
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  length_decimal_lat <- nchar(sub("^[^.]*", "", tmp_points$y))
  length_decimal_lon <- nchar(sub("^[^.]*", "", tmp_points$x))
  points <- points[which(length_decimal_lat>=1 & length_decimal_lon>=1),]
  npoints_end <- nrow(points)
  print(paste0(npoints_start - npoints_end, " points removed."))
  return(points)
}


#' Thinning distribution data to smooth sampling bias
#' @param points A data.frame of three columns containing species and coordinates
#' @param species A character string indicating name of column with species names
#' @param lat A character string indicating name of column with latitudes
#' @param lon character string indicating name of column with longitudes
#' @param n A number indicating how many points to keep in each cell after thinning
Thinning <- function(points, species="species", lat = "decimalLatitude", lon="decimalLongitude", n = 3) {
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  spp <- unique(tmp_points[,species])
  results <- list()
  for(species_index in 1:length(spp)) {
    coords <- tmp_points[tmp_points[,species]==spp[species_index],c("y","x")]
    coords <- coords[!duplicated(coords[,"x"]) & !duplicated(coords[,"y"]),]
    if(nrow(coords) > 1) {
      sp::coordinates(coords) <- ~ y + x
      raster::crs(coords) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
      r0 <- raster::raster(coords)
      raster::res(r0) <- 1 # cell resolution
      r0 <- raster::extend(r0, raster::extent(r0) + 5) 
      res <- cbind(spp[species_index], as.data.frame(dismo::gridSample(coords, r0, n))) # n = maximum number of points per cell
      colnames(res) <- c("species", "lat","lon")
      results[[species_index]] <- res
    } else {
      res <- cbind(spp[species_index],coords)
      colnames(res) <- c("species", "lat","lon")
      results[[species_index]] <- res
    }
  }
  results <- do.call(rbind, results)
  return(results)
}

#-----------------------------
# If local
dist_sample <- read.table("wcvp/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("wcvp/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

#-----------------------------
# Merge them in one big table
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")

# reference table for taxized names
#-----------------------------
# If local
# reference_table <- list.files("taxized_reference_tables", full.names = T)
# reference_table <- do.call(rbind, lapply(reference_table, read.csv))

# species_list <- unique(all_vars$taxon_name)
# saveRDS(species_list, file="species_list.Rdata")
# taxized_names <- resolveGBIF(species_list) # This function adjust the names to the GBIF taxonomic backbone
# # Make sure WCVP and GBIF communicate
# reference_table <- data.frame(wcvp_name = species_list, gbif_name = taxized_names) # you will need this table later
# write.csv(reference_table, file="taxized_reference_tables/reference_table.csv", row.names = F) # saving table that you will need later

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


