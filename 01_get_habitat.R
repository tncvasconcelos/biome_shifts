
# WWFload is taken from speciesgeocodeR; all credit goes to the original authors
WWFload <- function(x = NULL) {
  if (missing(x)) {
    x <- getwd()
  }
  download.file("http://assets.worldwildlife.org/publications/15/files/original/official_teow.zip",
                destfile = file.path(x, "wwf_ecoregions.zip"), quiet=TRUE)
  unzip(file.path(x, "wwf_ecoregions.zip"), exdir = file.path(x, "WWF_ecoregions"))
  file.remove(file.path(x, "wwf_ecoregions.zip"))
  wwf <- maptools::readShapeSpatial(file.path(x, "WWF_ecoregions", "official",
                                              "wwf_terr_ecos.shp"))
  return(wwf)
}


localityToBiome <- function (points, lat="lat",lon="lon") {
  #colnames(points) <- c("acceptedScientificName","key","decimalLatitude","decimalLongitude","basisOfRecord","issues")
  cat("Getting biome from locality data...")
  points[,lat] <-  as.numeric(points[,lat])
  points[,lon] <-  as.numeric(points[,lon])
  locations.spatial <- sp::SpatialPointsDataFrame(coords=points[,c(which(colnames(points)==lon), which(colnames(points)==lat))], data=points)
  wwf <- WWFload(tempdir())
  mappedregions <- sp::over(locations.spatial, wwf)
  biomes <- c("Tropical & Subtropical Moist Broadleaf Forests", "Tropical & Subtropical Dry Broadleaf Forests", "Tropical & Subtropical Coniferous Forests", "Temperate Broadleaf & Mixed Forests", "Temperate Conifer Forests", "Boreal Forests/Taiga", "Tropical & Subtropical Grasslands, Savannas & Shrubland", "Temperate Grasslands, Savannas & Shrublands", "Flooded Grasslands & Savannas", "Montane Grasslands & Shrublands", "Tundra", "Mediterranean Forests, Woodlands & Scrub", "Deserts & Xeric Shrublands", "Mangroves")
  points$eco_name <- mappedregions$ECO_NAME
  points$biome <- biomes[mappedregions$BIOME]
  return(points)
}


# getting biomes for each species
getBiomes <- function (points, species="species") {
  cat("Summarizing biome from locality data...")
  points <- as.data.frame(points) # not sure how to do it without transforming back to data.frame
  points <- subset(points, !is.na(points[,"biome"]))
  categories <- unique(points[,"biome"])
  taxa <- as.character(unique(points[,species]))
  result <- matrix(0, nrow=length(taxa), ncol=length(categories))
  rownames(result) <- taxa
  colnames(result) <- categories
  for (taxon_index in seq_along(taxa)) {
    for (category_index in seq_along(categories)) {
      x0 <- points[,species]==taxa[taxon_index]
      x1 <- points[,"biome"]==categories[category_index]
      result[taxon_index, category_index] <- length(which(x0 & x1))
    }
    cat(taxon_index, "\r")
  }
  return(result)
}

# rm(list=ls())
setwd("~/Desktop/biome_shifts/")
#################################################################################################
#--------------------------------------------
#--------------------------------------------
library(data.table)
library(maptools)
data("wrld_simpl")

# Load distribution data

all_point_files <- list.files("1_thinned_points", full.names = T)
thinned_points <- as.data.frame(fread(all_point_files))

#--------------------------------------
  # add information on main biome
  biomes_for_points <- localityToBiome(thinned_points, lat="lat",lon="lon")
  biomes_for_species <- getBiomes(biomes_for_points, species="species") # we want summaries for each species
  biomes_for_genera <- as.data.frame(biomes_for_genera)
  
  closed_canopy <- c("Tropical & Subtropical Moist Broadleaf Forests", "Tropical & Subtropical Dry Broadleaf Forests", "Tropical & Subtropical Coniferous Forests", "Temperate Broadleaf & Mixed Forests", "Temperate Conifer Forests", "Boreal Forests/Taiga")
  
  open_canopy <- c("Tropical & Subtropical Grasslands, Savannas & Shrubland", "Temperate Grasslands, Savannas & Shrublands", "Flooded Grasslands & Savannas", "Montane Grasslands & Shrublands", "Tundra","Deserts & Xeric Shrublands", "Mediterranean Forests, Woodlands & Scrub",  "Mangroves")
  
  summary_biome <- data.frame(genera=row.names(biomes_for_genera),closed_canopy=NA,open_canopy=NA,main_habitat=NA)
  for(genus_index in 1:nrow(biomes_for_genera)) {
    one_genus <- rownames(biomes_for_genera)[genus_index]
    subset_biomes <- biomes_for_genera[genus_index,]
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
  family <- rep(labels[family_index], nrow(summary_biome))
  all_results[[family_index]] <- cbind(family, summary_biome)

all_results <- do.call(rbind, all_results)
write.csv(all_results, file="datasets/myrtales_habitat.csv", row.names=F)
#View(all_results)
#table(all_results$main_habitat)


# adding tables for sensitivity analyses:
habitats <- read.csv("datasets/myrtales_habitat.csv")
habitats$cutoff_40 <- NA
habitats$cutoff_60 <- NA

for(i in 1:nrow(habitats)) {
  habitats$cutoff_40[i] <- ifelse(habitats$open_canopy[i]>0.4,"open","closed")
  habitats$cutoff_60[i] <- ifelse(habitats$open_canopy[i]>0.6,"open","closed")
}
write.csv(habitats, "datasets/myrtales_habitat.csv", row.names=F)

