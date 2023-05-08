#rm(list=ls())

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
  cat("\n")
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

thinned_points <- as.data.frame(fread("cleaned_points.csv"))

#--------------------------------------
  # add information on main biome
  biomes_for_points <- localityToBiome(thinned_points, lat="lat",lon="lon")
  biomes_for_species <- getBiomes(biomes_for_points, species="species") # we want summaries for each species
  #save(biomes_for_species, file="biomes_for_species.Rsave")
  
  biomes_for_species <- as.data.frame(biomes_for_species)
  write.csv(biomes_for_species, file="biomes_for_species.csv", row.names=F)
  
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
  
write.csv(summary_biome, file="summarized_habitat.csv", row.names=F)

#View(all_results)
#table(all_results$main_habitat)


# adding tables for sensitivity analyses:
habitats <- read.csv("summarized_habitat.csv")
habitats$cutoff_40 <- NA
habitats$cutoff_60 <- NA

for(i in 1:nrow(habitats)) {
  habitats$cutoff_40[i] <- ifelse(habitats$open_canopy[i]>0.4,"open","closed")
  habitats$cutoff_60[i] <- ifelse(habitats$open_canopy[i]>0.6,"open","closed")
}

# adding columns for species that share biomes
#habitats <- read.csv("summarized_habitat.csv")
habitats$area_open <- NA
habitats$area_closed <- NA

for(i in 1:nrow(habitats)) {
  habitats$area_open[i] <- ifelse(habitats$open_canopy[i]>0.25,1,0)
  habitats$area_closed[i] <- ifelse(habitats$closed_canopy[i]>0.25,1,0)
}

#table(rowSums(habitats[,c(5,6)]))
#   1     2 
# 27219  6416 
# 27219 endemic to one habitat
# 6416 occurring in both habitats

write.csv(habitats, "summarized_habitat.csv", row.names=F)

