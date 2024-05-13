# rm(list=ls())
setwd("/Users/tvasc/Desktop/biome_shifts")
# setwd("biome_shifts/")

library(ggplot2)
library(data.table)
library(maps)
library(ggthemes)
library(viridis)
library(sf)
library(rmapshaper)
library(maptools)

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

# Descriptive plots

wwf <- WWFload(tempdir())
wwf_biomes <- wwf

biomes <- c("Tropical & Subtropical Moist Broadleaf Forests", "Tropical & Subtropical Dry Broadleaf Forests", "Tropical & Subtropical Coniferous Forests", "Temperate Broadleaf & Mixed Forests", "Temperate Conifer Forests", "Boreal Forests/Taiga", "Tropical & Subtropical Grasslands, Savannas & Shrubland", "Temperate Grasslands, Savannas & Shrublands", "Flooded Grasslands & Savannas", "Montane Grasslands & Shrublands", "Tundra", "Mediterranean Forests, Woodlands & Scrub", "Deserts & Xeric Shrublands", "Mangroves")
wwf_biomes$biome <- biomes[wwf_biomes$BIOME]
closed_canopy <- c("Tropical & Subtropical Moist Broadleaf Forests", "Tropical & Subtropical Dry Broadleaf Forests", "Tropical & Subtropical Coniferous Forests", "Temperate Broadleaf & Mixed Forests", "Temperate Conifer Forests", "Boreal Forests/Taiga")
open_canopy <- c("Tropical & Subtropical Grasslands, Savannas & Shrubland", "Temperate Grasslands, Savannas & Shrublands", "Flooded Grasslands & Savannas", "Montane Grasslands & Shrublands", "Tundra","Deserts & Xeric Shrublands", "Mediterranean Forests, Woodlands & Scrub",  "Mangroves")
wwf_biomes$canopy <- NA
wwf_biomes$canopy[which(wwf_biomes$biome%in%closed_canopy)] <- "closed"
wwf_biomes$canopy[which(wwf_biomes$biome%in%open_canopy)] <- "open"

# table(wwf$canopy)
# closed   open 
# 7837   6514 

wwf <- sf::st_as_sf(wwf)
all_open <- subset(wwf, wwf_biomes$canopy=="open")
all_closed <- subset(wwf, wwf_biomes$canopy=="closed")

all_closed_simple <- ms_simplify(all_closed)
all_closed_simple <- st_union(all_closed_simple)

all_open_simple <- ms_simplify(all_open)
all_open_simple <- st_union(all_open_simple)

plot(all_open_simple, col="beige")
plot(all_closed_simple, col="lightgreen", add=T)

# tmp_map1 <- ggplot(data = all_open_simple) +
#   geom_sf(aes(fill = all_open_simple, lwd=0.5)) +
#   scale_fill_viridis_c(option = "plasma") +
#   theme_classic() 

#-----------------
# Points
thinned_points <- as.data.frame(fread("1_occurrence_data/cleaned_points.csv"))
coord <- thinned_points[,c("lon","lat")]
coordinates(coord) <- ~ lon + lat
plot(coord, col=adjustcolor("grey", alpha.f = 1), add=T, pch=16, cex=0.05)

