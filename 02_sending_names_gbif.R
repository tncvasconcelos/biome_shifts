
# rm(list=ls())
# setwd("~/Desktop/treebank/treebank")
setwd("~/Desktop/biome_shifts")
library(rgbif)
source("00_utility_functions.R")

trees <- load.trees("2_trees/")

all_tips_to_gbif <- c()
for(tree_index in 1:length(trees)) {
  tips_one_tree <- resolve.names(trees[[tree_index]]$tip.label)
  all_tips_to_gbif <- c(all_tips_to_gbif, tips_one_tree)
}

# Now we send a request to GBIF to download the points for this list of species 
user <- "" # username
pwd <- "" # password
email <- "@gmail.com" # email
rgbif::occ_download(rgbif::pred_in("scientificName", all_tips_to_gbif),
                    pred_in("basisOfRecord", 'PRESERVED_SPECIMEN'),
                    pred("hasCoordinate", TRUE),
                    format = "SIMPLE_CSV", user=user,pwd=pwd,email=email) # Sending request to GBIF  
