# rm(list=ls())
# setwd("~/Desktop/biome_shifts")

#--------------------------------------
library(rgbif)
source("00_utility_functions.R")

# Only trees that are reasonably sampled will be kept
all_trees <- load.trees("2_trees/")
summary_trees <- read.csv("treebank_info_simplified.csv")
summary_trees <- subset(summary_trees, summary_trees$large_and_well_sampled=="y")
trees <- subset(all_trees, names(all_trees) %in% summary_trees$label)

# Now we send a request to GBIF to download the points for this list of species 
user <- "" # username
pwd <- "" # password
email <- "thais.nogales@gmail.com" # email

# Get all synonyms
all_tips_to_gbif <- c()
for(tree_index in 1:length(trees)) {
  name <- trees[[tree_index]]$tip.label
  all_names <- list()
  for(i in 1:length(name)) {
    one_name <- name[i]
    taxon_info <- name_backbone(one_name)
    if("scientificName" %in% colnames(taxon_info)) {
      all_names_for_one_taxon <- taxon_info$scientificName
      synonyms <- name_usage(key = taxon_info$usageKey, data = "synonyms")$data
      if(nrow(synonyms)>0) {
        synonym_names <- synonyms$scientificName
        all_names_for_one_taxon <- c(all_names_for_one_taxon, unique(c(one_name, synonym_names)))
      }
    }
    all_names[[i]] <- all_names_for_one_taxon
    cat(i, "\r")
  }
  tips_one_tree_syn <- unlist(all_names)
  
  rgbif::occ_download(rgbif::pred_in("scientificName", tips_one_tree_syn),
                      #pred_in("basisOfRecord", 'PRESERVED_SPECIMEN'),
                      pred("hasCoordinate", TRUE),
                      format = "SIMPLE_CSV", user=user,pwd=pwd,email=email) # Sending request to GBIF  
  
  all_tips_to_gbif <- c(all_tips_to_gbif, tips_one_tree_syn)
  save(all_tips_to_gbif, file="1_occurrence_data/list_of_synonyms.Rsave")
}

#load("1_occurrence_data/list_of_synonyms.Rsave")

