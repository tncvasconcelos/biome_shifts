# rm(list=ls())
# setwd("~/Desktop/biome_shifts")

#--------------------------------------
library(rgbif)
source("00_utility_functions.R")

trees <- load.trees("2_trees/")

# sum(unlist(lapply(trees, Ntip)))
# [1] 11593

# Now we send a request to GBIF to download the points for this list of species 
user <- "" # username
pwd <- "" # password
email <- "@gmail.com" # email

# Get all GBIF synonyms
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
                      pred("hasCoordinate", TRUE),
                      format = "SIMPLE_CSV", user=user,pwd=pwd,email=email) # Sending request to GBIF  
  
  all_tips_to_gbif <- c(all_tips_to_gbif, tips_one_tree_syn)
  save(all_tips_to_gbif, file="1_occurrence_data/list_of_synonyms.Rsave")
}


# sum(2259058,
# 1080134,
# 30638,
# 2839886,
# 80223,
# 1986922,
# 11818494,
# 1533602,
# 9536019,
# 900717,
# 491703,
# 1285232,
# 7654275)

