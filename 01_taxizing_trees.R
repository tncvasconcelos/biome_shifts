# rm(list=ls())
# setwd("~/Desktop/biome_shifts")

#--------------------------------------
source("00_utility_functions.R")

all_trees <- load.trees("2_trees")

for(tree_index in 1:length(all_trees)){
  # 
  # all_trees[[tree_index]]$tip.label <- resolve.names(all_trees[[tree_index]]$tip.label)
  # one_tree <- all_trees[[tree_index]]
  # save(one_tree, file=paste0("2_trees/",names(all_trees)[tree_index],".Rsave"))
  # 
  one_tree <- all_trees[[tree_index]]
  names_to_solve <- one_tree$tip.label
  sources <- taxize::gna_data_sources()
  wcvp_name <- taxize::gna_verifier(names=names_to_solve,data_sources=sources$id[sources$title == "The World Checklist of Vascular Plants"],best_match_only=TRUE)$currentName
  one_tree$tip.label <- wcvp_name
  save(one_tree, file=paste0("2_trees/",names(all_trees)[tree_index],".Rsave"))
  
}




