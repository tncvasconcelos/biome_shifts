# rm(list=ls())
# setwd("biome_shifts")

#--------------------------------------
source("00_utility_functions.R")

all_trees <- load.trees("2_trees")

for(tree_index in 1:length(all_trees)){
  all_trees[[tree_index]]$tip.label <- resolve.names(all_trees[[tree_index]]$tip.label)
  one_tree <- all_trees[[tree_index]]
  save(one_tree, file=paste0("2_trees/",names(all_trees)[tree_index],".Rsave"))
}

