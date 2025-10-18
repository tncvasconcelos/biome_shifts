# rm(list=ls())
# setwd("~/Desktop/biome_shifts")

#--------------------------------------
source("00_utility_functions.R")

all_trees <- load.trees("2_trees")

# Harmonizing tip taxonomy to follow accepted names from The World Checklist of Vascular Plants
for(tree_index in 1:length(all_trees)){
  one_tree <- all_trees[[tree_index]]
  names_to_solve <- one_tree$tip.label
  sources <- taxize::gna_data_sources()
  wcvp_name <- taxize::gna_verifier(names=names_to_solve,data_sources=sources$id[sources$title == "The World Checklist of Vascular Plants"],best_match_only=TRUE)$currentName
  one_tree$tip.label <- wcvp_name
  save(one_tree, file=paste0("2_trees/",names(all_trees)[tree_index],".Rsave"))
}

# get.node.age <- function (phy) {
#   root.node <- length(phy$tip.label)+1
#   seq.nodes <- phy$edge
#   dists <- phy$edge.length
#   res <- numeric(max(phy$edge))
#   for (i in seq_len(nrow(seq.nodes))) {
#     res[seq.nodes[i, 2]] <- res[seq.nodes[i,1]] + dists[i]
#   }
#   ages <- abs(round(res,3)-round(max(res),3))
#   return(ages)
# }
# 
# # Tree ages
# for(i in 1:length(all_trees)) {
#     print(names(all_trees)[i])
#     print(max(get.node.age(all_trees[[i]])))
#     print(Ntip(all_trees[[i]]))
# }
