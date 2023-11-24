
# rm(list=ls())
setwd("~/Desktop/biome_shifts")

target_trees <- read.csv("all_trees_used.csv")
tree_labels <- target_trees$label

all_trees <- list.files("~/Desktop/treebank/treebank/4_cleaned_trees/")
all_trees <- all_trees[grep(paste(tree_labels, collapse="|"), all_trees)]

library(ape)

for(i in 1:length(all_trees)) {
    load(paste0("~/Desktop/treebank/treebank/4_cleaned_trees/",all_trees[i]))
    if(exists("one_tree")) {
      save(one_tree, file=paste0("2_trees/", all_trees[i]))
    }
  }



