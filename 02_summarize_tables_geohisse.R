# rm(list=ls())
setwd("~/Desktop/biome_shifts/")
source("00_utility_functions.R")

habitats <- read.csv("summarized_habitat.csv")
all_trees <- load.trees("2_trees")

colnames(habitats)[1] <- "species"
habitats$species <- simplify.names.taxize(habitats$species)
ntips_total <- data.frame(tree=names(all_trees),ntips=NA)
for(i in 1:length(all_trees)) {
  all_trees[[i]]$tip.label <- simplify.names.taxize(all_trees[[i]]$tip.label)
  habitats_subset <- subset(habitats, habitats$species %in% all_trees[[i]]$tip.label)
  if(nrow(habitats_subset)>0) {
    pruned_tree <- keep.tip(all_trees[[i]], habitats_subset$species)
    habitats_subset$species <- gsub(" ","_",habitats_subset$species)
    write.csv(habitats_subset, file=paste0("3_organized_datasets_geohisse/",names(all_trees)[i],"_area_score.csv"), row.names=F)
    write.tree(pruned_tree, file=paste0("4_organized_trees_geohisse/",names(all_trees)[i],".tre"))
    ntips_total[i,2] <- Ntip(pruned_tree)    
  }
}

#ntips_total$tree[is.na(ntips_total$ntips)]
#one_subset<-subset(ntips_total, ntips_total$ntips>=100)
#sum(one_subset$ntips)

