# rm(list=ls())
setwd("~/Desktop/biome_shifts/")
source("00_utility_functions.R")

habitats <- read.csv("1_occurrence_data/summarized_habitat.csv")
all_trees <- load.trees("2.1_trees")

#which(names(all_trees)=="Croton-Arevalo_et_al-2017")

colnames(habitats)[1] <- "species"
#habitats$species <- simplify.names.taxize(habitats$species)
ntips_total <- data.frame(tree=names(all_trees),ntips=NA)
for(i in 1:length(all_trees)) {
  #all_trees[[i]]$tip.label <- simplify.names.taxize(all_trees[[i]]$tip.label)
  habitats_subset <- subset(habitats, habitats$species %in% all_trees[[i]]$tip.label)
  if(nrow(habitats_subset)>0) {
    pruned_tree <- keep.tip(all_trees[[i]], habitats_subset$species)
    habitats_subset$species <- gsub(" ","_",habitats_subset$species)
    
    splitted <- strsplit(habitats_subset$species, "_")
    flag_genus_only <- c()
    for(j in 1:length(splitted)) {
      if(length(is.na(lapply(splitted[[j]], "[", 2)))==1) {
        flag_genus_only <- c(flag_genus_only, j)
      }
      if(grepl("^[[:upper:]]", splitted[[j]][2])) {
        flag_genus_only <- c(flag_genus_only, j)
      }
    }
    if(length(flag_genus_only)>0) {
      habitats_subset <- habitats_subset[-flag_genus_only,] 
    }
    
    save(habitats_subset, file=paste0("3_organized_datasets_geohisse/",names(all_trees)[i],"_area_score.Rsave"), row.names=F)
    save(pruned_tree, file=paste0("4_organized_trees_geohisse/",names(all_trees)[i],".Rsave"))
    
    ntips_total[i,2] <- Ntip(pruned_tree)    
  }
}

#ntips_total$tree[is.na(ntips_total$ntips)]
#one_subset<-subset(ntips_total, ntips_total$ntips>=100)
sum(ntips_total$ntips)
# [1] 19424
