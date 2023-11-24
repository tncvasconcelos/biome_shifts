# rm(list=ls())
setwd("~/Desktop/biome_shifts")

library(ape)
library(phytools)

#--------------------------------
# Loading trait data 
all_areas <- list.files("3_organized_datasets_geohisse/")
all_habitats <- list()
for(i in 1:length(all_areas)){
  label <- gsub("_area_score.Rsave","",all_areas[i])
  load(paste0("3_organized_datasets_geohisse/", all_areas[i]))
  if(exists("habitats_subset")) {
    all_habitats[[i]] <- habitats_subset
    names(all_habitats)[i] <- label
    rm("habitats_subset")  
  }  
}

#--------------------------------
# Loading tree data 
all_tree_files <- list.files("4_organized_trees_geohisse/")
all_trees <- list()
for(i in 1:length(all_tree_files)){
  label <- gsub(".Rsave","",all_tree_files[i])
  load(paste0("4_organized_trees_geohisse/", all_tree_files[i]))
  if(exists("pruned_tree")) {
    all_trees[[i]] <- pruned_tree
    names(all_trees)[i] <- label
    rm("pruned_tree")  
  }  
}


#--------------------------------
# Plot to visualize habitat distribution in the tree:
pdf("inspecting_trees.pdf")

for(i in 1:length(all_trees)) {
  one_tree <- all_trees[[i]]
  one_label <- names(all_trees)[i]
  one_habitat <- all_habitats[[which(names(all_habitats)==one_label)]]
  
  one_habitat$area <- NA
  
  for (j in 1:length(one_habitat$area)){
    if (one_habitat[j, "area_open"] == 1 & one_habitat[j, "area_closed"] == 1){
      one_habitat[j, "area"] = "widespread"
    }
    if (one_habitat[j, "area_open"] == 0 & one_habitat[j, "area_closed"] == 1){
      one_habitat[j, "area"] = "closed"
    }
    if (one_habitat[j, "area_open"] == 1 & one_habitat[j, "area_closed"] == 0){
      one_habitat[j, "area"] = "open"
    }
  }
  
  one_tree$tip.label <- gsub(" ","_",one_tree$tip.label)
  intersected_species <- intersect(one_habitat[,"species"], one_tree$tip.label)
  one_tree <- keep.tip(one_tree, intersected_species)
  one_habitat <- subset(one_habitat, one_habitat$species%in%intersected_species)
  
  group_traits <- one_habitat[one_habitat[,"species"] %in% one_tree$tip.label,]
  group_traits <- group_traits[order(match(group_traits[,"species"],one_tree$tip.label)),]
  
  mode <- group_traits[,"area"]
  names(mode) <- group_traits[,1]
  colors_states <- c("midnightblue", "goldenrod","darkred")
  
  plot(one_tree, show.tip.label=T, edge.width=0.2, adj=1, cex=0.05)
  par(fg="transparent")
  tiplabels(pie=to.matrix(mode, sort(unique(mode))),piecol=colors_states,cex= 0.2,lwd=0.2, frame = "n")
  par(fg="black")
  #tiplabels(pch=1, bg=tip.cols, adj=1, cex=0.1, width = 0.1)
  legend("topleft", legend=sort(unique(mode)), pt.bg = colors_states, pch=21, cex=0.8)
  axisPhylo()
  title(one_label)
}
dev.off()






