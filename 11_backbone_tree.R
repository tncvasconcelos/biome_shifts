setwd("~/Desktop/biome_shifts/")
library(ape)
library(phytools)
library(phangorn)

load.trees <- function(tree.dir) {
  tree_files <- list.files(tree.dir, full.names = T)
  all_trees <- list()
  for(i in 1:length(tree_files)) {
    load(tree_files[i])
    if(exists("one_tree")) {
      all_trees[[i]] <- one_tree
      names(all_trees)[i] <- gsub(paste0(c(paste0(tree.dir,"/"), ".Rsave"), collapse="|"),"", tree_files[i])
      rm("one_tree")
    }
  }
  return(all_trees)
}

fix.names.taxize <- function(focal_species_trees) {
  for(name_index in 1:length(focal_species_trees)){
    one_tmp_string <- focal_species_trees[name_index]
    if(any(grepl("[()]", one_tmp_string))){
      splitted_names <- strsplit(one_tmp_string," ")[[1]]
      begin_author <- which(grepl("[()]", splitted_names))[1]
      species_name <- paste0(splitted_names[1:(begin_author-1)], collapse=" ")
      author <- splitted_names[begin_author:length(splitted_names)]
      old_authors <- author[grep("[()]", author)]
      end_first_half <- floor(length(old_authors)/2)
      before <- old_authors[1:end_first_half]
      after <- old_authors[(end_first_half+1):(length(old_authors))]
      if(paste(before,collapse = " ") == paste(after, collapse=" ")) {
        author <- paste(author[1:(length(author)/2)], collapse=" ")
        focal_species_trees[name_index] <- paste0(species_name, " ", author, collapse=" ")
      } else {
        author <- paste(author, collapse=" ")
        focal_species_trees[name_index] <- paste0(species_name, " ", author, collapse=" ")
      }
    }
  }
  return(focal_species_trees)
}

simplify.names.taxize <- function(names) {
  results <- c()
  for(name_index in 1:length(names)){
    one_tmp_string <- names[name_index]
    if(grepl("UNMATCHED",one_tmp_string)) {
      next 
    } else {
      splitted_names <- strsplit(one_tmp_string," ")[[1]]
      genus <- splitted_names[1]
      epiphet <- splitted_names[2]
      if(is.na(epiphet)) {
        next
      } else {
        if(any(grepl("indet_sp",splitted_names))) {
          full_name <- "tip_to_drop" # indet species
        } else {
          if(stringr::str_detect(epiphet,"[[:upper:]]")) {
            full_name <- "tip_to_drop" # indet species
          } else {
            if(length(splitted_names) > 2) {
              complement <- splitted_names[3:length(splitted_names)]
              if(grepl("[()]", complement[1])) {
                full_name <- paste(c(genus, epiphet), collapse = " ")
              } else {
                if(stringr::str_detect(complement[1],"[[:upper:]]")) {
                  full_name <- paste(c(genus, epiphet), collapse = " ")
                } else {
                  complement <- subset(complement, !stringr::str_detect(complement,"[[:upper:]]"))
                  complement <- subset(complement, !grepl(paste(c("[()]","&","([0-9]+).*$","^ex$"), collapse="|"), complement))
                  if(length(complement)==0){
                    full_name <- paste(c(genus, epiphet), collapse = " ")
                  } else {
                    full_name <- paste(c(genus, epiphet, complement), collapse = " ")
                  }
                }
              } 
            }
          }
        }  
      }
    }
    results[name_index] <- full_name
  }
  return(results)
}

get.node.age <- function (phy) {
  root.node <- length(phy$tip.label)+1
  seq.nodes <- phy$edge
  dists <- phy$edge.length
  res <- numeric(max(phy$edge))
  for (i in seq_len(nrow(seq.nodes))) {
    res[seq.nodes[i, 2]] <- res[seq.nodes[i,1]] + dists[i]
  }
  ages <- abs(round(res,3)-round(max(res),3))
  return(ages)
} 


# most recent common ancestor of a bunch of tips
big_tree <- readRDS("taxized_GBMB.Rdata")
big_tree$tip.label <- unname(big_tree$tip.label)
big_tree$tip.label <- simplify.names.taxize(big_tree$tip.label)

clade_trees_files <- list.files("2.1_trees", full.names=T)
labels <- gsub(paste0(c("2.1_trees/",".Rsave"), collapse="|"),"", clade_trees_files)

all_trees <- list()
for(i in 1:length(clade_trees_files)){
  load(clade_trees_files[i])
  all_trees[[i]] <- one_tree$tip.label
  names(all_trees)[i] <- labels[i]
} 

tips_to_keep <- c()
for(i in 1:length(all_trees)){
  one_label <- names(all_trees)[i]
  one_tree_tip_labels <- all_trees[[i]]
  one_tree_tip_labels <- simplify.names.taxize(one_tree_tip_labels)
  tips_to_keep_tmp <- which(big_tree$tip.label %in% one_tree_tip_labels)
  big_tree$tip.label[tips_to_keep_tmp] <- one_label
  tips_to_keep <- c(tips_to_keep, tips_to_keep_tmp)
}

pruned_big_tree <- keep.tip(big_tree, tips_to_keep)
pruned_big_tree$node.label <- NULL
backbone_tree <- keep.tip(pruned_big_tree, which(!duplicated(pruned_big_tree$tip.label)))


plot(ladderize(backbone_tree), cex=0.6)
axisPhylo()
write.tree(backbone_tree, file="backbone_tree.tre")
