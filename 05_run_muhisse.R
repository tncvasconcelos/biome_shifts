# rm(list=ls())
setwd("~/biome_shifts")

library(ape)
library(hisse)
library(parallel)
library(partitions)

# Load trees
all_trees_files <- list.files("4_organized_trees_geohisse/")
all_trees <- list()
for(i in 1:length(all_trees_files)){
  load(paste0("4_organized_trees_geohisse/",all_trees_files[i]))
  pruned_tree$tip.label <- gsub(" ","_",pruned_tree$tip.label)
  all_trees[[i]] <- pruned_tree
  names(all_trees)[i] <- gsub(".Rsave","",all_trees_files[i])
}

# Load datasets
all_area_files <- list.files("3_organized_datasets_geohisse/")
all_areas <- list()
for(i in 1:length(all_area_files)){
  load(paste0("3_organized_datasets_geohisse/",all_area_files[i]))
  all_areas[[i]] <-  habitats_subset
  names(all_areas)[i] <- gsub("_area_score.Rsave","",all_area_files[i])
}


focal_clades <- gsub(".Rsave","",all_trees_files)
### Import data
state_list <- list()

for (group_index in 1:length(focal_clades)) {
  group = names(all_areas)[group_index] # name the group
  tree <- all_trees[[group_index]] # load tree file 
  dist <- all_areas[[group_index]] # load distribution file 
  
  # Preparing data - areas have to be as 0 (11 - widespread), 
  # 1 (10, endemic of first area) 
  # and 2 (01, endemic of second area)
  dist <- as.data.frame(dist)
  dist$area <- 0
  #   areas <- as.data.frame(rep(1, nrow(dist)))
  #   dist <- cbind(dist, areas)
  #   colnames(dist)[7] <- "area"
  
  for (i in 1:length(dist$area)){
    if (dist[i, "area_open"] == 1 & dist[i, "area_closed"] == 1){
      dist[i, "area"] = 0 
    }
    if (dist[i, "area_open"] == 0 & dist[i, "area_closed"] == 1){
      dist[i, "area"] = 1
    }
    if (dist[i, "area_open"] == 1 & dist[i, "area_closed"] == 0){
      dist[i, "area"] = 2
    }
  }
  state_list[[group_index]]<-dist[,c("species", "area")]
  names(state_list)[group_index] <- group
  # in the dataframe and tree
  included_species <- intersect(state_list[[group_index]][,1], tree$tip.label)
  state_list[[group_index]] <- state_list[[group_index]][match(state_list[[group_index]][,1], included_species), ]
  all_trees[[group_index]] <- keep.tip(all_trees[[group_index]], included_species)
  all_areas[[group_index]] <- all_areas[[group_index]][match(all_areas[[group_index]][,1], included_species), ]
}

# model partitions
# ef and turn
# 111
# 123
# 111 111
# 111 222
# 123 123
# 123 456
# q_mat
# 1<>2<>3
# 1<>2<>3 1<>2<>3
# 1<>2<>3 4<>5<>6

# single rate class models
ef_1 <- list(c(1,1,0,1),
           c(1,2,0,3))
turn_1 <- list(c(1,1,0,1),
             c(1,2,0,3))

index_1 <- expand.grid(1:length(ef_1), 1:length(turn_1), 1, "1R")

# two rate class models
ef_2 <- list(c(1,1,0,1,1,1,0,1),
             c(1,2,0,3,1,2,0,3),
             c(1,1,0,1,2,2,0,2),
             c(1,2,0,3,4,5,0,6))
turn_2 <- list(c(1,1,0,1,1,1,0,1),
               c(1,2,0,3,1,2,0,3),
               c(1,1,0,1,2,2,0,2),
               c(1,2,0,3,4,5,0,6))

index_2 <- expand.grid(1:length(ef_2), 1:length(turn_2), 1:2, "2R")

# make the q_mats
# one rate class
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=0, make.null = TRUE)
trans_rate_1 <- ParDrop(trans.rate, c(2,5,6,8))
# two rate classes
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=1, make.null = FALSE)
trans_rate_2cd <- ParDrop(trans.rate, c(2,5,6,8,10,13,14,16))
trans.rate <- TransMatMakerMuHiSSE(hidden.traits=1, make.null = TRUE)
trans_rate_2cid <- ParDrop(trans.rate, c(2,5,6,8))

trans_rate_2 <- list(trans_rate_2cid, trans_rate_2cd)

index <- as.data.frame(rbind(index_1, index_2))
index_list <- split(index, seq(nrow(index)))


ref_table <- read.csv("all_trees_used.csv")
# recalculating sfs
updated_sfs <- c()
for(group_index in 1:length(all_trees)) {
  tree <- all_trees[[group_index]] # load tree file 
  one_sf <- ref_table$sf_ingroup[which(ref_table$label==names(all_trees)[group_index])]
  n_tips_full <- ref_table$n_tips_ingroup[which(ref_table$label==names(all_trees)[group_index])]
  n_tips_now <- Ntip(tree)
  updated_sfs[group_index] <- round((n_tips_now*one_sf)/n_tips_full, 2)
  names(updated_sfs)[group_index] <- names(all_trees)[group_index]
}

index_row <- index_list[[13]]
# two rate model
run_single_model <- function(dat, phy, sf, index_row, ef_1, ef_2, turn_1, turn_2, trans_rate_1, trans_rate_2){
  f = rep(sf, 4)
  if(index_row$Var4 == "1R"){
    # run single rate class model
    ef <- ef_1[[index_row$Var1]]
    turn <- turn_1[[index_row$Var2]]
    trans_rate <- trans_rate_1
    out <- MuHiSSE(phy=phy, data=dat, f=f, turnover=turn,
                   eps=ef, hidden.states=FALSE,
                   trans.rate=trans_rate)
  }else{
    # run double rate class model
    ef <- ef_2[[index_row$Var1]]
    turn <- turn_2[[index_row$Var2]]
    trans_rate <- trans_rate_2[[index_row$Var3]]
    out <- MuHiSSE(phy=phy, data=dat, f=f, turnover=turn,
                   eps=ef, hidden.states=TRUE,
                   trans.rate=trans_rate)
  }
  return(out)
}
# Preparing data - areas have to be as 0 (11 - widespread), 
# 1 (10, endemic of first area) 
# and 2 (01, endemic of second area
for(i in seq_len(length(state_list))){
  dat <- state_list[[i]]
  dat_muhisse <- cbind(dat, dat[,2])
  dat_muhisse[dat_muhisse[,2] == 0, 3] <- 1
  dat_muhisse[dat_muhisse[,2] == 1, c(2,3)] <- 0
  dat_muhisse[dat_muhisse[,2] == 2, c(2,3)] <- 1
  phy <- all_trees[[i]]
  sf <- updated_sfs[i]
  res <- mclapply(index_list, function(x) run_single_model(dat_muhisse, phy, sf, x, ef_1, ef_2, turn_1, turn_2, trans_rate_1, trans_rate_2), mc.cores=36)
  save(res, file=paste0("5_results/results_", names(all_trees)[i], ".RData"))
}

