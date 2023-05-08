# rm(list=ls())
setwd("~/biome_shifts/")

library(ape)

# Load trees
all_trees_files <- list.files("4_organized_trees_geohisse/")
all_trees <- lapply(paste0("4_organized_trees_geohisse/",all_trees_files), read.tree)
names(all_trees) <- gsub(".tre","",all_trees_files)

# Load datasets
all_area_files <- list.files("3_organized_datasets_geohisse/")
all_areas <- lapply(paste0("3_organized_datasets_geohisse/",all_area_files), read.csv)
names(all_areas) <- gsub("_area_score.csv","",all_area_files)

# I still have to double check a lot of these because it looks like we lost a lot of data
# when filtering the points. Let's start with the following trees for a test:

pilot <- c("Acacia-Renner_et_al-2019","Lamiales-Fonseca-2021","Poaceae-Spriggs_et_al-2014",
"Carex-MartinBravo_et_al-2019" ,"Onagraceae-Freyman_&_Hohna-2018","Quercus-Hipp_et_al-2017",
"Viburnum-Landis_et_al-2021")

pilot_areas <- all_areas[pilot]
pilot_trees <- all_trees[pilot]


#### GeoHiSSE code ####

# For more information, see:
# Caetano, D. S., O'Meara, B. C., & Beaulieu, J. M. (2018). 
# Hidden state models improve state‐dependent diversification approaches, 
# including biogeographical models. Evolution, 72(11), 2308-2324.

library(hisse)
library(parallel)

### Import data
pilot_states <- list()
for (group_index in 1:length(pilot_trees)) {
  group = names(pilot_areas)[group_index] # name the group
  tree <- pilot_trees[[group_index]] # load tree file 
  dist <- pilot_areas[[group_index]] # load distribution file 
  
  # Preparing data - areas have to be as 0 (11 - widespread), 
  # 1 (10, endemic of first area) 
  # and 2 (01, endemic of second area)
  
  areas <- as.data.frame(rep(1, nrow(dist)))
  dist <- cbind(dist, areas)
  colnames(dist)[7] <- "area"
  
  for (i in 1:length(dist$area)){
    if (dist[i, "area_open"] == 1 && dist[i, "area_closed"]  == 1){
      dist[i, "area"] = 0 
    }
    if (dist[i, "area_open"] == 0 && dist[i, "area_closed"]  == 1){
      dist[i, "area"] = 1
    }
    if (dist[i, "area_open"] == 1 && dist[i, "area_closed"]  == 0){
      dist[i, "area"] = 2
    }
  }
  pilot_states[[group_index]]<-dist[,c("species", "area")]
  names(pilot_states)[group_index] <- group
}

#table(states$area) # check if species-richness in each range make sense

# 2 - cr "endemic"
# 1 - non-cr "endemic"
# 0 - widespread

# Load sampling fraction for the group
sf<-c(1,1,1) # e.g. if it's fully sampled 

#

# For each run, we will use one pilot_tree and the equivalent pilot_states

phy=tree
dat=states

# We used the same 18 models of Caetano et al. (2018) (plus a second set of models including jump dispersal) - see their original publication for more information

###############################################################################
## Mods 1 through 6
## Block of GeoSSE-like models.
## Here extirpation is linked to range reduction.
###############################################################################

###############################################################################
## Mods 7 through 12
## Block of GeoSSE+extinction models.
## Here extirpation is NOT linked to range reduction.
## Range reduction is different from the extinction of an endemic lineage.
###############################################################################

###############################################################################
## Mods 13 through 18
## Block of anagenetic geographic models (MuSSE).
## Here models emulate GeoSSE (or GeoHiSSE) but changes only happen along branches.
###############################################################################

#################### vvvvv JUMP vvvvv MODELS vvvvv ########################## 
###############################################################################
## Mods 19 through 24
## Block of GeoSSE-like models.
## Here extirpation is linked to range reduction.
###############################################################################

###############################################################################
## Mods 25 through 30
## Block of GeoSSE+extinction models.
## Here extirpation is NOT linked to range reduction.
## Range reduction is different from the extinction of an endemic lineage.
###############################################################################

###############################################################################
## Mods 31 through 36
## Block of anagenetic geographic models (MuSSE).
## Here models emulate GeoSSE (or GeoHiSSE) but changes only happen along branches.
###############################################################################

model_set <- list(
    list(
        ## Model 1 - Dispersal parameters vary only, no range-dependent diversification. 
        speciation <- c(1,1,1),
        extirpation <- c(1,1),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, include.jumps=FALSE, separate.extirpation=FALSE) 
    ),
    list(
        ## Model 2. Canonical GeoSSE model, range effect on diversification 
        speciation <- c(1,2,3),
        extirpation <- c(1,2),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, include.jumps=FALSE, separate.extirpation=FALSE) 
    ),
    list(
        ## Model 3. Heterogeneous diversification, not tied to range evolution.
        speciation <- c(1,1,1,2,2,2,3,3,3),
        extirpation <- c(1,1,2,2,3,3),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=2, make.null=TRUE, include.jumps=FALSE, separate.extirpation=FALSE)
    ),
    list(
        ## Model 4. Heterogeneous diversification, tied to range evolution. 
        ## Assumes 6 distinct diversification rates.
        speciation <- c(1,2,3,4,5,6),
        extirpation <- c(1,2,3,4),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, include.jumps=FALSE , 
                                            separate.extirpation=FALSE)

    ),
    list(
        ## Model 5. Heterogeneous diversification, not tied to range evolution. 
        ## Assumes 5 distinct diversification rates.
        speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3)),
        extirpation <- c(rep(1,2), rep(2,2), rep(3,2), rep(4,2), rep(5,2)),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=4, make.null=TRUE, 
                                            include.jumps=FALSE, separate.extirpation=FALSE) 
    ),
    list(
        ## Model 6. Heterogeneous diversification, not tied to range evolution. 
        ## Assumes two distinct diversification rates.
        speciation <- c(1,1,1,2,2,2),
        extirpation <- c(1,1,2,2),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, include.jumps=FALSE , 
                                            separate.extirpation=FALSE)
    ),
    list(
        ## Model 7 - Dispersal parameters vary only, no range-dependent diversification. 
        speciation <- c(1,1,1),
        extirpation <- c(1,1),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, include.jumps=FALSE
                                            , separate.extirpation=TRUE) 
    ),
    list(
        ## Model 8. GeoSSE model, with range effect on diversification
        speciation <- c(1,2,3),
        extirpation <- c(1,2),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, include.jumps=FALSE
                                            , separate.extirpation=TRUE) 
    ),
    list(
        ## Model 9. Heterogeneous diversification, not tied to range evolution.
        ## Assumes three distinct diversification rates.
        speciation <- c(1,1,1,2,2,2,3,3,3),
        extirpation <- c(1,1,2,2,3,3),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=2, make.null=TRUE,include.jumps=FALSE,
                                            separate.extirpation=TRUE)
    ),
    list(
        ## Model 10. Heterogeneous diversification, tied to range evolution.
        ## Assumes 6 distinct diversification rates.
        speciation <- c(1,2,3,4,5,6),
        extirpation <- c(1,2,3,4),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, include.jumps=FALSE, 
                                            separate.extirpation=TRUE) 
    ),
    list(
        ## Model 11. Heterogeneous diversification, not tied to range evolution. 
        ## Assumes 5 distinct diversification rates.
        speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3)),
        extirpation <- c(rep(1,2), rep(2,2), rep(3,2), rep(4,2), rep(5,2)),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=4, make.null=TRUE, include.jumps=FALSE,
                                            separate.extirpation=TRUE) 
    ),
    list(
        ## Model 12. Heterogeneous diversification, not tied to range evolution. 
        ## Assumes two distinct diversification rates.
        speciation <- c(1,1,1,2,2,2),
        extirpation <- c(1,1,2,2),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, include.jumps=FALSE , 
                                            separate.extirpation=TRUE)
    ),
    list(
        ## Model 13. Transitions only. No character effect on diversification
        speciation <- c(1,1,1),
        extirpation <- c(1,1,1),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, include.jumps=FALSE
                                            , separate.extirpation=TRUE) 
    ),
    list(
        ## Model 14. Character effect on diversification.
        speciation <- c(1,2,3),
        extirpation <- c(1,2,3),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, include.jumps=FALSE
                                            , separate.extirpation=TRUE) 
    ),
    list(
        ## Model 15. No character effect on diversification.
        speciation <- c(1,1,1,2,2,2,3,3,3),
        extirpation <- c(1,1,1,2,2,2,3,3,3),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=2, include.jumps=FALSE
                                            , separate.extirpation=TRUE, make.null=TRUE) 
    ),
    list(
        ## Model 16. Character effect on diversification, with a hidden state
        speciation <- c(1,2,3,4,5,6),
        extirpation <- c(1,2,3,4,5,6),
        trans.rate <- trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, include.jumps=FALSE
                                                        , separate.extirpation=TRUE)
    ),
    list(
        ## Model 17. No character effect on diversification, multiple shifts
        speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3)),
        extirpation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3)),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=4, include.jumps=FALSE,
                                            separate.extirpation=TRUE, make.null=TRUE) 
    ),
    list(
        ## Model 18. No character effect on diversification, multiple shifts.
        speciation <- c(rep(1,3), rep(2,3)),
        extirpation <- c(rep(1,3), rep(2,3)),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, include.jumps=FALSE, 
                                            separate.extirpation=TRUE, make.null=TRUE)
    ),
    list(
        ## Model 19 - Dispersal parameters vary only, no range-dependent diversification. 
        speciation <- c(1,1,1),
        extirpation <- c(1,1),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, include.jumps=TRUE 
                                            , separate.extirpation=FALSE) 
    ),
    list(
        ## Model 20. Canonical GeoSSE model, range effect on diversification 
        speciation <- c(1,2,3),
        extirpation <- c(1,2),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, include.jumps=TRUE
                                            , separate.extirpation=FALSE) 
    ),
    list(
        ## Model 21. Heterogeneous diversification, not tied to range evolution.
        ## Assumes three distinct diversification rates.
        ## Dispersion parameters across hidden areas are the same.
        speciation <- c(1,1,1,2,2,2,3,3,3),
        extirpation <- c(1,1,2,2,3,3),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=2, make.null=TRUE
                                            , include.jumps=TRUE, separate.extirpation=FALSE) 
    ),
    list(
        ## Model 22. Heterogeneous diversification, tied to range evolution. 
        ## Assumes 6 distinct diversification rates.
        speciation <- c(1,2,3,4,5,6),
        extirpation <- c(1,2,3,4),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, include.jumps=TRUE , 
                                            separate.extirpation=FALSE)
    ),
    list(
        ## Model 23. Heterogeneous diversification, not tied to range evolution. 
        ## Assumes 5 distinct diversification rates.
        speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3)) ,
        extirpation <- c(rep(1,2), rep(2,2), rep(3,2), rep(4,2), rep(5,2)) ,
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=4, make.null=TRUE, 
                                            include.jumps=TRUE, separate.extirpation=FALSE) 
    ),
    list(
        ## Model 24. Heterogeneous diversification, not tied to range evolution. 
        ## Assumes two distinct diversification rates.
        speciation <- c(1,1,1,2,2,2),
        extirpation <- c(1,1,2,2),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, include.jumps=TRUE , 
                                            separate.extirpation=FALSE)
    ),
    list(
        ## Model 25 - Dispersal parameters vary only, no range-dependent diversification. 
        speciation <- c(1,1,1),
        extirpation <- c(1,1),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, include.jumps=TRUE
                                            , separate.extirpation=TRUE) 
    ),
    list(
        ## Model 26. GeoSSE model, with range effect on diversification 
        speciation <- c(1,2,3),
        extirpation <- c(1,2),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, include.jumps=TRUE
                                            , separate.extirpation=TRUE) 
    ),
    list(
        ## Model 27. Heterogeneous diversification, not tied to range evolution.
        ## Assumes three distinct diversification rates.
        speciation <- c(1,1,1,2,2,2,3,3,3),
        extirpation <- c(1,1,2,2,3,3),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=2, make.null=TRUE,include.jumps=TRUE,
                                            separate.extirpation=TRUE)
    ),
    list(
        ## Model 28. Heterogeneous diversification, tied to range evolution.
        ## Assumes 6 distinct diversification rates.
        speciation <- c(1,2,3,4,5,6),
        extirpation <- c(1,2,3,4),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, include.jumps=TRUE, 
                                            separate.extirpation=TRUE) 
    ),
    list(
        ## Model 29. Heterogeneous diversification, not tied to range evolution. 
        ## Assumes 5 distinct diversification rates.
        speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3)) ,
        extirpation <- c(rep(1,2), rep(2,2), rep(3,2), rep(4,2), rep(5,2)) ,
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=4, make.null=TRUE, include.jumps=TRUE,
                                            separate.extirpation=TRUE) 
    ),
    list(
        ## Model 30. Heterogeneous diversification, not tied to range evolution. 
        ## Assumes two distinct diversification rates.
        speciation <- c(1,1,1,2,2,2),
        extirpation <- c(1,1,2,2),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, include.jumps=TRUE , 
                                            separate.extirpation=TRUE)
    ),
    list(
        ## Model 31. Transitions only. No character effect on diversification
        speciation <- c(1,1,1),
        extirpation <- c(1,1,1),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, include.jumps=TRUE
                                            , separate.extirpation=TRUE) 
    ),
    list(
        ## Model 32. Character effect on diversification.
        speciation <- c(1,2,3),
        extirpation <- c(1,2,3),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, include.jumps=TRUE
                                            , separate.extirpation=TRUE) 
    ),
    list(
        ## Model 33. No character effect on diversification.
        speciation <- c(1,1,1,2,2,2,3,3,3),
        extirpation <- c(1,1,1,2,2,2,3,3,3),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=2, include.jumps=TRUE
                                            , separate.extirpation=TRUE, make.null=TRUE) 
    ),
    list(
        ## Model 34. Character effect on diversification, with a hidden state
        speciation <- c(1,2,3,4,5,6),
        extirpation <- c(1,2,3,4,5,6),
        trans.rate <- trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, include.jumps=TRUE
                                                        , separate.extirpation=TRUE)
    ),
    list(
        ## Model 35. No character effect on diversification, multiple shifts
        speciation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3)),
        extirpation <- c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3)),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=4, include.jumps=TRUE,
                                            separate.extirpation=TRUE, make.null=TRUE) 
    ),
    list(
        ## Model 36. No character effect on diversification, multiple shifts.
        speciation <- c(rep(1,3), rep(2,3)),
        extirpation <- c(rep(1,3), rep(2,3)),
        trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, include.jumps=TRUE, 
                                            separate.extirpation=TRUE, make.null=TRUE)
    )
)

quickFunc <- function(par_list, dat, phy){
  hidden <- ifelse(dim(par_list[[3]])[1]>3, TRUE, FALSE) 
  res <- GeoHiSSE(phy, dat, f=sf, turnover=par_list[[1]], eps=par_list[[2]], hidden.states=hidden, trans.rate=par_list[[3]], assume.cladogenetic=FALSE)
  return(res)
}

for(i in seq_len(length(pilot_states))){
  dat <- pilot_states[[i]]
  phy <- pilot_trees[[i]]
  res <- mclapply(model_set, function(x) quickFunc(x, dat, phy), mc.cores=36)
  save(res, file=paste0("pilot_results/pilot_results_", names(pilot_trees)[i], ".RData"))
}

