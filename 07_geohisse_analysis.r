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


# 2 - cr "endemic"
# 1 - non-cr "endemic"
# 0 - widespread

# Load sampling fraction for the group
sf<-c(1,1,1) # e.g. if it's fully sampled 

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

# finished model sets for particular datsets
to_load <- dir("pilot_results/", full.names = TRUE)

# function for pulling out AIC and stuff like that
getModelRes <- function(model_res){
    if(class(model_res)[1] == "try-error"){
        out <- c(lnLik=NA, k=NA, AIC=NA)
    }else{
        out <- c(lnLik=model_res$loglik, k=max(model_res$index.par, na.rm=TRUE), AIC=model_res$AIC)
    }
    return(out)
}

getModelTable <- function(model_list){
    model_table <- as.data.frame(do.call(rbind, lapply(model_list, getModelRes)))
    model_table$dAIC <- model_table$AIC - min(model_table$AIC, na.rm=TRUE)
    model_table$AICwt <- exp(-0.5 * model_table$dAIC)/sum(exp(-0.5 * model_table$dAIC), na.rm=TRUE)
    return(model_table)
}

tables <- list()
for(i in 1:length(to_load)){
    load(to_load[i])
    tables[[i]] <- getModelTable(res)
}
names(tables) <- dir("pilot_results/")


# pilot_results_Acacia-Renner_et_al-2019.RData
        # AICwt = 80%
        ## Model 11. Heterogeneous diversification, not tied to range evolution. No jump dispersal.
        ## Assumes 5 distinct diversification rates.

# pilot_results_Carex-MartinBravo_et_al-2019.RData
        # AICwt = 60%
        ## Model 22. Heterogeneous diversification, tied to range evolution. Jump dispersal.
        ## Assumes 6 distinct diversification rates.
        # AICwt = 28% 
        ## Model 10. Heterogeneous diversification, tied to range evolution.
        ## Assumes 6 distinct diversification rates.

# pilot_results_Lamiales-Fonseca-2021.RData
        # AICwt = 100%
        ## Model 28. Heterogeneous diversification, tied to range evolution. Jump dispersal.
        ## Assumes 6 distinct diversification rates.

# pilot_results_Onagraceae-Freyman_&_Hohna-2018.RData
        # AICwt = 88%
        ## Model 12. Heterogeneous diversification, not tied to range evolution. 
        ## Assumes two distinct diversification rates.

# pilot_results_Poaceae-Spriggs_et_al-2014.RData
        # AICwt = 99%
        ## Model 10. Heterogeneous diversification, tied to range evolution.
        ## Assumes 6 distinct diversification rates.

# pilot_results_Quercus-Hipp_et_al-2017.RData
        # AICwt = 99%
        ## Model 4. Heterogeneous diversification, tied to range evolution. 
        ## Assumes 6 distinct diversification rates.

# pilot_results_Viburnum-Landis_et_al-2021.RData
        # AICwt = 32%
        ## Model 9. Heterogeneous diversification, not tied to range evolution.
        ## Assumes three distinct diversification rates.
        # AICwt = 25%
        ## Model 20. Canonical GeoSSE model, range effect on diversification. Jump dispersal.
        # AICwt = 12%
        ## Model 27. Heterogeneous diversification, not tied to range evolution.
        ## Assumes three distinct diversification rates.