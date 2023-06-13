getTransMatFromOutput <- function(q_mat, focal_model){
  mismatch_index <- stringr::str_length(rownames(q_mat)[1]) < 3
  for(i in 1:dim(q_mat)[1]){
    focal_row <- rownames(q_mat)[i]
    for(j in 1:dim(q_mat)[2]){
      focal_column <- colnames(q_mat)[j]
      if(focal_row == focal_column){
        next
      }else{
        if(mismatch_index){
          focal_area <- substr(rownames(q_mat)[1], 1, 1)
          focal_trans <- paste0("d", focal_area, focal_row, "_", focal_area, focal_column)
          par <- focal_model$solution[names(focal_model$solution) == focal_trans]
          q_mat[i,j] <- par
        }else{
          focal_trans <- paste0("d", focal_row, "_", focal_column)
          par <- focal_model$solution[names(focal_model$solution) == focal_trans]
          q_mat[i,j] <- par
        }
      }
    }
  }
  return(q_mat)
}

# rm(list=ls())
setwd("~/biome_shifts/")

library(ape)
library(hisse)
library(parallel)
require(dplyr)

# finished model sets for particular datsets
clades <- dir("5_results/") %>% gsub("results_", "", .) %>% gsub(".RData", "", .) %>% unique(.)
to_load_results <- dir("5_results/", full.names = TRUE)
to_load_recons <- dir("6_recons/", full.names = TRUE)
to_load_model_tables <- dir("tables/", full.names = TRUE)
load("model_set.rsave")

model_table <- read.csv(to_load_model_tables[1])
load(to_load_results[1])
load(to_load_recons[1])

?GeoHiSSE
aicwt <- model_table$AICwt
aicwt[is.na(aicwt)] <- 0 # failed runs get 0

nboot = 100
models_to_simulate <- sample(1:length(aicwt), size = 1, replace = TRUE, prob = aicwt)

focal_model <- res[[models_to_simulate]]
focal_recon <- recon[[models_to_simulate]]
focal_index <- model_set[[models_to_simulate]]

add_jumps <- models_to_simulate > 18
add_extinction <- any(models_to_simulate == c(7:18, 24:36))
hidden_states <- dim(focal_index[[3]])[1]/3
root_probs <- focal_recon$node.mat[1, -1]
root_state <- sample(names(root_probs), 1, prob = root_probs)
x0 <- gsub("\\(", "", root_state) %>% gsub("\\)", "", .)
par.areas <- paste0(c("01", "0", "1"), rep(LETTERS[1:hidden_states], each = 3)) #from internal geohisse sim code

geohisse_pars <- SimulateGeoHiSSE(hidden.traits = hidden_states - 1, return.GeoHiSSE_pars = TRUE)

geohisse_pars$q.01 <- getTransMatFromOutput(geohisse_pars$q.01, focal_model)
geohisse_pars$q.0 <- getTransMatFromOutput(geohisse_pars$q.0, focal_model)
geohisse_pars$q.1 <- getTransMatFromOutput(geohisse_pars$q.1, focal_model)
geohisse_pars$model.pars
focal_model$
focal_model$solution[grep("d01A", names(focal_model$solution))]
focal_model$solution[grep("tau01A", names(focal_model$solution))]
names(focal_model$solution)

TranslateParsMakerGeoHiSSE(k = hidden_states-1, add.extinction = add_extinction, add.jumps = add_jumps)

sim_data <- SimulateGeoHiSSE(pars = geohisse_pars, 
                             hidden.traits = hidden_states-1,
                             x0 = x0,
                             max.taxa = length(focal_model$phy$tip.label),
                             add.jumps = add_jumps,
                             add.extinction = add_extinction)

# taxa stopping point following caetano et al 2018

?tree.geosse


