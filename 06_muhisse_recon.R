# rm(list=ls())
setwd("~/biome_shifts/")

library(ape)
library(hisse)
library(parallel)

# finished model sets for particular datsets
to_load <- dir("5_results/", full.names = TRUE)

quick_check <- function(model_path){
  load(model_path)
  return(all(unlist(lapply(res, function(x) class(x) == "try-error"))))
}

failed_load <- to_load[sapply(to_load, quick_check)]
to_load <- to_load[!sapply(to_load, quick_check)]

load(to_load[2])
model_res <- res[[1]]

individual_recon <- function(model_res){
  if(class(model_res)[1] == "try-error"){
    recon <- NA
  }else{
    hidden.states <- ncol(model_res$trans.matrix )/4
    recon <- MarginReconMuHiSSE(phy = model_res$phy, data = model_res$data, hidden.states = hidden.states, f = model_res$f, pars = model_res$solution, root.type = model_res$root.type, root.p = model_res$root.p, n.cores = 4, get.tips.only = FALSE)
  }
  return(recon)
}

for(i in 1:length(to_load)){
  load(to_load[i])
  aic_weights <- GetAICWeights(res)
  good_res <- res[aic_weights > 1e-2]
  recon <- vector("list", length(res))
  recon[aic_weights > 1e-2] <- mclapply(good_res, individual_recon, mc.cores=length(good_res))
  file_name <- gsub("5_results/", "6_recons/", to_load[i])
  file_name <- gsub("results_", "recon_", file_name)
  save(recon, file=file_name)
}

