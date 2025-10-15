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

to_run <- c(
  to_load[grep("idx1.RData", to_load)][-c(1,3)],
  to_load[grep("idx4.RData", to_load)][-c(1,2,3)],
  to_load[grep("idx31.RData", to_load)[-c(1,3)]]
)

load(to_load[1])

individual_recon <- function(model_res){
  if(class(model_res)[1] == "try-error"){
    recon <- NA
  }else{
    hidden.states <- ncol(model_res$trans.matrix )/4
    recon <- MarginReconMuHiSSE(phy = model_res$phy, data = model_res$data, hidden.states = hidden.states, f = model_res$f, pars = model_res$solution, root.type = model_res$root.type, root.p = model_res$root.p, get.tips.only = FALSE)
  }
  return(recon)
}

run_recon_path <- function(to_load){
  load(to_load)
  file_name <- gsub("5_results/", "6_recons/", to_load)
  file_name <- gsub("res_state", "recon", file_name)
  recon <- try(individual_recon(res))
  save(recon, file=file_name)
  return(NULL)
}

num_cores <- 32
mclapply(to_run, run_recon_path, mc.cores = num_cores)


