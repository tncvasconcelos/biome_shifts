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

load(to_load[2])
failed_load <- to_load[sapply(to_load, quick_check)]
to_load <- to_load[!sapply(to_load, quick_check)]

individual_recon <- function(model_res){
  if(class(model_res)[1] == "try-error"){
    recon <- NA
  }else{
    hidden.states <- ncol(model_res$trans.matrix )/4
    recon <- MarginReconMuHiSSE(phy = model_res$phy, data = model_res$data, hidden.states = hidden.states, f = model_res$f, pars = model_res$solution, root.type = model_res$root.type, root.p = model_res$root.p, n.cores = 4, get.tips.only = FALSE)
  }
  return(recon)
}

convert_muhisse_pars <- function(model_res){
  est_pars <- model_res$solution[model_res$index.par < max(model_res$index.par)]
  trans_matrix <- model_res$trans.matrix
  hidden_states <- dim(trans_matrix)[1]/4
  trans_rates <- est_pars[grep("q", names(est_pars))]
  turnover <- eps <- setNames(rep(0, 4 * hidden_states), colnames(trans_matrix))
  turnover_tmp <- est_pars[grep("turnover", names(est_pars))]
  eps_tmp <- est_pars[grep("eps", names(est_pars))]
  if(hidden_states == 1){
    names(trans_rates) <- gsub("A", "", names(trans_rates))
    names(turnover_tmp) <- gsub("A", "", names(turnover_tmp))
    names(eps_tmp) <- gsub("A", "", names(eps_tmp))
    names(trans_rates) <- gsub("q", "", names(trans_rates))
    names(turnover_tmp) <- gsub("turnover", "", names(turnover_tmp))
    names(eps_tmp) <- gsub("eps", "", names(eps_tmp))
    from_to <- do.call(rbind, strsplit(names(trans_rates), "_"))
    for(i in 1:length(eps_tmp)){
      eps[grep(names(eps_tmp[i]), names(eps))] <- eps_tmp[i]
      turnover[grep(names(turnover_tmp[i]), names(eps))] <- turnover_tmp[i]
    }
  }else{
    names(trans_rates) <- gsub("q", "", names(trans_rates))
    names(turnover_tmp) <- gsub("turnover", "", names(turnover_tmp))
    names(eps_tmp) <- gsub("eps", "", names(eps_tmp))
    from_to <- do.call(rbind, strsplit(names(trans_rates), "_"))
    for(i in 1:length(eps_tmp)){
      eps[grep(names(eps_tmp[i]), names(eps))] <- eps_tmp[i]
      turnover[grep(names(turnover_tmp[i]), names(eps))] <- turnover_tmp[i]
    }
  }
  for(i in 1:dim(from_to)[1]){
    focal_trans <- from_to[i,]
    from_index <- grep(focal_trans[1], rownames(trans_matrix))
    to_index <- grep(focal_trans[2], rownames(trans_matrix))
    trans_matrix[from_index, to_index] <- trans_rates[i]
  }
  trans_matrix[trans_matrix == 0] <- NA
  return(list(turnover = turnover, eps = eps, trans_matrix = trans_matrix))
}


model_res <- res[[1]]

muhisse_pars <- convert_muhisse_pars(model_res)
max_t <- max(branching.times((model_res$phy)))
max_taxa <- Ntip(model_res$phy)
simulated.result <- SimulateHisse(muhisse_pars$turnover, muhisse_pars$eps, muhisse_pars$trans_matrix, max.t = max_t, max.taxa = max_taxa, x0=c(0.3, 0.5, 0, 0.2), nstart = 2)

# for(i in 1:length(to_load)){
#   load(to_load[i])
#   aic_weights <- GetAICWeights(res)
#   good_res <- res[aic_weights > 1e-2]
#   recon <- vector("list", length(res))
#   recon[aic_weights > 1e-2] <- mclapply(good_res, individual_recon, mc.cores=length(good_res))
#   file_name <- gsub("5_results/", "6_recons/", to_load[i])
#   file_name <- gsub("results_", "recon_", file_name)
#   save(recon, file=file_name)
# }



# par(mfcol=c(1,2))
# plot(SimToPhylo(simulated.result$results, include.extinct=TRUE))
# nodelabels()
# plot(SimToPhylo(simulated.result$results, include.extinct=FALSE))
# nodelabels()
