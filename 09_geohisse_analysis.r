# rm(list=ls())
setwd("~/biome_shifts/")

library(ape)
library(hisse)
library(parallel)
require(dplyr)

# finished model sets for particular datsets
clades <- dir("5_results/") %>% gsub("results_", "", .) %>% gsub(".RData", "", .) %>% unique(.)
to_load_recon <- dir("6_recons/", full.names = TRUE)
to_load_results <- dir("5_results/", full.names = TRUE)

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
for(i in 1:length(to_load_results)){
  load(to_load_results[i])
  tables[[i]] <- getModelTable(res)
  write.csv(tables[[i]], file=paste0("tables/", clades[i], "_model_table.csv"), row.names=FALSE)
}
names(tables) <- clades

# a function for summarizing muhisse output
convertGeohisseFit2Pars <- function(res){
  tmp <- res
  pars <- tmp$solution
  data.castor <- tmp$data[,2]
  names(data.castor) <- tmp$data[,1]
  data.castor <- data.castor[tmp$phy$tip.label]
  
  NPstates <- length(unique(data.castor))
  Nstates <- NPstates * dim(tmp$trans.matrix)[1]/3
  
  proxy_map <- rep(1:NPstates, length.out = Nstates)
  
  # from the cache we get ALL parameters of interest
  transition_matrix <- pars[grep("d", names(pars))]
  birth_rates <- pars[grep("tau", names(pars))]
  death_rates <- pars[grep("ef", names(pars))]
  
  # get the s we are interested
  transition_matrix <- transition_matrix[transition_matrix > 0]
  birth_rates <- birth_rates[birth_rates > 0]
  death_rates <- death_rates[death_rates > 0]
  
  # get the indexes of all possible transition rates and compare it to the requested trans.rate
  from <- strsplit(names(transition_matrix), "_")
  from <- lapply(from, function(x) gsub("d", "", x = x[1]))
  from <- do.call(rbind, from)
  
  to <- strsplit(names(transition_matrix), "_") 
  to <- lapply(to, function(x) x[2]) 
  to <- do.call(rbind, to)
  
  transition_IndexMatrix <- paste(from, to, sep = "")
  trans.rate_IndexMatrix <- which(!(tmp$trans.matrix == 0 | is.na(tmp$trans.matrix)), arr.ind = TRUE)
  if(NPstates == Nstates){
    transition_IndexMatrix <- gsub(c("A|B|C|D|E|F|G|H|I|J"), "", transition_IndexMatrix)
  }
  rows2Match <- rownames(tmp$trans.matrix)[trans.rate_IndexMatrix[,1]]
  rows2Match <- gsub("\\(", "", rows2Match)
  rows2Match <- gsub("\\)", "", rows2Match)
  
  cols2Match <- colnames(tmp$trans.matrix)[trans.rate_IndexMatrix[,2]] 
  cols2Match <- gsub("\\(", "", cols2Match) 
  cols2Match <- gsub("\\)", "", cols2Match)
  
  trans.rate2Match <- paste(rows2Match, cols2Match, sep = "")
  neededTransitions <- transition_matrix[match(trans.rate2Match, transition_IndexMatrix)]
  full_mat <- matrix(0, dim(tmp$trans.matrix)[1], dim(tmp$trans.matrix)[2])
  
  for(i in 1:length(neededTransitions)){
    full_mat[trans.rate_IndexMatrix[i, 1], trans.rate_IndexMatrix[i, 2]] <- neededTransitions[i]
  }
  mat_reducer <- proxy_map + (rep(0:((Nstates/NPstates)-1), each = NPstates) * 3)
  full_mat <- full_mat[c(mat_reducer), c(mat_reducer)]
  colnames(full_mat) <- rownames(full_mat) <- colnames(tmp$trans.matrix)[mat_reducer]
  pars.matrix <- tmp$trans.matrix[c(mat_reducer), c(mat_reducer)]
  
  return(list(
    nHiddenStates = Nstates/NPstates,
    nObservedStates  = NPstates,
    parameter_matrix = pars.matrix, 
    turnover = birth_rates,
    extinction_fraction = death_rates,
    transition_matrix = full_mat
  ))
}

# a function for getting model averaged parameters
get_full_clade_result <- function(clade, model_files, recon_files){
  load(model_files[grep(clade, model_files)])
  model_table <- getModelTable(res)
  load(recon_files[grep(clade, recon_files)])
  return(list(res = res, recon = recon, model_table = model_table))
}

# get the full results for each clade
getTipPars <- function(result_i, recon_i){
  summ_results <- convertGeohisseFit2Pars(result_i)
  turn <- summ_results$turnover
  rate <- rowSums(summ_results$transition_matrix)
  n_max <- summ_results$nHiddenStates * summ_results$nObservedStates
  recon_mat <- recon_i$tip.mat[,c(2:(n_max+1))]
  tip_turn <- colSums(turn * t(recon_mat))
  tip_rate <- colSums(rate * t(recon_mat))
  return(list(tip_turn = tip_turn, tip_rate = tip_rate))
}

# full_result <- get_full_clade_result(clades[1], to_load_results, to_load_recon)
get_model_avg_params <- function(full_result, return_summ = FALSE){
  # remove try errors
  try_error_index <- is.na(full_result$model_table[,1])
  result <- full_result$res[!try_error_index]
  model_table <- full_result$model_table[!try_error_index,]
  recon <- full_result$recon[!try_error_index]
  if(return_summ){
    return(lapply(result, convertGeohisseFit2Pars))
  }
  all_tip_rates <- mapply(getTipPars, result, recon)
  tip_avg_turn <- tip_avg_rate <- 0
  for(i in 1:ncol(all_tip_rates)){
    tip_avg_turn <-  tip_avg_turn + all_tip_rates[,i][[1]] * model_table$AICwt[i]
    tip_avg_rate <- tip_avg_rate + all_tip_rates[,i][[2]] * model_table$AICwt[i]
  }
  return(data.frame(tip_avg_turn = tip_avg_turn, tip_avg_rate = tip_avg_rate))
}

# Preparing data - areas have to be as 0 (11 - widespread), 
# 1 (10, endemic of first area) 
# and 2 (01, endemic of second area)

# png("figures/pilot_turnover_rate.png", width = 1000, height = 1000)
# par(mfrow=c(7,1))
for(i in 1:length(clades)){
  full_result <- get_full_clade_result(clades[i], to_load_results, to_load_recon)
  if(all(is.na(full_result$model_table))){
    next
  }
  dat <- full_result$res[[1]]$data[,2]+1
  phy <- full_result$res[[1]]$phy
  tip_correlations <- get_model_avg_params(full_result, return_summ = FALSE)
  rownames(tip_correlations) <- full_result$res[[1]]$data[,1]
  # hist(tip_correlations$tip_avg_turn, xlab = "Average turnover", main = clades[i], col = dat)
  # hist(tip_correlations$tip_avg_rate, xlab = "Average rate", main = clades[i], col = dat)
  phy$root.edge <- NULL
  res_plm <- phylolm::phylolm(tip_avg_turn ~ tip_avg_rate, data = tip_correlations, phy = phy)
  sig <- summary(res_plm)
  to_add <- ifelse(sig$coefficients[2,4] < 0.05, "*", "")
  png(paste0("plots/plot_", clades[i], ".png"))
  plot(x = tip_correlations$tip_avg_rate, y = tip_correlations$tip_avg_turn, xlab = "Average rate", ylab = "Average turnover", main = paste0(clades[i], " ", to_add), pch = 16, col = dat)
  # res_lm <- lm(tip_correlations$tip_avg_turn ~ tip_correlations$tip_avg_rate)
  abline(res_plm, col = "red")
  dev.off()
}
# dev.off()
