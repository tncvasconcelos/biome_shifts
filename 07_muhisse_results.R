# rm(list=ls())
setwd("~/biome_shifts/")
source("00_utility_functions.R")

library(ape)
library(hisse)
library(parallel)
library(BSDA)

# finished model sets for particular datsets
to_load <- dir("5_results/", full.names = TRUE)

quick_check <- function(model_path){
  load(model_path)
  return(all(unlist(lapply(res, function(x) class(x) == "try-error"))))
}

failed_load <- to_load[sapply(to_load, quick_check)]
to_load <- to_load[!sapply(to_load, quick_check)]

load(to_load[1])
res

# the model ids that are the "main" structures are 1 (dull null), 4 (bisse), and 31 (hisse null) 
clade_i <- 1
model_table <- c()
all_res <- list()
for(i in 1:13){
  clade_i <- i
  clade_res <- list()
  load(to_load[grep("idx1.RData", to_load)][clade_i])
  clade_res[[1]] <- res
  load(to_load[grep("idx4.RData", to_load)][clade_i])
  clade_res[[2]] <- res
  # load(to_load[grep("idx15.RData", to_load)][clade_i])
  # clade_res[[3]] <- res
  load(to_load[grep("idx31.RData", to_load)][clade_i])
  clade_res[[3]] <- res
  model_table <- rbind(model_table, setNames(GetAICWeights(clade_res), c("M1", "M2","M3")))
  all_res[[i]] <- clade_res
}

output_table <- list()
slope_table <- c()
for(i in 1:13){
  # LRT
  dLik_31 <- all_res[[i]][[1]]$loglik - all_res[[i]][[3]]$loglik
  dLik_32 <- all_res[[i]][[2]]$loglik - all_res[[i]][[3]]$loglik
  p_31 <- pchisq(q = -2*dLik_31, df = 7, lower.tail = FALSE)
  p_32 <- pchisq(q = -2*dLik_32, df = 3, lower.tail = FALSE)
  
  if(p_31 > 0.05 | p_32 > 0.05){
    next
  }
  
  # Pars
  quickConvert <- function(turn_rate, eps_rate, index = 3){
    tmp <- c()
    for(i in 1:length(turn_rate)){
      tmp[i] <- convertBetweenPars(c(NA, NA, NA, turn_rate[i], eps_rate[i]))[index]
    }
    return(tmp)
  }
  M1_pars <- all_res[[i]][[1]]$solution[all_res[[i]][[1]]$index.par < 7]
  M2_pars <- all_res[[i]][[2]]$solution[all_res[[i]][[2]]$index.par < 11]
  M3_pars <- all_res[[i]][[3]]$solution[all_res[[i]][[3]]$index.par < 14]
  
  tA <- mean(M3_pars[1:3])
  tB <- mean(M3_pars[14:16])
  fA <- mean(M3_pars[4:6])
  fB <- mean(M3_pars[17:19])
  dA <- mean(quickConvert(M3_pars[1:3], M3_pars[4:6]))
  dB <- mean(quickConvert(M3_pars[14:16], M3_pars[17:19]))
  qA <- mean(M3_pars[7:10])
  qB <- mean(M3_pars[20:23])
  qA_ln <- exp(mean(log(M3_pars[7:10])))
  qB_ln <- exp(mean(log(M3_pars[20:23])))
  
  tm <- (qA - qB)/(tA - tB)
  dm <- (qA - qB)/(dA - dB)
  tm_ln <- (qA_ln - qB_ln)/(tA - tB)
  dm_ln <- (qA_ln - qB_ln)/(dA - dB)
  
  slope_table <- rbind(slope_table, (c(t = tm, d = dm)))
  clade_i <- gsub(".*res_state(.*)_idx.*", "\\1", to_load[i])
  
  data.frame(clade = clade_i, 
model = "M1", loglik = all_res[[i]][[1]]$loglik, AIC = all_res[[i]][[1]]$AICc)
}

SIGN.test(slope_table[,1], md = 0)
SIGN.test(slope_table[,2], md = 0)




output_table <- data.frame()
for(i in 1:13){
  clade_i <- gsub(".*res_state(.*)_idx.*", "\\1", to_load[i])
  
  M1_pars <- all_res[[i]][[1]]$solution[all_res[[i]][[1]]$index.par < 7]
  M2_pars <- all_res[[i]][[2]]$solution[all_res[[i]][[2]]$index.par < 11]
  M3_pars <- all_res[[i]][[3]]$solution[all_res[[i]][[3]]$index.par < 14]
  
  all_names <- names(M3_pars)
  
  pad_pars <- function(par_vec){
    res <- setNames(rep(NA, length(all_names)), all_names)
    res[names(par_vec)] <- par_vec
    return(res)
  }
  
  M1_full <- pad_pars(M1_pars)
  M2_full <- pad_pars(M2_pars)
  M3_full <- pad_pars(M3_pars)
  
  tab_M1 <- data.frame(
    clade = clade_i,
    model = "M1",
    loglik = all_res[[i]][[1]]$loglik,
    AICc = all_res[[i]][[1]]$AICc,
    t(M1_full),
    check.names = FALSE
  )
  
  tab_M2 <- data.frame(
    clade = clade_i,
    model = "M2",
    loglik = all_res[[i]][[2]]$loglik,
    AICc = all_res[[i]][[2]]$AICc,
    t(M2_full),
    check.names = FALSE
  )
  
  tab_M3 <- data.frame(
    clade = clade_i,
    model = "M3",
    loglik = all_res[[i]][[3]]$loglik,
    AICc = all_res[[i]][[3]]$AICc,
    t(M3_full),
    check.names = FALSE
  )
  
  output_table <- rbind(output_table, tab_M1, tab_M2, tab_M3)
}

head(output_table)
