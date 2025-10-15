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

nm_idx1 <- to_load[grep("idx1.RData", to_load)]
all_res <- vector("list", length(nm_idx1))
model_table_list <- vector("list", length(nm_idx1))

for(i in seq_along(nm_idx1)){
  clade_name <- sub("_idx1\\.RData$","", basename(nm_idx1[i]))
  files <- to_load[grepl(paste0("^", clade_name, "_idx[0-9]+\\.RData$"), basename(to_load))]
  if(length(files) == 0) next
  idxs <- as.integer(sub(".*_idx([0-9]+)\\.RData$","\\1", files))
  ord <- order(idxs)
  files <- files[ord]; idxs <- idxs[ord]
  clade_res <- lapply(files, function(f){ load(f); res })
  names(clade_res) <- paste0("M", idxs)
  all_res[[i]] <- clade_res
  model_table_list[[i]] <- setNames(GetAICWeights(clade_res), names(clade_res))
}

# combine AIC weight vectors into a rectangular data.frame (NA where a model is missing)
all_models <- sort(unique(unlist(lapply(model_table_list, names))))
model_table <- t(sapply(model_table_list, function(x){
  v <- setNames(rep(NA_real_, length(all_models)), all_models)
  if(!is.null(x)) v[names(x)] <- x
  v
}))
rownames(model_table) <- sub("_idx1\\.RData$","", basename(nm_idx1))
model_table <- as.data.frame(model_table)
model_table <- model_table[,paste0("M", 1:36)]


quickConvert2 <- function(turn_rate, eps_rate, index = 3){
  tmp <- c()
  for(i in 1:length(turn_rate)){
    tmp[i] <- convertBetweenPars(c(NA, NA, NA, turn_rate[i], eps_rate[i]))[index]
  }
  return(tmp)
}

slope_table <- c()
summary_table <- data.frame()
nm_idx1 <- to_load[grep("idx1.RData", to_load)]
for(i in 1:13){
  best_model <- names(which.max(model_table[i,]))
  # Pars
  focal_model <- all_res[[i]][[best_model]]
  M_pars <- focal_model$solution[focal_model$index.par < max(focal_model$index.par)]
  if(dim(focal_model$trans.matrix)[1] <= 4){
    next
  }
  
  tA <- mean(M_pars[1:3])
  tB <- mean(M_pars[14:16])
  fA <- mean(M_pars[4:6])
  fB <- mean(M_pars[17:19])
  dA <- mean(quickConvert2(M_pars[1:3], M_pars[4:6]))
  dB <- mean(quickConvert2(M_pars[14:16], M_pars[17:19]))
  qA <- mean(M_pars[7:10])
  qB <- mean(M_pars[20:23])
  qA_ln <- exp(mean(log(M_pars[7:10])))
  qB_ln <- exp(mean(log(M_pars[20:23])))
  
  tm <- (qA - qB)/(tA - tB)
  dm <- (qA - qB)/(dA - dB)
  tm_ln <- (qA_ln - qB_ln)/(tA - tB)
  dm_ln <- (qA_ln - qB_ln)/(dA - dB)
  clade_i <- gsub(".*res_state(.*)_idx.*", "\\1", nm_idx1[i])
  
  tmp_table <- data.frame(
    clade = clade_i,
    nTip = Ntip(all_res[[i]][[1]]$phy), 
    bestModel = best_model, 
    nRateClass = dim(focal_model$trans.matrix)[1]/4,
    AICwt = max(model_table[i,]),
    meanTransA = qA,
    meanTurnA = tA,
    meanNetDivA = dA,
    meanTransB = qB,
    meanTurnB = tB,
    meanNetDivB = dB,
    slopeTurn = tm,
    slopeNetDiv = dm
  )
  
  summary_table <- rbind(summary_table, tmp_table)
  slope_table <- rbind(slope_table, (c(t = tm, d = dm)))
}

summary_table[,-(1:4)] <- round(summary_table[,-(1:4)], 3)
summary_table <- summary_table[order(
  -summary_table$slopeTurn,
  -summary_table$slopeNetDiv
), ]

summary_table$clade <- gsub("-.*", "", summary_table$clade)
summary_table$clade[12] <- "Mimosoids"
write.csv(summary_table, "tables/summary_table.csv", row.names = FALSE)

SIGN.test(slope_table[,1], md = 0)
SIGN.test(slope_table[,2], md = 0)
