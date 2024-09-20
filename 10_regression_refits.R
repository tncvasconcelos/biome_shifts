# rm(list=ls())
setwd("~/biome_shifts/")
source("00_utility_functions.R")

library(ape)
library(hisse)
library(parallel)
library(phytools)
library(phylolm)
library(phytools)
library(RColorBrewer)

##############################
### load data
##############################
phy_bb <- read.tree("backbone_tree.tre")
phy_bb$tip.label <- gsub("-.*", "", phy_bb$tip.label)
phy_bb$tip.label <- gsub("_.*", "", phy_bb$tip.label)
phy_bb$tip.label <- clade_names[match(phy_bb$tip.label, gsub(" .*", "", clade_names))]
phy_bb$tip.label[is.na(match(phy_bb$tip.label, gsub(" .*", "", clade_names)))]
phy_bb <- drop.tip(phy_bb, which(is.na(phy_bb$tip.label)))

# models
to_load <- dir("5_results/", full.names = TRUE)
to_load <- to_load[-grep("Quercus", to_load)]
to_load <- to_load[-grep("Araceae", to_load)]

all_model_list <- list()
for(i in 1:length(to_load)){
  load(to_load[i])
  all_model_list[[i]] <- res
}

clade_names <- unlist(lapply(strsplit(to_load, "results_"), function(x) gsub(".RData", "", x[2])))
clade_names <- unlist(lapply(strsplit(clade_names, "-"), function(x) x[[1]]))
clade_names <- unlist(lapply(strsplit(clade_names, "_"), function(x) x[[1]]))
# clade_names <- paste0(clade_names, " (", n_tips, ")")
clade_names_res <- clade_names
names(all_model_list) <- clade_names_res

# recons
to_load <- dir("6_recons/", full.names = TRUE)
to_load <- to_load[-grep("Quercus", to_load)]
to_load <- to_load[-grep("Araceae", to_load)]

all_recon_list <- list()
for(i in 1:length(to_load)){
  load(to_load[i])
  all_recon_list[[i]] <- recon
}

length(all_recon_list) # 49 clades
length(all_recon_list[[1]]) # 36 models

model_par_list <- list()
for(clade_i in 1:length(all_recon_list)){
  clade_names[clade_i]
  recon_index <- !unlist(lapply(all_recon_list[[clade_i]], is.null))
  recon_AIC_wt <- GetAICWeights(all_model_list[[clade_i]][recon_index])
  model_pars <- lapply(all_model_list[[clade_i]][recon_index], convert_muhisse_pars)
  tip_labels <- all_model_list[[clade_i]][[1]]$phy$tip.label
  tip_states <- rowSums(all_model_list[[clade_i]][[1]]$data[,c(2,3)])
  tip_states <- setNames(tip_states, all_model_list[[clade_i]][[1]]$data[,1])
  tip_states <- tip_states[tip_labels]
  print(recon_AIC_wt[which.max(recon_AIC_wt)])
  model_par_list[[clade_i]] <- model_pars[[which.max(recon_AIC_wt)]]
}
names(model_par_list) <- clade_names


quickConvert <- function(turn_rate, eps_rate){
  return(convertBetweenPars(c(NA, NA, NA, turn_rate, eps_rate)))
}

getTipRateMat <- function(model_par){
  if(length(model_par$turnover) > 4){
    focal_mat_a <- model_par$trans_matrix[c(1,2,4),c(1,2,4)]
    focal_mat_b <- model_par$trans_matrix[c(5,6,8),c(5,6,8)]
    trans_rates <- c(rowSums(focal_mat_a, na.rm = TRUE), rowSums(focal_mat_b, na.rm = TRUE))
    turn_rates <- model_par$turnover[c(1,2,4,5,6,8)]
    eps_rates <- model_par$eps[c(1,2,4,5,6,8)]
    div_rate_mat <- mapply(quickConvert, turn_rates, eps_rates)
    tip_recon_mat <- all_recon_list[[clade_i]][recon_index][[2]]$tip.mat[,c(2,3,5,6,7,9)]
    rate_mat <- rbind(trans_rates, div_rate_mat)
    tip_rate_mat <- apply(rate_mat, 1, function(x) colSums(x * t(tip_recon_mat)))
  }else{
    focal_mat_a <- model_par$trans_matrix[c(1,2,4),c(1,2,4)]
    trans_rates <- c(rowSums(focal_mat_a, na.rm = TRUE))
    turn_rates <- model_par$turnover[c(1,2,4)]
    eps_rates <- model_par$eps[c(1,2,4)]
    div_rate_mat <- mapply(quickConvert, turn_rates, eps_rates)
    tip_recon_mat <- all_recon_list[[clade_i]][recon_index][[2]]$tip.mat[,c(2,3,5)]
    rate_mat <- rbind(trans_rates, div_rate_mat)
    rate_mat <- cbind(rate_mat, matrix(NA, 6, 3))
    tip_rate_mat <- apply(rate_mat, 1, function(x) colSums(x * t(tip_recon_mat)))
  }
  return(list(tip_rate_mat=tip_rate_mat, rate_mat=rate_mat))
}

model_mats <- lapply(model_par_list, getTipRateMat)
tip_rate_mats <- lapply(model_mats, "[[", "tip_rate_mat")

rate_mats <- lapply(model_mats, "[[", "rate_mat")
rate_mats_std <- lapply(rate_mats, function(x) apply(x, 1, function(y) (y - mean(y))/sd(y)))
rate_mat_df <- cbind(clade = rep(clade_names, each = 6), 
  par_name = rep(rownames(rate_mats[[1]]), length(clade_names)),
  as.data.frame(do.call(rbind, rate_mats)))
rownames(rate_mat_df) <- NULL


total_log_diff <- lapply(rate_mats, 
  function(x) log(rowMeans(x[,c(1:3)], na.rm = TRUE)) - log(rowMeans(x[,c(4:6)], na.rm = TRUE)))
total_log_diff <- as.data.frame(do.call(rbind, total_log_diff))
total_log_diff$trans_rates[total_log_diff$trans_rates == 0] <- NaN
total_log_diff$net.div[total_log_diff$net.div==0] <- NaN
total_log_diff$turn[total_log_diff$turn==0] <- NaN


##################### PLOT1

phy_bb$tip.label <- gsub(" .*", "", phy_bb$tip.label)
total_log_diff <- total_log_diff[phy_bb$tip.label,]

fit_r = phylolm(net.div ~ trans_rates, data=total_log_diff, phy=phy_bb, boot = 1000)
fit_t = phylolm(turn ~ trans_rates, data=total_log_diff, phy=phy_bb, boot = 1000)



# scatter plot of diff
all_fits_t <- list()
all_fits_r <- list()
for(j in 1:1000){
  print(j)
  no_to_swap <- sample(1:length(rate_mats), 1)
  to_swap <- sample(1:length(rate_mats), no_to_swap, FALSE)
  d_nets <- d_turns <- d_trans <- c()
  rate_mats2 <- rate_mats
  for(i in 1:length(rate_mats)){
    if(i %in% to_swap){
      rate_mats2[[i]] <- cbind(rate_mats2[[i]][,4:6], rate_mats2[[i]][,1:3])
    }else{
      rate_mats2[[i]] <- rate_mats[[i]]
    }
  }
  total_log_diff <- lapply(rate_mats2, 
    function(x) log(rowMeans(x[,c(4:6)], na.rm = TRUE)) - log(rowMeans(x[,c(1:3)], na.rm = TRUE)))
  total_log_diff <- as.data.frame(do.call(rbind, total_log_diff))
  total_log_diff$trans_rates[total_log_diff$trans_rates == 0] <- NaN
  total_log_diff$net.div[total_log_diff$net.div==0] <- NaN
  total_log_diff$turn[total_log_diff$turn==0] <- NaN
  total_log_diff <- total_log_diff[phy_bb$tip.label,]
  fit_r = phylolm(net.div ~ trans_rates, data=total_log_diff, phy=phy_bb, boot = 0)
  fit_t = phylolm(turn ~ trans_rates, data=total_log_diff, phy=phy_bb, boot = 0)
  all_fits_t[[j]] <- fit_t
  all_fits_r[[j]] <- fit_r
}


all_ps_t <- do.call(rbind, lapply(all_fits_t, function(x) summary(x)$coefficients[,4]))
all_est_t <- do.call(rbind, lapply(all_fits_t, function(x) summary(x)$coefficients[,1]))

all_ps_r <- do.call(rbind, lapply(all_fits_r, function(x) summary(x)$coefficients[,4]))
all_est_r <- do.call(rbind, lapply(all_fits_r, function(x) summary(x)$coefficients[,1]))


median(all_ps_t[,2])
median(all_est_t[,2])

median(all_ps_r[,2])
median(all_est_r[,2])

sum(all_ps_t[,2] <= 0.15)

hist(all_ps[,2])
abline(v = 0.05, col = "red")
abline(v = median(all_ps_t[,2]), col = "blue")
