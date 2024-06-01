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

plot_data <- read.csv("all_par_table.csv")
plot_data <- plot_data[plot_data$Group != "Quercus",]
plot_data <- plot_data[plot_data$Group != "Araceae",]

##############################
### correcting tip labels
##############################

phy_bb <- read.tree("backbone_tree.tre")
dat_names <- unlist(lapply(strsplit(paste(plot_data$Group), " "), function(x) x[[1]]))
phy_bb$tip.label <- unlist(lapply(strsplit(phy_bb$tip.label, "-"), function(x) x[[1]]))
phy_bb$tip.label <- unlist(lapply(strsplit(phy_bb$tip.label, "_"), function(x) x[[1]]))
to_drop <- phy_bb$tip.label[!phy_bb$tip.label %in% dat_names]
phy_bb <- drop.tip(phy_bb, to_drop)
# plotTree(phy_bb, fsize=0.75, ftype="i")
clade_names_sorted <- tip_names <- phy_bb$tip.label


##############################
### regression of mean rateclass rate against turnover
##############################
rate_class_a_rate <- cbind(plot_data$f00_t01_a,
                           plot_data$f01_t00_a,
                           plot_data$f01_t11_a,
                           plot_data$f11_t01_a)
rate_class_b_rate <- cbind(plot_data$f00_t01_b,
                           plot_data$f01_t00_b,
                           plot_data$f01_t11_b,
                           plot_data$f11_t01_b)
rate_class_a_turn <- cbind(plot_data$tn_00_a,
                           plot_data$tn_01_a,
                           plot_data$tn_11_a)
rate_class_b_turn <- cbind(plot_data$tn_00_b,
                           plot_data$tn_01_b,
                           plot_data$tn_11_b)
obs_00_turn <- cbind(plot_data$tn_00_a,
                     plot_data$tn_00_b)
obs_01_turn <- cbind(plot_data$tn_01_a,
                     plot_data$tn_01_b)
obs_11_turn <- cbind(plot_data$tn_11_a,
                     plot_data$tn_11_b)


# scatter plot of diff
lm_dat <- data.frame(row.names = gsub(" .*", "", plot_data[,1]),
                     d_trans = log(rowMeans(rate_class_b_rate)) - 
                       log(rowMeans(rate_class_a_rate)),
                     d_turns = log(rowMeans(rate_class_b_turn)) - 
                       log(rowMeans(rate_class_a_turn)))

init_fit = phylolm(d_turns ~ d_trans , data=lm_dat, phy=phy_bb, boot = 1000)
summary(init_fit)

plot(lm_dat)
abline(init_fit, col="red")
abline(a = 0, b = 1, type = 2)

all_fits <- list()
for(j in 1:10000){
  print(j)
  no_to_swap <- sample(1:49, 1)
  to_swap <- sample(1:49, no_to_swap, FALSE)
  d_turns <- d_trans <- c()
  for(i in 1:49){
    if(i %in% to_swap){
      d_trans <- (c(d_trans, log(rowMeans(rate_class_a_rate)[i]) - 
                      log(rowMeans(rate_class_b_rate)[i])))
      d_turns <- (c(d_turns, log(rowMeans(rate_class_a_turn)[i]) - 
                      log(rowMeans(rate_class_b_turn)[i])))
    }else{
      d_trans <- (c(d_trans, log(rowMeans(rate_class_b_rate)[i]) - 
                      log(rowMeans(rate_class_a_rate)[i])))
      d_turns <- (c(d_turns, log(rowMeans(rate_class_b_turn)[i]) - 
                      log(rowMeans(rate_class_a_turn)[i])))
    }
  }
  
  lm_dat_tmp <- data.frame(row.names = gsub(" .*", "", plot_data[,1]),
                       d_trans = d_trans,
                       d_turns = d_turns)
  fit_tmp = phylolm(d_turns ~ d_trans, data=lm_dat_tmp, phy=phy_bb, boot = 0)
  # plot(lm_dat, ylim=c(-6, 6), xlim=c(-6, 6))
  # abline(fit)
  # abline(v=0)
  # abline(h=0)
  all_fits[[j]] <- fit_tmp
  # Sys.sleep(.2)
}


tmp <- summary(fit)
tmp$coefficients
all_ps <- do.call(rbind, lapply(all_fits, function(x) summary(x)$coefficients[,4]))
all_est <- do.call(rbind, lapply(all_fits, function(x) summary(x)$coefficients[,1]))

colMeans(all_est)
hist(all_est[,2])
length(which(all_ps[,2] <= 0.05))

hist(all_ps[,2])
abline(v = 0.05, col = "red")
