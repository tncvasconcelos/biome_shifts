# rm(list=ls())
# setwd("~/Desktop/biome_shifts")
setwd("~/biome_shifts/")
source("00_utility_functions.R")

library(ape)
library(hisse)
library(parallel)
library(phytools)
library(phylolm)

plot_data <- read.csv("all_par_table.csv")

##############################
### Creating tip rates table
##############################
# models
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

# for bar charts
bar_dat <- data.frame(row.names = gsub(" .*", "", plot_data[,1]),
                      trans_a = log(rowMeans(rate_class_a_rate)),
                      trans_b = log(rowMeans(rate_class_b_rate)),
                      turns_a = log(rowMeans(rate_class_a_turn)),
                      turns_b = log(rowMeans(rate_class_b_turn)))

# for obs turnover
obs_dat <- data.frame(row.names = gsub(" .*", "", plot_data[,1]),
                      turns_00 = log(rowMeans(obs_00_turn)),
                      turns_01 = log(rowMeans(obs_01_turn)),
                      turns_11 = log(rowMeans(obs_01_turn)))


lm_dat <- lm_dat[phy_bb$tip.label,]
bar_dat <- bar_dat[phy_bb$tip.label,]
obs_dat <- obs_dat[phy_bb$tip.label,]

# plot(lm_dat)
# boxplot(bar_dat)
# fit model
fit = phylolm(d_turns ~ d_trans, data=lm_dat, phy=phy_bb, boot = 1000)
summary(fit)

# determine outliers
outliers <- unique(c(which(lm_dat$d_trans < quantile(lm_dat$d_trans, probs = 0.05)),
                     which(lm_dat$d_trans > quantile(lm_dat$d_trans, probs = 0.95)),
                     which(lm_dat$d_turns < quantile(lm_dat$d_turns, probs = 0.05)),
                     which(lm_dat$d_turns > quantile(lm_dat$d_turns, probs = 0.95))))
outlier_points <- lm_dat[outliers,c(1,2)]

# begin custom regression plot
dev.off()
pdf("plots/trans-turn-regression-color.pdf", width = 10, height=7)
plot.new()
# Set up the full plot area
par(fig=c(0, 1, 0, 1), mar=c(.1,.1,.1,.1), new=TRUE)
plot(type="n", x = lm_dat$d_trans, y = lm_dat$d_turns, bty='n',
     xaxt='n', yaxt='n', xlab = "", ylab = "", 
     xlim=range(lm_dat$d_trans))
grid()
# axes lines
lines(x = round(c(min(lm_dat$d_trans), max(lm_dat$d_trans))), y = c(0, 0), lwd=2)
lines(x = c(0, 0), y = round(c(min(lm_dat$d_turns), max(lm_dat$d_turns))), lwd=2)

# custom x
ticks_x <- round(seq(from = min(lm_dat$d_trans), to = max(lm_dat$d_trans+1), by = 2))
ticks_x <- ticks_x[-c(3,4)]
segments(x0 = ticks_x, y0 = -0.05, x1 = ticks_x, y1 = 0.05) 
text(x = ticks_x, y = -0.2, labels = round(ticks_x, 2), srt = 0, adj = 0.5)
text(x = max(lm_dat$d_trans), y = -0.6, labels = expression(Delta * log(transition~rates)), srt = 0, adj = 1, cex = 1.25)

# custom y
ticks_y <- round(seq(from = min(lm_dat$d_turns), to = max(lm_dat$d_turns), by = 2))
ticks_y <- ticks_y[-c(3,4,5)]
segments(x0 = -0.05, y0 = ticks_y, x1 = 0.05, y1 = ticks_y)
text(x = 0.2, y = ticks_y, labels = round(ticks_y, 2), srt = 0, adj = 0, xpd = TRUE)
text(x = -0.4, y = max(lm_dat$d_turns), labels = expression(Delta * log(turnover~rates)), srt = 90, adj = 1, cex = 1.25)

# prepare CI data 
X <- cbind(1, 
           seq(from=min(lm_dat$d_trans)*2, to=max(lm_dat$d_trans)*2, length.out=100))
Y_hat <- X %*% fit$coefficients  # Predicted values
Var_Y_hat <- diag(X %*% fit$vcov %*% t(X))
t_value <- qt(0.975, df = fit$n - length(fit$coefficients))
CI_lower <- Y_hat - t_value * sqrt(Var_Y_hat + fit$sigma2)
CI_upper <- Y_hat + t_value * sqrt(Var_Y_hat + fit$sigma2)
CI_results <- data.frame(Y_hat = as.vector(Y_hat), 
                         CI_lower = as.vector(CI_lower), 
                         CI_upper = as.vector(CI_upper))
X_plot <- X[,2]
sorted_indices <- order(X_plot) 
X_plot <- X_plot[sorted_indices]
sorted_Y_hat <- CI_results$Y_hat[sorted_indices]
sorted_CI_lower <- CI_results$CI_lower[sorted_indices]
sorted_CI_upper <- CI_results$CI_upper[sorted_indices]
mat_data <- cbind(sorted_Y_hat, sorted_CI_lower, sorted_CI_upper)

# plotting
matlines(X_plot, mat_data, lty = c(1, 2, 2), col = c("#737373", "#6BAED6", "#6BAED6"), lwd = c(2,1,1))
polygon(c(X_plot, rev(X_plot)), c(sorted_CI_lower, rev(sorted_CI_upper)), col = rgb(0.6784314, 0.8470588, 0.9019608, .25), border = NA)

# add points
palette <- brewer.pal(12, "Set3")
colors <- rep(palette, length.out = 51)
# Define the corners of the rectangle
x1 <- 2; y1 <- -6   # Bottom left
x2 <- 6; y2 <- -6  # Bottom right
x3 <- 2; y3 <- -2  # Top left
x4 <- 6; y4 <- -2 # Top right
n_col <- 5  # Columns
n_row <- 10 # Rows
x_seq <- seq(x1, x2, length.out = n_col)
y_seq <- seq(y1, y3, length.out = n_row)
grid_points <- expand.grid(x = x_seq, y = y_seq)
for(i in seq_along(clade_names_sorted)){
  points(x = lm_dat$d_trans[i], y = lm_dat$d_turns[i],
         pch = 21, bg = colors[i], cex=3)
  text(x = lm_dat$d_trans[i], y = lm_dat$d_turns[i],
       label=i, cex = .85)
  # points(x = grid_points$x[i], y = grid_points$y[i],
  #        pch = 21, bg = colors[i], cex=1)
  # text(x = grid_points$x[i], y = grid_points$y[i],
  #      label=clade_names_sorted[i], cex = .25)
}

dev.off()
# add outlier titles
# for(i in 1:nrow((outlier_points))){
#   text(x = outlier_points[i,1], y = outlier_points[i,2] + 0.3, labels = rownames(outlier_points)[i])
# }

# loc<- locator(1)
# print(loc)
# outlier_points$image_x <- c(-1.7, -3.28, -4, 4, 2, 6, 1, -2.3, -1.2, 1.7)
# outlier_points$image_y <- c(2.5, -2.8, 2, 3.2, 5.3, 4.4, -5.1, -5.2, -3.2, 2.8)
# outlier_points$names <- rownames(outlier_points)
# placePolygon <- function(x_init, y_init, x_poly, y_poly, pointName, side=1){
#   x_init <- as.numeric(x_init)
#   y_init <- as.numeric(y_init)
#   x_poly <- as.numeric(x_poly)
#   y_poly <- as.numeric(y_poly)
#   left <- x_poly - side / 2
#   bottom <- y_poly - side / 2
#   right <- x_poly + side / 2
#   top <- y_poly + side / 2
#   segments(x_init, y_init, x_poly, y_poly, col = "black", lty = 2)
#   polygon(c(left, left, right, right), y = c(top, bottom, bottom, top), col = "white")
#   text(x_poly, top + 0.05 * (top - bottom), labels = pointName, pos = 3)
# }
# apply(outlier_points, 1, function(x) placePolygon(x[1], x[2], x[3], x[4], x[5], 1))
# 


# inlaid boxplot
# par(fig=c(0.495, 0.99, 0.05, 0.38), new=TRUE, mar=c(.1, .1, .1, .1))
# dev.off()
pdf("plots/trans-turn-boxplot.pdf", height=5, width=8)
boxplot(bar_dat, main="", xlab="", ylab="", axes=FALSE, outline = FALSE,
        col=c("#b10026", "#fc9272", "#034e7b", "#74a9cf"),
        ylim=range(bar_dat)*c(1.1, 1),
        xlim=c(-.5,5))
# box()
ticks_y <- round(seq(from = min(bar_dat), to = max(bar_dat), by = 2))
segments(x0 = 0, y0 = ticks_y, x1 = 0.2, y1 = ticks_y)
segments(x0 = 0.1, y0 = min(ticks_y), x1 = 0.1, y1 = max(ticks_y))
text(x = -0.1, y = ticks_y, labels = round(ticks_y, 2), 
     srt = 0, adj = 1, xpd = TRUE, cex=0.75)
text(x = -0.35, y = mean(range(bar_dat)), labels = expression(log(rate)), srt = 90, adj = 0.5, cex = 1)
text(x = 1, y = min(bar_dat)*1.1, labels = expression(q[A]))
text(x = 2, y = min(bar_dat)*1.1, labels = expression(q[B]))
text(x = 3, y = min(bar_dat)*1.1, labels = expression(tau[A]))
text(x = 4, y = min(bar_dat)*1.1, labels = expression(tau[B]))
points(x = rep(1, dim(bar_dat)[1]) + rnorm(dim(bar_dat)[1], sd = 0.05), 
       y = bar_dat$trans_a, pch = 21, bg="#b10026")
points(x = rep(2, dim(bar_dat)[1]) + rnorm(dim(bar_dat)[1], sd = 0.05), 
       y = bar_dat$trans_b, pch = 21, bg="#fc9272")
points(x = rep(3, dim(bar_dat)[1]) + rnorm(dim(bar_dat)[1], sd = 0.05), 
       y = bar_dat$turns_a, pch = 21, bg="#034e7b")
points(x = rep(4, dim(bar_dat)[1]) + rnorm(dim(bar_dat)[1], sd = 0.05), 
       y = bar_dat$turns_b, pch = 21, bg="#74a9cf")
dev.off()
# Reset to full plot area
# par(fig=c(0, 1, 0, 1), new=TRUE)
# rasterImage(img, left, bottom, right, top)

#### observed state transi

pdf("plots/turns-obs-boxplot.pdf", height=5, width=10)
brew_col <- brewer.pal(9, "Set1")[c(3,6,5)]
cols <- setNames(brew_col, c("closed", "widespread", "open"))
boxplot(obs_dat, main="", xlab="", ylab="", axes=FALSE, outline = FALSE,
        col=cols, xlim=c(0, 4), ylim=range(obs_dat)*c(1.1, 1))
# box()
ticks_y <- round(seq(from = min(obs_dat), to = max(obs_dat), by = 2))
segments(x0 = 0.32, y0 = ticks_y, x1 = 0.18, y1 = ticks_y)
segments(x0 = 0.25, y0 = min(ticks_y), x1 = 0.25, y1 = max(ticks_y))
text(x = 0.15, y = ticks_y, labels = round(ticks_y, 2), 
     srt = 0, adj = 1, xpd = TRUE, cex=0.75)
text(x = 0, y = mean(range(obs_dat)), labels = expression(log(tau)), srt = 90, adj = 0.5, cex = 1)
text(x = 1, y = min(obs_dat)*1.15, labels = expression(C), cex=1)
text(x = 2, y = min(obs_dat)*1.15, labels = expression(W), cex=1)
text(x = 3, y = min(obs_dat)*1.15, labels = expression(O), cex=1)
points(x = rep(1, dim(obs_dat)[1]) + rnorm(dim(obs_dat)[1], sd = 0.05), 
       y = obs_dat$turns_00, pch = 21, bg=cols[1])
points(x = rep(2, dim(obs_dat)[1]) + rnorm(dim(obs_dat)[1], sd = 0.05), 
       y = obs_dat$turns_01, pch = 21, bg=cols[2])
points(x = rep(3, dim(obs_dat)[1]) + rnorm(dim(obs_dat)[1], sd = 0.05), 
       y = obs_dat$turns_11, pch = 21, bg=cols[3])
dev.off()


##############################
### model avg ancestral states
##############################

# finished model sets for particular datsets
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

##############################
### load tree
##############################

phy_bb <- read.tree("backbone_tree.tre")
phy_bb$tip.label <- gsub("-.*", "", phy_bb$tip.label)
phy_bb$tip.label <- gsub("_.*", "", phy_bb$tip.label)
phy_bb$tip.label <- clade_names[match(phy_bb$tip.label, gsub(" .*", "", clade_names))]
phy_bb$tip.label[is.na(match(phy_bb$tip.label, gsub(" .*", "", clade_names)))]
phy_bb <- drop.tip(phy_bb, which(is.na(phy_bb$tip.label)))

##############################
### begin analysis
##############################

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
total_log_diff <- total_log_diff[phy_bb$tip.label,]


##################### PLOT1


pdf("plots/trans-turn-net-div-regression-color.pdf", width = 10, height=7)
par(fig=c(0, 1, 0, 1), mar=c(.5,.5,.5,.5), new=TRUE, mfrow=c(2,1))
trans_rates <- total_log_diff$trans_rates
div_rates <- total_log_diff$net.div
fit = phylolm(net.div ~ trans_rates, data=total_log_diff, phy=phy_bb, boot = 1000)
summary(fit)

# Call:
#   phylolm(formula = net.div ~ trans_rates, data = total_log_diff, 
#     phy = phy_bb, boot = 1000)
# 
# AIC logLik 
# 96.32 -45.16 
# 
# Raw residuals:
#   Min      1Q  Median      3Q     Max 
# -2.0630 -1.2386 -0.1461  1.3241  2.5119 
# 
# Mean tip height: 135.758
# Parameter estimate(s) using ML:
#   sigma2: 0.02396223 
# 
# Coefficients:
#               Estimate    StdErr   t.value lowerbootCI upperbootCI p.value
# (Intercept) -0.064810  0.656694 -0.098692   -1.274527      1.1104  0.9223
# trans_rates -0.018468  0.098798 -0.186923   -0.200141      0.1576  0.8534
# 
# R-squared: 0.001586	Adjusted R-squared: -0.0438 
# 
# sigma2: 0.02396223
# bootstrap mean: 0.02239889 (on raw scale)
# 0.02144528 (on log scale, then back transformed)
# bootstrap 95% CI: (0.01168352,0.03718905)
# 
# Parametric bootstrap results based on 1000 fitted replicates

# Set up the full plot area
plot(type="n", x = fit$X[,2], y = fit$y, bty='n',
  xaxt='n', yaxt='n', xlab = "", ylab = "", 
  xlim=range(fit$X[,2]))
grid()
# axes lines
lines(x = round(c(min(fit$X[,2]), max(fit$X[,2]))), y = c(0, 0), lwd=2)
lines(x = c(0, 0), y = round(c(min(fit$X[,2]), max(fit$X[,2]))), lwd=2)
tick_size <- 0.1

# custom x
ticks_x <- round(seq(from = min(fit$X[,2]), to = max(fit$X[,2], na.rm = TRUE), by = 2))
ticks_x <- ticks_x[ticks_x != 0]
segments(x0 = ticks_x, y0 = -tick_size, x1 = ticks_x, y1 = tick_size) 
text(x = ticks_x, y = -tick_size*3, labels = round(ticks_x, 2), srt = 0, adj = 0.5)
text(x = max(fit$X[,2]), y = -tick_size * 6, labels = expression(Delta * log(transition~rates)), srt = 0, adj = 1, cex = 1)

# custom y
ticks_y <- round(seq(from = min(fit$y), to = max(fit$y), by = 2))
ticks_y <- ticks_y[ticks_y!=0]
segments(x0 = -tick_size, y0 = ticks_y, x1 = tick_size, y1 = ticks_y)
text(x = -tick_size*4, y = ticks_y, labels = round(ticks_y, 2), srt = 0, adj = 0, xpd = TRUE)
text(x = -tick_size*6, y = max(fit$y), labels = expression(Delta * log(net~div~rates)), srt = 90, adj = 1, cex = 1)

# prepare CI data 
X <- cbind(1, 
  seq(from=min(fit$X[,2])*2, to=max(fit$X[,2])*2, length.out=100))
Y_hat <- X %*% fit$coefficients  # Predicted values
Var_Y_hat <- diag(X %*% fit$vcov %*% t(X))
t_value <- qt(0.975, df = fit$n - length(fit$coefficients))
CI_lower <- Y_hat - t_value * sqrt(Var_Y_hat + fit$sigma2)
CI_upper <- Y_hat + t_value * sqrt(Var_Y_hat + fit$sigma2)
CI_results <- data.frame(Y_hat = as.vector(Y_hat), 
  CI_lower = as.vector(CI_lower), 
  CI_upper = as.vector(CI_upper))
X_plot <- X[,2]
sorted_indices <- order(X_plot) 
X_plot <- X_plot[sorted_indices]
sorted_Y_hat <- CI_results$Y_hat[sorted_indices]
sorted_CI_lower <- CI_results$CI_lower[sorted_indices]
sorted_CI_upper <- CI_results$CI_upper[sorted_indices]
mat_data <- cbind(sorted_Y_hat, sorted_CI_lower, sorted_CI_upper)

# plotting
matlines(X_plot, mat_data, lty = c(1, 2, 2), col = c("blue", "#6BAED6", "#6BAED6"), lwd = c(2,1,1))
polygon(c(X_plot, rev(X_plot)), c(sorted_CI_lower, rev(sorted_CI_upper)), col = rgb(0.6784314, 0.8470588, 0.9019608, .25), border = NA)

# add points
palette <- RColorBrewer::brewer.pal(12, "Set3")
colors <- rep(palette, length.out = 51)
for(i in seq_along(phy_bb$tip.label)){
  points(x = trans_rates[i], y = div_rates[i],
    pch = 21, bg = colors[i], cex=4)
  text(x = trans_rates[i], y = div_rates[i],
    label=i, cex = 1.5)
}

##################### PLOT2
trans_rates <- total_log_diff$trans_rates
div_rates <- total_log_diff$turn
fit = phylolm(turn ~ trans_rates, data=total_log_diff, phy=phy_bb, boot = 1000)
summary(fit)

# Call:
#   phylolm(formula = turn ~ trans_rates, data = total_log_diff, 
#     phy = phy_bb, boot = 1000)
# 
# AIC logLik 
# 116.82 -55.41 
# 
# Raw residuals:
#   Min      1Q  Median      3Q     Max 
# -5.4750 -1.1323 -0.1888  1.3494  2.6908 
# 
# Mean tip height: 135.758
# Parameter estimate(s) using ML:
#   sigma2: 0.04668805 
# 
# Coefficients:
#               Estimate    StdErr   t.value lowerbootCI upperbootCI p.value  
# (Intercept) -0.423386  0.909908 -0.465306   -2.041813      1.4683 0.64609  
# trans_rates  0.277288  0.133809  2.072265    0.030342      0.5336 0.04963 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-squared: 0.1573	Adjusted R-squared: 0.1207 
# 
# sigma2: 0.04668805
# bootstrap mean: 0.04299894 (on raw scale)
# 0.04126994 (on log scale, then back transformed)
# bootstrap 95% CI: (0.02200801,0.06910136)
# 
# Parametric bootstrap results based on 1000 fitted replicates

# Set up the full plot area
plot(type="n", x = fit$X[,2], y = fit$y, bty='n',
  xaxt='n', yaxt='n', xlab = "", ylab = "", 
  xlim=range(fit$X[,2]))
grid()
# axes lines
lines(x = round(c(min(fit$X[,2]), max(fit$X[,2]))), y = c(0, 0), lwd=2)
lines(x = c(0, 0), y = round(c(min(fit$X[,2]), max(fit$X[,2]))), lwd=2)
tick_size <- .2
# custom x
ticks_x <- round(seq(from = min(fit$X[,2]), to = max(fit$X[,2], na.rm = TRUE), by = 2))
ticks_x <- ticks_x[ticks_x!=0]
segments(x0 = ticks_x, y0 = -tick_size, x1 = ticks_x, y1 = tick_size) 
text(x = ticks_x, y = -tick_size*3, labels = round(ticks_x, 2), srt = 0, adj = 0.5)
text(x = max(fit$X[,2]), y = -tick_size * 6, labels = expression(Delta * log(transition~rates)), srt = 0, adj = 1, cex = 1)

# custom y
ticks_y <- round(seq(from = min(fit$y), to = max(fit$y), by = 2))
segments(x0 = -tick_size*.5, y0 = ticks_y, x1 = tick_size*.5, y1 = ticks_y)
text(x = -tick_size*2, y = ticks_y, labels = round(ticks_y, 2), srt = 0, adj = 0, xpd = TRUE)
text(x = -tick_size*4, y = min(fit$y), labels = expression(Delta * log(turn~rates)), srt = 90, adj = 0, cex = 1)

# prepare CI data 
X <- cbind(1, 
  seq(from=min(fit$X[,2])*2, to=max(fit$X[,2])*2, length.out=100))
Y_hat <- X %*% fit$coefficients  # Predicted values
Var_Y_hat <- diag(X %*% fit$vcov %*% t(X))
t_value <- qt(0.975, df = fit$n - length(fit$coefficients))
CI_lower <- Y_hat - t_value * sqrt(Var_Y_hat + fit$sigma2)
CI_upper <- Y_hat + t_value * sqrt(Var_Y_hat + fit$sigma2)
CI_results <- data.frame(Y_hat = as.vector(Y_hat), 
  CI_lower = as.vector(CI_lower), 
  CI_upper = as.vector(CI_upper))
X_plot <- X[,2]
sorted_indices <- order(X_plot) 
X_plot <- X_plot[sorted_indices]
sorted_Y_hat <- CI_results$Y_hat[sorted_indices]
sorted_CI_lower <- CI_results$CI_lower[sorted_indices]
sorted_CI_upper <- CI_results$CI_upper[sorted_indices]
mat_data <- cbind(sorted_Y_hat, sorted_CI_lower, sorted_CI_upper)

# plotting
matlines(X_plot, mat_data, lty = c(1, 2, 2), col = c("blue", "#6BAED6", "#6BAED6"), lwd = c(2,1,1))
polygon(c(X_plot, rev(X_plot)), c(sorted_CI_lower, rev(sorted_CI_upper)), col = rgb(0.6784314, 0.8470588, 0.9019608, .25), border = NA)

# add points
palette <- RColorBrewer::brewer.pal(12, "Set3")
colors <- rep(palette, length.out = 51)
for(i in seq_along(phy_bb$tip.label)){
  points(x = trans_rates[i], y = div_rates[i],
    pch = 21, bg = colors[i], cex=4)
  text(x = trans_rates[i], y = div_rates[i],
    label=i, cex = 1.5)
}
dev.off()

##################### FIN

######### REMOVING OUTLEIRS
z_index <- abs(apply(total_log_diff, 2, function(x)( x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE))) > 2
total_log_diff[z_index] <- NA
trans_rates <- total_log_diff$trans_rates
div_rates <- total_log_diff$turn
fit = phylolm(turn ~ trans_rates, data=total_log_diff, phy=phy_bb, boot = 1000)
summary(fit)

# Call:
#   phylolm(formula = turn ~ trans_rates, data = total_log_diff, 
#     phy = phy_bb, boot = 1000)
# 
# AIC logLik 
# 109.71 -51.85 
# 
# Raw residuals:
#   Min      1Q  Median      3Q     Max 
# -5.5546 -1.3558 -0.2258  1.4759  2.6981 
# 
# Mean tip height: 135.758
# Parameter estimate(s) using ML:
#   sigma2: 0.05063446 
# 
# Coefficients:
#   Estimate    StdErr   t.value lowerbootCI upperbootCI p.value
# (Intercept) -0.401621  0.970294 -0.413916   -2.244797      1.2963  0.6831
# trans_rates  0.259218  0.188179  1.377513   -0.093784      0.6403  0.1829
# 
# R-squared: 0.08287	Adjusted R-squared: 0.0392 
# 
# sigma2: 0.05063446
# bootstrap mean: 0.04661747 (on raw scale)
# 0.04444588 (on log scale, then back transformed)
# bootstrap 95% CI: (0.02236352,0.07813109)
# 
# Parametric bootstrap results based on 1000 fitted replicates

plot(type="n", x = fit$X[,2], y = fit$y, bty='n',
  xaxt='n', yaxt='n', xlab = "", ylab = "", 
  xlim=range(fit$X[,2]))
grid()
# axes lines
lines(x = round(c(min(fit$X[,2]), max(fit$X[,2]))), y = c(0, 0), lwd=2)
lines(x = c(0, 0), y = round(c(min(fit$X[,2]), max(fit$X[,2]))), lwd=2)
tick_size <- .2
# custom x
ticks_x <- round(seq(from = min(fit$X[,2]), to = max(fit$X[,2], na.rm = TRUE), by = 2))
ticks_x <- ticks_x[ticks_x!=0]
segments(x0 = ticks_x, y0 = -tick_size, x1 = ticks_x, y1 = tick_size) 
text(x = ticks_x, y = -tick_size*3, labels = round(ticks_x, 2), srt = 0, adj = 0.5)
text(x = max(fit$X[,2]), y = -tick_size * 6, labels = expression(Delta * log(transition~rates)), srt = 0, adj = 1, cex = 1)

# custom y
ticks_y <- round(seq(from = min(fit$y), to = max(fit$y), by = 2))
segments(x0 = -tick_size*.5, y0 = ticks_y, x1 = tick_size*.5, y1 = ticks_y)
text(x = -tick_size*2, y = ticks_y, labels = round(ticks_y, 2), srt = 0, adj = 0, xpd = TRUE)
text(x = -tick_size*4, y = min(fit$y), labels = expression(Delta * log(turn~rates)), srt = 90, adj = 0, cex = 1)

# prepare CI data 
X <- cbind(1, 
  seq(from=min(fit$X[,2])*2, to=max(fit$X[,2])*2, length.out=100))
Y_hat <- X %*% fit$coefficients  # Predicted values
Var_Y_hat <- diag(X %*% fit$vcov %*% t(X))
t_value <- qt(0.975, df = fit$n - length(fit$coefficients))
CI_lower <- Y_hat - t_value * sqrt(Var_Y_hat + fit$sigma2)
CI_upper <- Y_hat + t_value * sqrt(Var_Y_hat + fit$sigma2)
CI_results <- data.frame(Y_hat = as.vector(Y_hat), 
  CI_lower = as.vector(CI_lower), 
  CI_upper = as.vector(CI_upper))
X_plot <- X[,2]
sorted_indices <- order(X_plot) 
X_plot <- X_plot[sorted_indices]
sorted_Y_hat <- CI_results$Y_hat[sorted_indices]
sorted_CI_lower <- CI_results$CI_lower[sorted_indices]
sorted_CI_upper <- CI_results$CI_upper[sorted_indices]
mat_data <- cbind(sorted_Y_hat, sorted_CI_lower, sorted_CI_upper)

# plotting
matlines(X_plot, mat_data, lty = c(1, 2, 2), col = c("blue", "#6BAED6", "#6BAED6"), lwd = c(2,1,1))
polygon(c(X_plot, rev(X_plot)), c(sorted_CI_lower, rev(sorted_CI_upper)), col = rgb(0.6784314, 0.8470588, 0.9019608, .25), border = NA)

# add points
palette <- RColorBrewer::brewer.pal(12, "Set3")
colors <- rep(palette, length.out = 51)
for(i in seq_along(phy_bb$tip.label)){
  points(x = trans_rates[i], y = div_rates[i],
    pch = 21, bg = colors[i], cex=4)
  text(x = trans_rates[i], y = div_rates[i],
    label=i, cex = 1.5)
}

##################### ##################### ##################### ##################### ##################### ##################### 
##################### ##################### boxplot by hidden rate class
##################### ##################### ##################### ##################### ##################### ##################### 
box_raw_rates <- do.call(rbind, lapply(rate_mats, 
  function(x) c(log(rowMeans(x[,c(1:3)], na.rm = TRUE)), log(rowMeans(x[,c(4:6)], na.rm = TRUE)))))
bar_dat <- box_raw_rates[,c(1,7,5,11,4,10)]
# inlaid boxplot
# par(fig=c(0.495, 0.99, 0.05, 0.38), new=TRUE, mar=c(.1, .1, .1, .1))
# dev.off()
pdf("plots/trans-turn-boxplot.pdf", height=8, width=5)

bar_dat <- bar_dat[!apply(bar_dat, 1, function(x) any(na.omit(x) < -10)),]

boxplot(bar_dat, main="", xlab="", ylab="", axes=FALSE, outline = FALSE,
  col=c("#b10026", "#fc9272", "#034e7b", "#74a9cf", "#FFCC00", "#FFEE8C"),
  ylim=range(bar_dat, na.rm = TRUE)*c(1.1, 1),
  xlim=c(-.5,7))
# box()
ticks_y <- round(seq(from = min(bar_dat, na.rm = TRUE), to = max(bar_dat, na.rm = TRUE), by = 2))
segments(x0 = 0, y0 = ticks_y, x1 = 0.2, y1 = ticks_y)
segments(x0 = 0.1, y0 = min(ticks_y), x1 = 0.1, y1 = max(ticks_y))
text(x = -0.1, y = ticks_y, labels = round(ticks_y, 2), 
  srt = 0, adj = 1, xpd = TRUE, cex=0.75)
text(x = -0.35, y = mean(range(bar_dat, na.rm = TRUE)), labels = expression(log(rate)), srt = 90, adj = 0.5, cex = 1)
text(x = 1, y = min(bar_dat, na.rm = TRUE)*1.1, labels = expression(q[A]))
text(x = 2, y = min(bar_dat, na.rm = TRUE)*1.1, labels = expression(q[B]))
text(x = 3, y = min(bar_dat, na.rm = TRUE)*1.1, labels = expression(tau[A]))
text(x = 4, y = min(bar_dat, na.rm = TRUE)*1.1, labels = expression(tau[B]))
text(x = 5, y = min(bar_dat, na.rm = TRUE)*1.1, labels = expression(r[A]))
text(x = 6, y = min(bar_dat, na.rm = TRUE)*1.1, labels = expression(r[B]))
points(x = rep(1, dim(bar_dat)[1]) + rnorm(dim(bar_dat)[1], sd = 0.05), 
  y = bar_dat[,1], pch = 21, bg="#b10026")
points(x = rep(2, dim(bar_dat)[1]) + rnorm(dim(bar_dat)[1], sd = 0.05), 
  y = bar_dat[,2], pch = 21, bg="#fc9272")
points(x = rep(3, dim(bar_dat)[1]) + rnorm(dim(bar_dat)[1], sd = 0.05), 
  y = bar_dat[,3], pch = 21, bg="#034e7b")
points(x = rep(4, dim(bar_dat)[1]) + rnorm(dim(bar_dat)[1], sd = 0.05), 
  y = bar_dat[,4], pch = 21, bg="#74a9cf")
points(x = rep(5, dim(bar_dat)[1]) + rnorm(dim(bar_dat)[1], sd = 0.05), 
  y = bar_dat[,5], pch = 21, bg="#FFCC00")
points(x = rep(6, dim(bar_dat)[1]) + rnorm(dim(bar_dat)[1], sd = 0.05), 
  y = bar_dat[,6], pch = 21, bg="#FFEE8C")

dev.off()
# Reset to full plot area# 

##################### ##################### ##################### ##################### 
##################### ##################### boxplot by observed state
##################### ##################### ##################### ##################### 

box_obs_raw_rates <- do.call(rbind, lapply(rate_mats, 
  function(x) c(log(rowMeans(x[,c(1,4)], na.rm = TRUE)), log(rowMeans(x[,c(2,5)], na.rm = TRUE)), log(rowMeans(x[,c(3,6)], na.rm = TRUE)))))
obs_dat <- box_obs_raw_rates[,c(1,7,13, 5,11,17, 4,10, 16)]
z_index <- abs(apply(obs_dat, 2, function(x)( x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE))) > 2
obs_dat[z_index] <- NA
# obs_dat <- obs_dat[apply(obs_dat, 1, function(x) !any(na.omit(x) < -10)),]

pdf("plots/turns-net-obs-boxplot.pdf", height=5, width=12)
brew_col <- RColorBrewer::brewer.pal(9, "Set1")[c(3,6,5)]
cols <- setNames(brew_col, c("closed", "widespread", "open"))

x_index <- c(seq(from = 1, to = 3, length.out = 3), seq(from = 5, to = 7, length.out = 3))

par(mar = c(5, 5, 5, 5))
plot(1, type = "n", xlim = c(0,8),
  ylim = range(obs_dat, na.rm = TRUE)*c(1.1, 1), xaxt = "n", xlab = "Rates by observed biome", ylab = expression(log(rate)),
  bty = "n")

axis(1, at = x_index[1:3],
  labels = c(NA, NA, NA))
axis(1, at = mean(x_index[1:3]),
  labels = expression(tau), tick = FALSE)

axis(1, at = x_index[4:6],
  labels = c(NA, NA, NA))
axis(1, at = mean(x_index[4:6]),
  labels = expression(r), tick=FALSE)


boxplot(obs_dat[,-c(1:3)], 
  at = x_index,
  col = rep(cols, 3),
  border = "black", 
  ylim=range(obs_dat, na.rm = TRUE)*c(1.1, 1), 
  outline=FALSE,
  add=TRUE, 
  boxwex = 0.5, axes=FALSE)

points(x = rep(x_index[1], dim(obs_dat)[1]) + rnorm(dim(obs_dat)[1], sd = 0.05), 
  y = obs_dat[,4], pch = 21, bg=cols[1])
points(x = rep(x_index[2], dim(obs_dat)[1]) + rnorm(dim(obs_dat)[1], sd = 0.05), 
  y = obs_dat[,5], pch = 21, bg=cols[2])
points(x = rep(x_index[3], dim(obs_dat)[1]) + rnorm(dim(obs_dat)[1], sd = 0.05), 
  y = obs_dat[,6], pch = 21, bg=cols[3])

points(x = rep(x_index[4], dim(obs_dat)[1]) + rnorm(dim(obs_dat)[1], sd = 0.05), 
  y = obs_dat[,7], pch = 21, bg=cols[1])
points(x = rep(x_index[5], dim(obs_dat)[1]) + rnorm(dim(obs_dat)[1], sd = 0.05), 
  y = obs_dat[,8], pch = 21, bg=cols[2])
points(x = rep(x_index[6], dim(obs_dat)[1]) + rnorm(dim(obs_dat)[1], sd = 0.05), 
  y = obs_dat[,9], pch = 21, bg=cols[3])


legend("topright", legend = names(cols), fill = cols, border = "black", 
  box.lwd = 1, title = "Biome type:", cex = 0.85, pt.cex = 0.85)
dev.off()

### trans rates
plot_dat <- read.csv("all_par_table.csv")
trns_dat <- plot_dat[,c(2:5)] + plot_dat[,c(6:9)]
rownames(trns_dat) <- plot_dat$Group

z_index <- abs(apply(log(trns_dat), 2, function(x)( x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE))) > 2
trns_dat[z_index] <- NA

pdf("plots/rates-obs-boxplot.pdf", height=5, width=8)
plot(1, type = "n", xlim = c(0,5),
  ylim = range(log(trns_dat), na.rm = TRUE)*c(1.1, 1), xaxt = "n", xlab = "Rates by observed biome", ylab = expression(log(rate)),
  bty = "n")

boxplot(log(trns_dat), 
  col = cols[c(1,2,2,3)],
  border = "black", 
  outline=FALSE,
  boxwex = 0.5, 
  axes = FALSE,
  add = TRUE)

axis(1, at = 1:4,
  labels = c(expression(q[cw]),expression(q[wc]),expression(q[wo]),expression(q[ow])))


points(x = 1 + rnorm(dim(trns_dat)[1], sd = 0.05), 
  y = log(trns_dat)[,1], pch = 21, bg=cols[1])
points(x = 2 + rnorm(dim(trns_dat)[1], sd = 0.05), 
  y = log(trns_dat)[,2], pch = 21, bg=cols[2])
points(x = 3 + rnorm(dim(trns_dat)[1], sd = 0.05), 
  y = log(trns_dat)[,3], pch = 21, bg=cols[2])
points(x = 4 + rnorm(dim(trns_dat)[1], sd = 0.05), 
  y = log(trns_dat)[,4], pch = 21, bg=cols[3])

legend("topright", legend = names(cols), fill = cols, border = "black", 
  box.lwd = 1, title = "Biome type:", cex = 0.85, pt.cex = 0.85)

dev.off()

### TTESTS

phy_bb_2 <- keep.tip(phy_bb, rownames(obs_dat))
run_ttest <- function(phy, dat_1, dat_2){
  na_index <- is.na(dat_1) | is.nan(dat_1) | is.na(dat_2) | is.nan(dat_2)
  dat_1 <- dat_1[!na_index]
  dat_2 <- dat_2[!na_index]
  phy <- keep.tip(phy, names(dat_1))
  return(phyl.pairedttest(phy, dat_1, dat_2))
}


closed_turn <- setNames(obs_dat[,4], rownames(obs_dat))
wide_turn <- setNames(obs_dat[,5], rownames(obs_dat))
open_turn <- setNames(obs_dat[,6], rownames(obs_dat))

closed_div <- setNames(obs_dat[,7], rownames(obs_dat))
wide_div <- setNames(obs_dat[,8], rownames(obs_dat))
open_div <- setNames(obs_dat[,9], rownames(obs_dat))

run_ttest(phy_bb, closed_turn, wide_turn)
run_ttest(phy_bb, closed_turn, open_turn)
run_ttest(phy_bb, open_turn, wide_turn)

run_ttest(phy_bb, closed_div, wide_div)
run_ttest(phy_bb, closed_div, open_div)
run_ttest(phy_bb, open_div, wide_div)

