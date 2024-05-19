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
fit = phylolm(d_turns ~ d_trans - 1, data=lm_dat, phy=phy_bb, boot = 1000)
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

# regression line and CI
slope <- coef(fit)
slope_ci <- confint(fit)

x_range <- seq(from=min(lm_dat$d_trans)*2, to=max(lm_dat$d_trans)*2, length.out=100)
y_main = x_range * slope
y_lower = x_range * slope_ci[1]
y_upper = x_range * slope_ci[2]

polygon(x=c(x_range, rev(x_range)), y=c(y_lower, rev(y_upper)), col=rgb(0.6784314, 0.8470588, 0.9019608, .25), border=NA)

abline(a = 0, b = slope, col = "#737373", lwd = 2)
abline(a = 0, b = slope_ci[1], col = "#6BAED6", lty = "dashed")  
abline(a = 0, b = slope_ci[2], col = "#6BAED6", lty = "dashed")  

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
outlier_points$image_x <- c(-1.7, -3.28, -4, 4, 2, 6, 1, -2.3, -1.2, 1.7)
outlier_points$image_y <- c(2.5, -2.8, 2, 3.2, 5.3, 4.4, -5.1, -5.2, -3.2, 2.8)
outlier_points$names <- rownames(outlier_points)
placePolygon <- function(x_init, y_init, x_poly, y_poly, pointName, side=1){
  x_init <- as.numeric(x_init)
  y_init <- as.numeric(y_init)
  x_poly <- as.numeric(x_poly)
  y_poly <- as.numeric(y_poly)
  left <- x_poly - side / 2
  bottom <- y_poly - side / 2
  right <- x_poly + side / 2
  top <- y_poly + side / 2
  segments(x_init, y_init, x_poly, y_poly, col = "black", lty = 2)
  polygon(c(left, left, right, right), y = c(top, bottom, bottom, top), col = "white")
  text(x_poly, top + 0.05 * (top - bottom), labels = pointName, pos = 3)
}
apply(outlier_points, 1, function(x) placePolygon(x[1], x[2], x[3], x[4], x[5], 1))
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

pdf("plots/turns-obs-boxplot.pdf", height=5, width=8)
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
# some labeling
# n_tips <- unlist(lapply(all_model_list, function(x) Ntip(x[[1]]$phy)))
clade_names <- unlist(lapply(strsplit(to_load, "results_"), function(x) gsub(".RData", "", x[2])))
clade_names <- unlist(lapply(strsplit(clade_names, "-"), function(x) x[[1]]))
clade_names <- unlist(lapply(strsplit(clade_names, "_"), function(x) x[[1]]))
# clade_names <- paste0(clade_names, " (", n_tips, ")")
clade_names_res <- clade_names
names(all_model_list) <- clade_names_res

to_load <- dir("6_recons/", full.names = TRUE)
to_load <- to_load[-grep("Quercus", to_load)]
to_load <- to_load[-grep("Araceae", to_load)]

all_recon_list <- list()
for(i in 1:length(to_load)){
  load(to_load[i])
  all_recon_list[[i]] <- recon
}
# make sure recons and results are in the same order
clade_names <- unlist(lapply(strsplit(to_load, "recon_"), function(x) gsub(".RData", "", x[2])))
clade_names <- unlist(lapply(strsplit(clade_names, "-"), function(x) x[[1]]))
clade_names <- unlist(lapply(strsplit(clade_names, "_"), function(x) x[[1]]))
clade_names_recon <- clade_names
all(match(clade_names_recon, clade_names_res) == 1:49)

# the calculation
aic_table <- do.call(rbind, lapply(all_model_list, GetAICWeights))
rownames(aic_table) <- clade_names
head(aic_table)

#dev.off()
#i = 1
# get model average recon
all_mod_avg_recon <- list()
for(i in 1:length(all_recon_list)){
  focal_recon_list <- all_recon_list[[i]]
  focal_aic_vec <- aic_table[i,][!unlist(lapply(focal_recon_list, is.null))]
  focal_aic_vec <- focal_aic_vec/sum(focal_aic_vec)
  focal_recon_list <- focal_recon_list[!unlist(lapply(focal_recon_list, is.null))]
  focal_mod_avg_recon <- matrix(0, dim(focal_recon_list[[1]]$node.mat)[1], 4)
  for(j in 1:length(focal_aic_vec)){
    tmp <- (focal_recon_list[[j]]$node.mat[,2:5] + 
      focal_recon_list[[j]]$node.mat[,6:9]) *
      focal_aic_vec[j]
    focal_mod_avg_recon <-  focal_mod_avg_recon + tmp
  }
  focal_mod_avg_recon <- focal_mod_avg_recon[,-3]
  colnames(focal_mod_avg_recon) <- c("closed", "widespread", "open")
  all_mod_avg_recon[[i]] <- list(phy = ladderize(focal_recon_list[[1]]$phy),
                                 recon = focal_mod_avg_recon)
  write.csv(focal_mod_avg_recon, row.names = FALSE,
            file = paste0("tables/", clade_names_recon[i], "-ASR.csv"))
  
}

##############################
### CLADE TABLE
##############################

names(all_mod_avg_recon) <- clade_names_recon
clade_table <- read.csv("clade_size_age_major_group.csv")
clade_table <- clade_table[clade_table$sp != "Quercus",]
clade_table <- clade_table[clade_table$sp != "Araceae",]
rownames(clade_table) <- clade_table$sp
clade_table <- clade_table[clade_names_recon,]
ages <- unlist(lapply(all_mod_avg_recon, function(x) max(branching.times(x$phy))))
ntips <- unlist(lapply(all_mod_avg_recon, function(x) Ntip(x$phy)))
brew_col <- brewer.pal(9, "Set1")[c(3,6,5)]
cols <- setNames(brew_col, c("closed", "widespread", "open"))

##############################
### plotting ancestral states for 51!!!!! clades...
##############################

##### set focal clade
focal_clade <- "monocot"
clade_index <- clade_table[,1] == focal_clade
names(clade_index) <- clade_names_recon
focal_ages <- sort(ages[clade_index], TRUE)
heights <- ntips[clade_index]/sum(ntips[clade_index]) + .1
heights <- heights[names(focal_ages)]
nclades <- length(which(clade_index))
pdf(file = paste0("plots/", focal_clade,"-asr-plot.pdf"), height = 20, width = 5)
layout(matrix(1:nclades, ncol = 1), heights = heights)
par(mar=c(.1,.1,.1,.1))
for(i in 1:nclades){
  ii <- which(names(focal_ages)[i] == names(clade_index))
  ntip <- Ntip(all_mod_avg_recon[[ii]]$phy)
  nnode <- Nnode(all_mod_avg_recon[[ii]]$phy)
  H <- max(branching.times(all_mod_avg_recon[[ii]]$phy))
  plot.phylo(all_mod_avg_recon[[ii]]$phy, 
             edge.color = "grey",
             x.lim = c(0, max(ages)), 
             plot = FALSE)
  box()
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  phylogram.plot(lastPP$edge, lastPP$Ntip, lastPP$Nnode, max(ages) - H + lastPP$xx, lastPP$yy, horizontal = TRUE, edge.color = "grey")
  node_labels_custom(pie = rep(1, ntip),
                     cex = .125,
                     piecol = "#FFFFFFFF",
                     offset = max(ages) - H)
  tmp_recon <- all_mod_avg_recon[[ii]]$recon + 1e-8 # a trick for polgon display
  tmp_recon <- tmp_recon/rowSums(tmp_recon)
  node_labels_custom(pie = tmp_recon,
             cex = .1,
             piecol = cols,
             border = NA,
             offset = max(ages) - H)
  node_labels_custom(node = ntip+1,
                     cex = .2,
                     pie = tmp_recon[1,],
                     piecol = cols,
                     offset = max(ages) - H)
  focal_dat <- all_model_list[[ii]][[1]]$data
  rownames(focal_dat) <- focal_dat[,1]
  focal_dat <- focal_dat[all_mod_avg_recon[[ii]]$phy$tip.label,]
  for(j in 1:ntip){
    tip_col <- cols[sum(focal_dat[j,c(2,3)]) + 1]
    segments(x0=max(ages) + 1, x1=max(ages) + 2, y0=j, y1=j, lwd = 1, col = tip_col)
  }
  text(x = 0, y = 2, labels = clade_names_recon[ii], adj=0)
}
dev.off()
# legend("bottomleft", legend=colnames(tmp), pch = 21, pt.bg = c("yellow", "orange", "red"), bty="n")

##############################
### ONE BIG ancestral states plot!!
##############################
phy_bb <- read.tree("backbone_tree.tre")
dat_names <- unlist(lapply(strsplit(paste(plot_data$Group), " "), function(x) x[[1]]))
phy_bb$tip.label <- unlist(lapply(strsplit(phy_bb$tip.label, "-"), function(x) x[[1]]))
phy_bb$tip.label <- unlist(lapply(strsplit(phy_bb$tip.label, "_"), function(x) x[[1]]))
to_drop <- phy_bb$tip.label[!phy_bb$tip.label %in% dat_names]
phy_bb <- drop.tip(phy_bb, to_drop)
phy_bb <- ladderize(phy_bb,right = F)
tip_names <- phy_bb$tip.label
phy_bb <- force.ultrametric(phy_bb)

names(all_mod_avg_recon) <- clade_names_recon
# clade_table <- read.csv("clade_size_age_major_group.csv")
# rownames(clade_table) <- clade_table$sp
# all(clade_names_recon %in% phy_bb$tip.label)
# all(phy_bb$tip.label %in% clade_names_recon)
clade_names_sorted <- phy_bb$tip.label
data.frame(no=1:49, sp=clade_names_sorted)
#some meta data
clade_table <- clade_table[clade_names_recon,]
ages <- unlist(lapply(all_mod_avg_recon, function(x) max(branching.times(x$phy))))
ntips <- unlist(lapply(all_mod_avg_recon, function(x) Ntip(x$phy)))
brew_col <- brewer.pal(9, "Set1")[c(3,6,5)]
cols <- setNames(brew_col, c("closed", "widespread", "open"))
# bb_ages <- branching.times(phy_bb)
# clade_inits <- c()
# for (i in seq_along(clade_names_recon)) {
#   tip_label <- clade_names_recon[i]
#   tip_index <- which(phy_bb$tip.label == tip_label)
#   node_label <- phy_bb$edge[tip_index,1]
#   init_age <- bb_ages[grep(node_label, names(bb_ages))]
#   clade_inits <- c(clade_inits, init_age)
# }
# names(clade_inits) <- clade_names_recon
# clade_inits <- max(bb_ages) - clade_inits

for (i in seq_along(clade_names_recon)) {
  tip_label <- clade_names_recon[i]
  tip_index <- which(phy_bb$tip.label == tip_label)
  subtree <- all_mod_avg_recon[[tip_label]]$phy
  subtree$node.label <- apply(all_mod_avg_recon[[tip_label]]$recon, 1, which.max)
  subtree <- ladderize(subtree, right=FALSE)
  edge_len  <- (phy_bb$edge.length[phy_bb$edge[,2] == tip_index]) - ages[i]
  if(edge_len < 0){
    edge_len <- (phy_bb$edge.length[phy_bb$edge[,2] == tip_index])
    subtree$edge.length <- subtree$edge.length/ages[i]*(edge_len-1)
    edge_len <- 1
  }
  phy_bb$edge.length[phy_bb$edge[,2] == tip_index] <- edge_len
  phy_bb <- bind.tree(phy_bb, subtree, where = tip_index)
}

C <- vcv(phy_bb)
dC <- diag(C)
H <- max(branching.times(phy_bb))
# make everything ultrametric
for(i in 1:Ntip(phy_bb)){
  to_add <- H - dC[i]
  index <- phy_bb$edge[,2] == i
  phy_bb$edge.length[index] <- phy_bb$edge.length[index] + to_add
}

bb_node_index <- Ntip(phy_bb) + which(is.na(phy_bb$node.label))
edges_to_shade <- match(phy_bb$edge[,1], bb_node_index)
edges_to_shade <- !is.na(edges_to_shade)
edge_lty <- rep(1, dim(phy_bb$edge)[1])
edge_lty[edges_to_shade] <- 3
edge_col <- rep("black", dim(phy_bb$edge)[1])
for(i in 1:nrow(phy_bb$edge)){
  focal_edge <- phy_bb$edge[i,]
  focal_state <- phy_bb$node.label[focal_edge[1] - Ntip(phy_bb)]
  if(!is.na(focal_state)){
    edge_col[i] <- cols[focal_state]
  }
}

# # # # # # # #actualy plotting
pdf(file="plots/all_clade_asr.pdf", width = 10, height = 10)

H <- 150
root_length <- H - max(branching.times(phy_bb))
phy_bb$root.edge <- root_length
plot.phylo_custom(x = phy_bb, 
                  show.tip.label = FALSE, 
                  type = "fan", 
                  edge.width = 0.25,
                  edge.lty =  edge_lty, 
                  edge.color = c(edge_col), 
                  no.margin = TRUE,
                  open.angle = 15, 
                  rotate.tree = 8.5,
                  root.edge = TRUE,
                  plot=FALSE,
                  x.lim=c(-165, 165),
                  y.lim=c(-165, 165))

# Define shades of grey
greys <- brewer.pal(9, "Greys")[1:6]

# Draw concentric filled circles (segments)
# geo_segments <- setNames(c(145, 100, 66, 23, 2.6), 
#                          nm = c("jurassic", "early cretaceous", "late cretaceous", "paleogene", "neogene"))
geo_segments <- setNames(c(145, 66, 23, 2.6, 0), 
                         nm = c("jurassic", "cretaceous", "paleogene", "neogene", "quaternary"))
geo_segments <- 150 - rev(geo_segments)

draw_circle(geo_segments[1], greys[5])
draw_circle(geo_segments[2], greys[4])
draw_circle(geo_segments[3], greys[3])
draw_circle(geo_segments[4], greys[2])
draw_circle(geo_segments[5], greys[1])

primary_ticks <- H - seq(from = 0, to = H, by = 50)
secondary_ticks <- H -seq(from = 0, to = H, by = 10)
tertiary_ticks <- H -seq(from = 0, to = H, by = 1)
labels <- seq(from = 0, to = H, by = 50)

segments(y0=-1, x0 = 0,
         y1=-1, x1 = max(tertiary_ticks))

segments(y0=-1, x0 = primary_ticks,
         y1=-5, x1 = primary_ticks)

segments(y0=-1, x0 = secondary_ticks,
         y1=-3.5, x1 = secondary_ticks)

segments(y0=-1, x0 = tertiary_ticks,
         y1=-2, x1 = tertiary_ticks)

text(y = -9, x = primary_ticks, labels = labels, cex = 0.5)
tmp <- data.frame(init = geo_segments[-c(1,5)], 
                  end = c(geo_segments[3:4],0))
x_pos <- rowMeans(tmp)
text(y = 2, x = x_pos, labels = names(x_pos), cex = 0.6)

# add polygons
pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip_coord <- data.frame(x = pp$xx[1:Ntip(phy_bb)], 
                        y = pp$yy[1:Ntip(phy_bb)])
# Generate a palette
palette <- brewer.pal(12, "Set3")
colors <- rep(palette, length.out = 51)
# colors <- viridis::rocket(51)
for (i in seq_along(clade_names_sorted)) {
  Sys.sleep(.2)
  tip_label <- clade_names_sorted[i]
  subtree <- all_mod_avg_recon[[tip_label]]$phy
  sub_tips <- subtree$tip.label[subtree$edge[,2][subtree$edge[,2] <= Ntip(subtree)]]
  focal_tips <- match(sub_tips, phy_bb$tip.label)
  tip_states <- rowSums(all_model_list[[tip_label]][[1]]$data[,c(2,3)]) + 1
  names(tip_states) <- all_model_list[[tip_label]][[1]]$data[,1]
  # if(length(tip_states)!=length(focal_tips)){
  #   print(tip_label)
  #   next
  # }
  tip_coords <- getPolygonCoords(tip_coord, 1.005, 1.02, focal_tips, TRUE)
  tip_coords$col <- cols[tip_states]
  poly_coords <- getPolygonCoords(tip_coord, 1.0275, 1.05, focal_tips, FALSE, i!=29)
  num_coords <- find_centroid(poly_coords[,1], poly_coords[,2])
  num_coords <- num_coords * 161/sqrt(sum(num_coords^2))
  for(j in 1:nrow(tip_coords)){
    tmp <- tip_coords[j,]
    segments(tmp[1,1], tmp[1,2], tmp[1,3], tmp[1,4], tmp[1,5], lwd=0.25)
  }
  polygon(poly_coords, col = colors[i], border = 'black')
  points(num_coords[1], num_coords[2], pch = 21, bg = colors[i], cex = 1.75)
  text(num_coords[1], num_coords[2], labels = i, cex = 0.5, adj = 0.5)
}

plot.phylo_custom(x = phy_bb,
                  show.tip.label = FALSE,
                  type = "fan",
                  edge.width = 0.25,
                  edge.lty =  edge_lty,
                  edge.color = c(edge_col),
                  no.margin = TRUE,
                  open.angle = 15,
                  root.edge = TRUE,
                  rotate.tree = 8.5,
                  add=TRUE)
dev.off()
##############################
### ancestral states by time plot!!
##############################

pdf("plots/ancestral_states_through_time.pdf", height=4, width=9)
bins <- seq(from=0, to=130, by=5)
recon_by_time <- list()
for(i in seq_along(clade_names_recon)){
  focal_tips <- rowSums(all_model_list[[i]][[1]]$data[,c(2,3)]) + 1
  names(focal_tips) <- all_model_list[[i]][[1]]$data[,1]
  focal_recon <- all_mod_avg_recon[[i]]$recon
  focal_subtree <- all_mod_avg_recon[[i]]$phy
  focal_times <- branching.times(focal_subtree)
  # remove root state
  focal_recon <- focal_recon[-1,]
  focal_times <- focal_times[-1]
  out <- c()
  for(j in 2:length(bins)){
    focal_bin <- bins[c(j,j-1)]
    index <- focal_times > focal_bin[2] & focal_bin[1] < focal_times
    focal_states <- colSums(matrix(focal_recon[index,], ncol = 3))
    out <- rbind(out, focal_states)
  }
  curr <- c(length(which(focal_tips == 1)),
            length(which(focal_tips == 2)),
            length(which(focal_tips == 3)))
  out <- rbind(curr, out)
  rownames(out) <- bins
  recon_by_time[[i]] <- out
}
names(recon_by_time) <- clade_names_recon

brew_col <- brewer.pal(9, "Set1")[c(3,6,5)]
cols <- setNames(brew_col, c("closed", "widespread", "open"))

dat_abs <- Reduce("+", recon_by_time)
dat_prop <- t(dat_abs/rowSums(dat_abs))
dat_prop <- dat_prop[,!is.nan(colSums(dat_prop))]
dat_prop <- dat_prop[,ncol(dat_prop):1]

dat_scaled <- t(log10(dat_abs+1))/4
dat_scaled <- dat_scaled[,colSums(dat_scaled)!=0]
dat_scaled <- dat_scaled[,ncol(dat_scaled):1]

###### plotting

#dev.off()

cols_mod <- make_less_vibrant(cols, 0.75, 1)
par(mar=c(0,0,.5,0))
barplot(dat_prop, 
        col=cols_mod, 
        axes=FALSE, ann=FALSE, axisnames=FALSE, 
        space = 0,
        ylim = c(-.3, 1.1), 
        xlim = c(-3.5, ncol(dat_prop)+3.5))

points(x = 1:ncol(dat_scaled)-.5, 
       y = dat_scaled[1,], 
       pch = 16, col = "white", cex = 2)

points(x = 1:ncol(dat_scaled)-.5, 
       y = dat_scaled[2,], 
       pch = 16, col = "white", cex = 2)

points(x = 1:ncol(dat_scaled)-.5, 
       y = dat_scaled[3,], 
       pch = 16, col = "white", cex = 2)

points(x = 1:ncol(dat_scaled)-.5, 
       y = dat_scaled[1,], 
       pch = 21, bg = cols[1], cex=1.5)

points(x = 1:ncol(dat_scaled)-.5, 
       y = dat_scaled[2,], 
       pch = 21, bg = cols[2], cex=1.5)

points(x = 1:ncol(dat_scaled)-.5, 
       y = dat_scaled[3,], 
       pch = 21, bg = cols[3], cex=1.5)


# axes lines
lines(x = c(0, ncol(dat_prop)), y = c(-.05,-.05), lwd=2)
# custom x
ticks_x <- 0:ncol(dat_prop)
segments(x0 = ticks_x, y0 = -0.05, x1 = ticks_x, y1 = -0.06) 
text(x = ticks_x, y = -0.075, 
     labels = c(c("110"), colnames(dat_prop)), 
     srt = 45, adj = 1, cex=0.7)
text(x = mean(ticks_x), y = -.15, labels = "Age (MY)", 
     srt = 0, adj = 0.5, xpd = TRUE)

# custom y
# left side
ticks_y <- c(0, 2.5, 5, 7.5, 10)/10
segments(x0 = -1, y0 = min(ticks_y), 
         x1 = -1, y1 = max(ticks_y))
segments(x0 = -1, y0 = ticks_y, 
         x1 = -1.1, y1 = ticks_y)
text(x = -1.15, y = ticks_y, labels = ticks_y, 
     srt = 0, adj = 1, xpd = TRUE)
text(x = -3, y = mean(ticks_y), labels = "Proportion", 
     srt = 90, adj = 0.5, xpd = TRUE)

# right side
ticks_y <- seq(from=0, to=10, length.out=5)/10
segments(x0 = ncol(dat_prop)+1, y0 = min(ticks_y), 
         x1 = ncol(dat_prop)+1, y1 = max(ticks_y))
segments(x0 = ncol(dat_prop)+1, y0 = ticks_y, 
         x1 = ncol(dat_prop)+1.1, y1 = ticks_y)
labels <- c(0,10,100,1000,10000)
text(x = ncol(dat_prop)+1.15, y = ticks_y, labels = labels, 
     srt = 0, adj = 0, xpd = TRUE)
text(x = ncol(dat_prop)+3.4, y = mean(ticks_y), labels = "Number of \nLineages", 
     srt = 90, adj = 0.5, xpd = TRUE)

dev.off()
# states <- phy_bb$node.label
# 
# times <- times[!is.na(states)]
# states <- states[!is.na(states)]
# 
# bins <- seq(from=0, to=130, by=2.5)
# out <- c()
# for(i in 2:length(bins)){
#   focal_bin <- bins[c(i,i-1)]
#   index <- times > focal_bin[2] & focal_bin[1] < times
#   focal_states <- states[index]
#   curr <- c(length(which(focal_states == 1)),
#             length(which(focal_states == 2)),
#             length(which(focal_states == 3)))
#   names(curr) <- c(1,2,3)
#   out <- rbind(out, curr)
# }
# rownames(out) <- bins
# 
# scaled_out <- apply(out, 1, function(x) x/sum(x))
# out[out == 0] <- NA
# out <- out + 1e-10
# par(mfrow=c(1,2))
# barplot(scaled_out, col=cols, horiz = T)
# barplot(log(t(out)) + 1, col=cols, horiz = T)


