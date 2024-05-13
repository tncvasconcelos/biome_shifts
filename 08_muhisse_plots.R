# rm(list=ls())
setwd("~/biome_shifts/")
source("00_utility_functions.R")

library(ape)
library(hisse)
library(parallel)
library(phytools)
library(phylolm)
library(phytools)

##############################
### load data
##############################

plot_data <- read.csv("all_par_table.csv")

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
tip_names <- phy_bb$tip.label


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

lm_dat <- lm_dat[phy_bb$tip.label,]
bar_dat <- bar_dat[phy_bb$tip.label,]

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
plot.new()
# Set up the full plot area
par(fig=c(0, 1, 0, 1), mar=c(.1,.1,.1,.1), new=TRUE)
plot(type="n", x = lm_dat$d_trans, y = lm_dat$d_turns, bty='n',
     xaxt='n', yaxt='n', xlab = "", ylab = "", 
     xlim=range(lm_dat$d_trans)+c(0, max(lm_dat$d_trans)*.2))
grid()
# axes lines
lines(x = round(c(min(lm_dat$d_trans), max(lm_dat$d_trans))), y = c(0, 0))
lines(x = c(0, 0), y = round(c(min(lm_dat$d_turns), max(lm_dat$d_turns))))

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

abline(a = 0, b = slope, col = "darkgrey", lwd = 2)
abline(a = 0, b = slope_ci[1], col = "lightblue", lty = "dashed")  
abline(a = 0, b = slope_ci[2], col = "lightblue", lty = "dashed")  

# add points
points(x = lm_dat$d_trans, y = lm_dat$d_turns,
       pch = 21, bg = rgb(1, 0, 0, 0.5), cex=1.5)

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
par(fig=c(0.495, 0.99, 0.05, 0.38), new=TRUE, mar=c(.1, .1, .1, .1))
# dev.off()
boxplot(bar_dat, main="", xlab="", ylab="", axes=FALSE, outline = FALSE,
        xlim=c(0, dim(bar_dat)[2]+1.5), ylim=range(bar_dat)*c(1.2,1),
        col=c("#b10026", "#fc9272", "#034e7b", "#74a9cf"))
box()
ticks_y <- round(seq(from = min(bar_dat), to = max(bar_dat), by = 2))
segments(x0 = 4.7, y0 = ticks_y, x1 = 4.9, y1 = ticks_y)
segments(x0 = 4.8, y0 = min(ticks_y), x1 = 4.8, y1 = max(ticks_y))
text(x = 5, y = ticks_y, labels = round(ticks_y, 2), 
     srt = 0, adj = 0, xpd = TRUE, cex=0.75)
text(x = 5.4, y = mean(range(bar_dat)), labels = expression(log(rate)), srt = 90, adj = 0.5, cex = 1)
text(x = 1, y = min(bar_dat)*1.18, labels = expression(q[A]))
text(x = 2, y = min(bar_dat)*1.18, labels = expression(q[B]))
text(x = 3, y = min(bar_dat)*1.18, labels = expression(tau[A]))
text(x = 4, y = min(bar_dat)*1.18, labels = expression(tau[B]))
points(x = rep(1, dim(bar_dat)[1]) + rnorm(dim(bar_dat)[1], sd = 0.05), 
       y = bar_dat$trans_a, pch = 21, bg="#b1002690")
points(x = rep(2, dim(bar_dat)[1]) + rnorm(dim(bar_dat)[1], sd = 0.05), 
       y = bar_dat$trans_b, pch = 21, bg="#fc927290")
points(x = rep(3, dim(bar_dat)[1]) + rnorm(dim(bar_dat)[1], sd = 0.05), 
       y = bar_dat$turns_a, pch = 21, bg="#034e7b90")
points(x = rep(4, dim(bar_dat)[1]) + rnorm(dim(bar_dat)[1], sd = 0.05), 
       y = bar_dat$turns_b, pch = 21, bg="#74a9cf90")

# Reset to full plot area
# par(fig=c(0, 1, 0, 1), new=TRUE)
# rasterImage(img, left, bottom, right, top)

##############################
### model avg ancestral states
##############################

# finished model sets for particular datsets
to_load <- dir("5_results/", full.names = TRUE)

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

to_load <- dir("6_recons/", full.names = TRUE)

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
all(match(clade_names_recon, clade_names_res) == 1:51)

# the calculation
aic_table <- do.call(rbind, lapply(all_model_list, GetAICWeights))
rownames(aic_table) <- clade_names
head(aic_table)

dev.off()
i = 1
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
### plotting ancestral states for 51!!!!! clades...
##############################

names(all_mod_avg_recon) <- clade_names_recon
clade_table <- read.csv("clade_size_age_major_group.csv")
rownames(clade_table) <- clade_table$sp
clade_table <- clade_table[clade_names_recon,]
ages <- unlist(lapply(all_mod_avg_recon, function(x) max(branching.times(x$phy))))
ntips <- unlist(lapply(all_mod_avg_recon, function(x) Ntip(x$phy)))
cols <- setNames(c('darkgreen', 'greenyellow', 'orange'), 
                 c("closed", "widespread", "open"))

##### set focal clade
focal_clade <- "asterids"
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
tip_names <- phy_bb$tip.label

names(all_mod_avg_recon) <- clade_names_recon
clade_table <- read.csv("clade_size_age_major_group.csv")
rownames(clade_table) <- clade_table$sp
all(clade_names_recon %in% phy_bb$tip.label)
all(phy_bb$tip.label %in% clade_names_recon)

#some meta data
clade_table <- clade_table[clade_names_recon,]
ages <- unlist(lapply(all_mod_avg_recon, function(x) max(branching.times(x$phy))))
ntips <- unlist(lapply(all_mod_avg_recon, function(x) Ntip(x$phy)))
cols <- setNames(c('darkgreen', 'greenyellow', 'orange'), 
                 c("closed", "widespread", "open"))
bb_ages <- branching.times(phy)
clade_inits <- c()
for (i in seq_along(clade_names_recon)) {
  tip_label <- clade_names_recon[i]
  tip_index <- which(phy_bb$tip.label == tip_label)
  node_label <- phy_bb$edge[tip_index,1]
  init_age <- bb_ages[grep(node_label, names(bb_ages))]
  clade_inits <- c(clade_inits, init_age)
}
names(clade_inits) <- clade_names_recon
clade_inits <- max(bb_ages) - clade_inits

for (i in seq_along(clade_names_recon)) {
  tip_label <- clade_names_recon[i]
  tip_index <- which(phy_bb$tip.label == tip_label)
  subtree <- all_mod_avg_recon[[tip_label]]$phy
  subtree$node.label <- apply(all_mod_avg_recon[[tip_label]]$recon, 1, which.max)
  edge_len  <- (phy_bb$edge.length[phy_bb$edge[,2] == tip_index]) - ages[i]
  phy_bb$edge.length[phy_bb$edge[,2] == tip_index] <- edge_len
  phy_bb <- bind.tree(phy_bb, subtree, where = tip_index)
}

bb_node_index <- Ntip(phy_bb) + which(is.na(phy_bb$node.label))
edges_to_shade <- match(phy_bb$edge[,1], bb_node_index)
edges_to_shade <- !is.na(edges_to_shade)
edge_lty <- rep(1, dim(phy_bb$edge)[1])
edge_lty[edges_to_shade] <- 2
plot.phylo(x = phy_bb, 
           show.tip.label = FALSE, 
           type = "fan", 
           edge.width = 0.25,
           edge.lty =  edge_lty, 
           edge.color = c("#444444", "darkgrey")[edge_lty], 
           no.margin = TRUE,
           open.angle = 10,
           rotate.tree = 275)
axisPhylo()
# abline(v = 0)
scale_init <- -max(branching.times(phy_bb)) + 150
scale_fin <- -max(branching.times(phy_bb))
segments(x0=0, y0 = -scale_init,
         x1=0, y1 = scale_fin)

nodelabels(pch = 16, 
           col = cols[phy_bb$node.label], 
           cex = 0.2)


which(is.na(phy_bb$node.label))