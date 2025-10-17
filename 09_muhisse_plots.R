# rm(list=ls())
setwd("~/biome_shifts/")
source("00_utility_functions.R")

library(ape)
library(hisse)

# finished model sets for particular datsets
to_load <- dir("6_recons/", full.names = TRUE)

recon_list <- list()
res_list <- list()

for(i in seq_along(to_load)){
  load(to_load[i])
  recon_list[[i]] <- recon
  file_name <- gsub("6_recons/", "5_results/", to_load[i])
  file_name <- gsub("recon", "res_state", file_name)
  load(file_name)
  res_list[[i]] <- res
}

clades <- sub(".*recon([^-_]+)[-_].*", "\\1", to_load)
names(res_list) <- names(recon_list) <- clades


recon_list$Antirrhineae

draw_full_ring <- function(r_in, r_out, col) {
  a <- seq(0, 2 * pi, length.out = 100)
  x_out <- r_out * cos(a)
  y_out <- r_out * sin(a)
  x_in <- r_in * cos(a)
  y_in <- r_in * sin(a)
  polygon(c(x_out, rev(x_in)), c(y_out, rev(y_in)), col = col, border = NA)
}
ED <- data.frame(
  S = c(2.6, 23.03, 66, 145),
  E = c(0, 2.6, 23.03, 66),
  C = c("grey100", "grey100", "grey100", "grey100"),
  stringsAsFactors = FALSE
)
CM <- c(
  "2" = "#009900", "3" = "#FFA500", "4" = "#FFFF00", 
  "6" = "#009900", "7" = "#FFA500", "8" = "#FFFF00"
)
h_list <- unlist(lapply(recon_list, function(x) max(branching.times(x$phy))))
par(mfrow=c(2,6), mar = c(0,0,1.5,0))
for(i in ord){
  phy <- ladderize(recon_list[[i]]$phy, right = FALSE)
  m_st <- apply(recon_list[[i]]$tip.mat[,-1], 1, which.max)
  m_st <- setNames(m_st, recon_list[[i]]$phy$tip.label[recon_list[[i]]$tip.mat[,1]])
  all_pars <- getPars(res_list[[i]])
  H <- h_list[i]
  plot.new()
  plot.window(xlim = c(-H * 1.05, H * 1.05), ylim = c(-H * 1.05, H * 1.05), asp = 1)
  for (j in rev(1:nrow(ED))) {
    T_s <- ED$S[j]
    T_e <- ED$E[j]
    R_in <- H - T_s
    R_out <- H - T_e
    R_d_in <- max(0, R_in)
    R_d_out <- min(H, R_out)
    if (R_d_out > R_d_in) {
      draw_full_ring(R_d_in, R_d_out, ED$C[j])
    }
  }
  focal_recon <- recon_list[[i]]$node.mat[,c(3:5,7:9)]
  turnover <- all_pars[grep("turnover", names(all_pars))]
  V <- colSums(t(focal_recon) * turnover)
  min_V <- min(V)
  max_V <- max(V)
  range_V <- max_V - min_V
  norm_V <- (V - min_V) / range_V
  # C_pal <- colorRampPalette(c("#4477AA", "#FFFFFF", "#BB3333"))
  C_pal <- colorRampPalette(c("grey20", "grey80"))
  C_all <- C_pal(100)[floor(norm_V * 99) + 1]
  plot.phylo_custom(phy,
    show.tip.label = FALSE, 
    type = "fan", 
    open.angle = 15, 
    add = TRUE,
    edge.color = C_all[phy$edge[,1] - Ntip(phy)]
  )
  # nodelabels(pie = focal_recon, piecol = CM, cex = 0.25)
  P <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  N <- length(phy$tip.label)
  X0 <- P$xx[1:N]
  Y0 <- P$yy[1:N]
  Tip_o <- phy$tip.label
  Angles <- atan2(Y0, X0)
  R_e <- H * 1.05 
  X1 <- R_e * cos(Angles)
  Y1 <- R_e * sin(Angles)
  R_e <- H * 1.025
  X0 <- R_e * cos(Angles)
  Y0 <- R_e * sin(Angles)
  
  q_rates <- all_pars[grep("q", names(all_pars))]
  q_rates <- c(A = mean(q_rates[1:4]), B = mean(q_rates[5:8]))
  focal_tips <- recon_list[[i]]$tip.mat[,c(3:5,7:9)]
  focal_tips <- data.frame(A = rowSums(focal_tips[,c(1:3)]), 
    B = rowSums(focal_tips[,c(4:6)]))
  V <- colSums(t(focal_tips) * q_rates)
  min_V <- min(V)
  max_V <- max(V)
  range_V <- max_V - min_V
  norm_V <- round(V - min_V, 10) / range_V
  C_pal <- colorRampPalette(c("#4477AA", "#BB3333"))
  T_cols <- C_pal(100)[floor(norm_V * 99) + 1]
  # T_cols <- CM[as.character(m_st[Tip_o])]
  segments(x0 = X0, y0 = Y0, x1 = X1, y1 = Y1, col = T_cols, lwd = 500/N)
  title(clades[i])
  
  R_e <- H * 1.015 
  X1 <- R_e * cos(Angles)
  Y1 <- R_e * sin(Angles)
  R_e <- H * 1.00
  X0 <- R_e * cos(Angles)
  Y0 <- R_e * sin(Angles)
  T_cols <- CM[as.character(m_st[Tip_o])]
  segments(x0 = X0, y0 = Y0, x1 = X1, y1 = Y1, col = T_cols, lwd = 500/N)
  
  Y_s <- -0.05 * H
  segments(x0 = 0, x1 = H, y0 = Y_s, y1 = Y_s, col = "black")
  Tick_10 <- seq(0, floor(H/10) * 10, by = 10)
  Tick_5 <- setdiff(seq(0, floor(H/5) * 5, by = 5), Tick_10)
  Tick_1 <- setdiff(seq(0, floor(H) * 1, by = 1), c(Tick_10, Tick_5))
  X_1 <- H - Tick_1
  X_5 <- H - Tick_5
  X_10 <- H - Tick_10
  segments(x0 = X_1, x1 = X_1, 
    y0 = Y_s - 0.002 * H, y1 = Y_s + 0.002 * H, 
    col = "black")
  segments(x0 = X_5, x1 = X_5, 
    y0 = Y_s - 0.004 * H, y1 = Y_s + 0.004 * H, 
    col = "black")
  segments(x0 = X_10, x1 = X_10, 
    y0 = Y_s - 0.007 * H, y1 = Y_s + 0.007 * H, 
    col = "black")
  text(x = X_10, y = Y_s, labels = Tick_10, 
    pos = 1, cex = 0.6, col = "black", srt = 45, offset = 0.5)
}

