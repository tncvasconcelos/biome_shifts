library(ape)
library(TreeSim)
setwd("~/biome_shifts/")

set.seed(1500)
phy_a <- sim.bd.taxa.age(25, 1, 1,  0.9, age = 10, mrca = TRUE)[[1]]
phy_b <- sim.bd.taxa.age(n = 15, numbsim = 1, lambda = 0.5,  mu = 0.45, age = 10, mrca = TRUE)[[1]]
phy_c <- sim.bd.taxa.age(n = 10, numbsim = 1, lambda = 0.1,  mu = 0.09, age = 10, mrca = TRUE)[[1]]
phy_a <- ladderize(phy_a)
phy_b <- ladderize(phy_b)
phy_c <- ladderize(phy_c)

tip_cols <- c("orange", "lightgreen", "darkgreen")
edge_cols <- c("#b10026", "#034e7b")
# model outcome one
pdf("plots/MO1_PhyA.pdf", onefile = FALSE, width = 3)
plot.phylo(phy_a, show.tip.label = FALSE, edge.color = edge_cols, edge.width = 5, x.lim = c(0, 11))
tiplabels(pch = 16, col = tip_cols[1], offset = 0.5, cex = 2)
dev.off()
pdf("plots/MO1_PhyB.pdf", onefile = FALSE, height = 5, width = 3)
plot(phy_b, show.tip.label = FALSE, edge.color = edge_cols, edge.width = 5, x.lim = c(0, 11))
tiplabels(pch = 16, col = tip_cols[2], offset = 0.5, cex = 2)
dev.off()
pdf("plots/MO1_PhyC.pdf", onefile = FALSE, height = 4, width = 3)
plot(phy_c, show.tip.label = FALSE, edge.color = edge_cols, edge.width = 5, x.lim = c(0, 11))
tiplabels(pch = 16, col = tip_cols[3], offset = 0.5, cex = 2)
dev.off()


# model outcome two
pdf("plots/MO2_PhyA.pdf", onefile = FALSE, width = 3)
plot.phylo(phy_a, show.tip.label = FALSE, edge.color = edge_cols[1], edge.width = 5, x.lim = c(0, 11))
tiplabels(pch = 16, col = tip_cols, offset = 0.5, cex = 2)
dev.off()
pdf("plots/MO2_PhyB.pdf", onefile = FALSE, height = 5, width = 3)
plot(phy_b, show.tip.label = FALSE, edge.color = edge_cols, edge.width = 5, x.lim = c(0, 11))
tiplabels(pch = 16, col = tip_cols, offset = 0.5, cex = 2)
dev.off()
pdf("plots/MO2_PhyC.pdf", onefile = FALSE, height = 4, width = 3)
plot(phy_c, show.tip.label = FALSE, edge.color = edge_cols[2], edge.width = 5, x.lim = c(0, 11))
tiplabels(pch = 16, col = tip_cols, offset = 0.5, cex = 2)
dev.off()

# model outcome three
pdf("plots/MO3_PhyA.pdf", onefile = FALSE, width = 3)
plot.phylo(phy_a, show.tip.label = FALSE, edge.color = edge_cols[2], edge.width = 5, x.lim = c(0, 11))
tiplabels(pch = 16, col = tip_cols, offset = 0.5, cex = 2)
dev.off()
pdf("plots/MO3_PhyB.pdf", onefile = FALSE, height = 5, width = 3)
plot(phy_b, show.tip.label = FALSE, edge.color = edge_cols, edge.width = 5, x.lim = c(0, 11))
tiplabels(pch = 16, col = tip_cols, offset = 0.5, cex = 2)
dev.off()
pdf("plots/MO3_PhyC.pdf", onefile = FALSE, height = 4, width = 3)
plot(phy_c, show.tip.label = FALSE, edge.color = edge_cols[1], edge.width = 5, x.lim = c(0, 11))
tiplabels(pch = 16, col = tip_cols, offset = 0.5, cex = 2)
dev.off()


