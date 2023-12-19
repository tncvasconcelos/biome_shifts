# rm(list=ls())
setwd("~/biome_shifts/")
source("00_utility_functions.R")

library(ape)
library(hisse)
library(parallel)
library(phytools)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggtree)
library(phylolm)
require(aplot)
library(gridExtra)
library(ggplotify)

##############################
### MODEL SUPPORT GENERALLY 
##############################
# finished model sets for particular datsets
to_load <- dir("5_results/", full.names = TRUE)

quick_check <- function(model_path){
  load(model_path)
  return(all(unlist(lapply(res, function(x) class(x) == "try-error"))))
}

failed_load <- to_load[sapply(to_load, quick_check)]
to_load <- to_load[!sapply(to_load, quick_check)]

all_model_list <- list()
for(i in 1:length(to_load)){
  load(to_load[i])
  all_model_list[[i]] <- res
}

# some labeling
n_tips <- unlist(lapply(all_model_list, function(x) Ntip(x[[1]]$phy)))
clade_names <- unlist(lapply(strsplit(to_load, "results_"), function(x) gsub(".RData", "", x[2])))
clade_names <- unlist(lapply(strsplit(clade_names, "-"), function(x) x[[1]]))
clade_names <- unlist(lapply(strsplit(clade_names, "_"), function(x) x[[1]]))
clade_names <- paste0(clade_names, " (", n_tips, ")")
# the calculation
aic_table <- do.call(rbind, lapply(all_model_list, GetAICWeights))
which_models <- do.call(rbind, lapply(all_model_list[[1]], test_hypotheses))
result_matrix <- aic_table %*% which_models
rownames(result_matrix) <- clade_names
head(result_matrix)

# Add row names as a column
result_matrix <- as.data.frame(result_matrix)
result_matrix <- result_matrix[,-2]
result_matrix$RowNames <- rownames(result_matrix)
# Reshape the data from wide to long format
result_matrix_long <- result_matrix %>%
  pivot_longer(cols = -RowNames, names_to = "Variable", values_to = "Value")

result_matrix_long$Variable <- factor(result_matrix_long$Variable, unique(result_matrix_long$Variable))
colnames(result_matrix_long)[1] <- "id"
# Create the heatmap
#D01B1B, #FF4242, #FFFFFF, #e7f9ff, #95D2EC and #47abd8. 
phy_bb <- read.tree("backbone_tree.tre")
phy_bb$tip.label <- gsub("-.*", "", phy_bb$tip.label)
phy_bb$tip.label <- gsub("_.*", "", phy_bb$tip.label)
phy_bb$tip.label <- clade_names[match(phy_bb$tip.label, gsub(" .*", "", clade_names))]
phy_bb$tip.label[is.na(match(phy_bb$tip.label, gsub(" .*", "", clade_names)))]
phy_bb <- drop.tip(phy_bb, which(is.na(phy_bb$tip.label)))

a <- ggtree(phy_bb) +
  geom_tiplab() +
  coord_cartesian(xlim = c(0, 180)) +
  ggtitle("a) Backbone Phylog10eny")

b <- ggplot(result_matrix_long, aes(x = Variable, y = id, fill = Value)) +
  ggtitle("b) Support for rate heterogeneity") +
  geom_tile() +
  scale_fill_viridis_c(option = "E", name = expression(sum(AIC[wt]))) +
  labs(x = NULL, y = NULL) +
  theme_tree2()

# ggplot(result_matrix_long, aes(x = Variable, y = id, fill = Value)) +
#   ggtitle("b) Support for rate heterogeneity") +
#   geom_tile() +
#   scale_fill_viridis_c(option = "E", name = expression(sum(AIC[wt]))) +
#   labs(x = NULL, y = NULL)

ab <- b %>% insert_left(a, width = 3)  
ab

# Print the plot
ggsave("plots/08_aic_support_plot.pdf", ab)

##############################
### TRANSITIONS VERSUS TURNOVER 
##############################

##############################
### BY RATE CLASS
##############################
# look at the ones where our hypothesis is potentially supportable and see whether rates are higher in one than the other. need to consider how to measure "rates"
both_differ_model_list <- all_model_list[result_matrix$both_differ > 0.5]
both_differ_clade_names <- clade_names[result_matrix$both_differ > 0.5]
both_differ_model_table <- lapply(both_differ_model_list, get_model_table)
plot_data_het <- cbind(Group = both_differ_clade_names, as.data.frame(do.call(rbind, lapply(both_differ_model_list, evaluate_difference))))
plot_data_het <- as.data.frame(plot_data_het)

# big clades
big_model_list <- all_model_list[n_tips > 100]
big_clade_names <- clade_names[n_tips > 100]
big_model_table <- lapply(big_model_list, get_model_table)
plot_data_big <- cbind(Group = big_clade_names, as.data.frame(do.call(rbind, lapply(big_model_list, evaluate_difference))))
plot_data_big <- as.data.frame(plot_data_big)


# plot data for all clades
plot_data_all <- cbind(Group = clade_names, as.data.frame(do.call(rbind, lapply(all_model_list, evaluate_difference))))
plot_data_all <- as.data.frame(plot_data_all)

##############################
### set which to analyze
##############################
library(gridExtra)
library(aplot)
library(ggplotify)
plot_data <- plot_data_big
plot_data_long <- pivot_longer(plot_data, cols = -Group)

# backbone phylog10eny
phy_bb <- read.tree("backbone_tree.tre")
dat_names <- unlist(lapply(strsplit(paste(plot_data$Group), " "), function(x) x[[1]]))
phy_bb$tip.label <- unlist(lapply(strsplit(phy_bb$tip.label, "-"), function(x) x[[1]]))
phy_bb$tip.label <- unlist(lapply(strsplit(phy_bb$tip.label, "_"), function(x) x[[1]]))
to_drop <- phy_bb$tip.label[!phy_bb$tip.label %in% dat_names]
phy_bb <- drop.tip(phy_bb, to_drop)
plot(phy_bb)
tip_names <- phy_bb$tip.label

# scatter plot of diff
lm_dat <- data.frame(row.names = gsub(" .*", "", plot_data[,1]),
                     d_trans = log10(plot_data$trans_b) - log10(plot_data$trans_a),
                     d_turns = log10(plot_data$turn_b) - log10(plot_data$turn_a))
lm_dat <- lm_dat[match(phy_bb$tip.label, rownames(lm_dat)),]
lm_dat$trans_00
fit = phylolm(d_turns ~ d_trans, data=lm_dat, phy=phy_bb, boot = 1000)
summary(fit)
lm_dat$predicted <- predict(fit, lm_dat)

a <- ggplot(lm_dat, aes(x = lm_dat$d_trans, y = lm_dat$d_turns)) +
  geom_point(aes(color = rownames(lm_dat)), size = 3) +
  labs(
    x = expression(log_10(Delta * transition~rates)),
    y = expression(log_10(Delta * turnover~rates))) +
  theme_minimal() +
  ggtitle("a) Regression of differences") +
  geom_line(aes(x = lm_dat$d_trans, y = lm_dat$predicted), color = "blue") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme(legend.position = "none") +
  geom_text(aes(x = lm_dat$d_trans, y = lm_dat$d_turns, label = rownames(lm_dat)), hjust = 0.4, vjust = 1, size = 3)

b1 <- ggtree(phy_bb) +
  geom_tiplab(size = 2) +
  coord_cartesian(xlim = c(0, 180)) +
  ggtitle("b) Rate differences by clade")

lm_dat$id <- rownames(lm_dat)
lm_dat_long <- pivot_longer(lm_dat, cols = -id)
filtered_data <- lm_dat_long[lm_dat_long$name %in% c("d_trans", "d_turns"), ]

b2 <- ggplot(filtered_data, aes(x=id, y=value, fill=name)) + 
  geom_bar(stat="identity", position="dodge") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),      # Remove axis titles
    axis.text.y = element_blank(),       # Remove axis text
    plot.title = element_blank(),      # Remove plot title
    legend.position = "none"
  ) +
  coord_flip()

# ggplot(result_matrix_long, aes(x = Variable, y = id, fill = Value)) +
#   ggtitle("b) Support for rate heterogeneity") +
#   geom_tile() +
#   scale_fill_viridis_c(option = "E", name = expression(sum(AIC[wt]))) +
#   labs(x = NULL, y = NULL)

b <- b2 %>% insert_left(b1, width = 0.75)

c <- ggplot(filtered_data, aes(y = (value), x = name, fill=name)) +
  ggtitle("c) Distribution of rate differences") + 
  geom_violin() +
  theme_minimal()

bc <- grid.arrange(as.grob(b), as.grob(c), ncol = 2, widths = c(1,0.5))

abc <- grid.arrange(a, bc, ncol = 1)
ggsave("plots/08_hidden_state_plot.pdf", abc)

##############################
### EXAMINING TURNOVER BY OBSERVED STATE
##############################

rownames(plot_data) <- gsub(" .*", "", plot_data$Group)
plot_data$Group <- gsub(" .*", "", plot_data$Group)
plot_data <- plot_data[match(phy_bb$tip.label, rownames(plot_data)),]

head(plot_data[,-1])
apply(plot_data[,-1], 2, median)

plot_data_long <- pivot_longer(plot_data, cols = -Group) 

fit_00 = phylolm(log10(turn_00) ~ log10(trans_00), data=plot_data, phy=phy_bb)
fit_01 = phylolm(log10(turn_01) ~ log10(trans_01), data=plot_data, phy=phy_bb)
fit_11 = phylolm(log10(turn_11) ~ log10(trans_11), data=plot_data, phy=phy_bb)

summary(fit_00)
summary(fit_01)
summary(fit_11)

subset_data <- plot_data_long[grep("turn", plot_data_long$name),]
subset_data$value <- log10(subset_data$value)

ind_01 <- subset_data$name == "turn_01"
ind_11 <- subset_data$name == "turn_11"
ind_00 <- subset_data$name == "turn_00"

ttest_01_00 <- phyl.pairedttest(phy_bb, 
                                x1 = setNames(subset_data$value[ind_01], subset_data$Group[ind_01]),
                                x2 = setNames(subset_data$value[ind_00], subset_data$Group[ind_00]))
ttest_01_11 <- phyl.pairedttest(phy_bb, 
                                x1 = setNames(subset_data$value[ind_01], subset_data$Group[ind_01]),
                                x2 = setNames(subset_data$value[ind_11], subset_data$Group[ind_11]))
ttest_00_11 <- phyl.pairedttest(phy_bb, 
                                x1 = setNames(subset_data$value[ind_00], subset_data$Group[ind_00]),
                                x2 = setNames(subset_data$value[ind_11], subset_data$Group[ind_11]))


a <- ggplot(subset_data, aes(x = name, y = log10(value), fill = name)) +
  geom_violin() +
  ggtitle("a) comparison of observed state turnover rates") +
  theme_minimal()

b <- ggplot(plot_data, aes(x = log10(trans_00), y = log10(turn_00))) +
  geom_point(aes(color = rownames(plot_data)), size = 3) +
  labs(x = expression(log10(transition~rates)),y = expression(log10(turnover~rates))) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_smooth(method = "lm") +
  ggtitle("b) regression of 00") +
  theme(legend.position = "none")

c <- ggplot(plot_data, aes(x = log10(trans_01), y = log10(turn_01))) +
  geom_point(aes(color = rownames(plot_data)), size = 3) +
  labs(x = expression(log10(transition~rates)),y = expression(log10(turnover~rates))) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_smooth(method = "lm") +
  ggtitle("c) regression of 01") +
  theme(legend.position = "none")

d <- ggplot(plot_data, aes(x = log10(trans_11), y = log10(turn_11))) +
  geom_point(aes(color = rownames(plot_data)), size = 3) +
  labs(x = expression(log10(transition~rates)),y = expression(log10(turnover~rates))) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_smooth(method = "lm") +
  ggtitle("d) regression of 11") +
  theme(legend.position = "none")

bcd <- grid.arrange(b,c,d, ncol = 1)
abcd <- grid.arrange(a, bcd, ncol = 2, widths = c(0.5, 1))

ggsave("plots/08_observed_state_plot.pdf", abcd)

lm_dat$trans_00 <- plot_data$trans_00
lm_dat$trans_01 <- plot_data$trans_01
lm_dat$trans_11 <- plot_data$trans_11

fit = phylolm(d_turns ~ d_trans + trans_00 + trans_01 + trans_11, data=lm_dat, phy=phy_bb)
summary(fit)

# lm_dat$predicted <- predict(fit, lm_dat)

# ggplot(lm_dat, aes(x = lm_dat$d_trans, y = lm_dat$d_turns)) +
#   geom_point(aes(color = rownames(lm_dat)), size = 3) +
#   labs(
#     x = expression(Delta * log10(transition~rates)),
#     y = expression(Delta * log10(turnover~rates))) +
#   theme_minimal() +
#   geom_line(aes(x = lm_dat$d_trans, y = lm_dat$predicted), color = "blue") +
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = 0) +
#   theme(legend.position = "none") +
#   geom_text(aes(x = lm_dat$d_trans, y = lm_dat$d_turns, label = rownames(lm_dat)), hjust = 0.4, vjust = 1, size = 3)

##############################
### EXAMINING SHIFTS OUT OF THE WIDESPREAD STATE
##############################

# plot data for heterogenous clades
plot_data_het <- cbind(Group = both_differ_clade_names, as.data.frame(do.call(rbind, lapply(both_differ_model_list, function(x) evaluate_difference(x, type = "trans")))))
plot_data_het <- as.data.frame(plot_data_het)

# plot data for all clades
plot_data_all <- cbind(Group = clade_names, as.data.frame(do.call(rbind, lapply(all_model_list, function(x) evaluate_difference(x, type = "trans")))))
plot_data_all <- as.data.frame(plot_data_all)

plot_data <- plot_data_all

# backbone phylog10eny
phy_bb <- read.tree("backbone_tree.tre")
dat_names <- unlist(lapply(strsplit(paste(plot_data$Group), " "), function(x) x[[1]]))
phy_bb$tip.label <- unlist(lapply(strsplit(phy_bb$tip.label, "-"), function(x) x[[1]]))
phy_bb$tip.label <- unlist(lapply(strsplit(phy_bb$tip.label, "_"), function(x) x[[1]]))
to_drop <- phy_bb$tip.label[!phy_bb$tip.label %in% dat_names]
phy_bb <- drop.tip(phy_bb, to_drop)
plot(phy_bb)
tip_names <- phy_bb$tip.label

boxplot(log10(plot_data[,-1]))

# scatter plot of diff

rownames(plot_data) <- gsub(" .*", "", plot_data$Group)
plot_data$Group <- gsub(" .*", "", plot_data$Group)
plot_data <- plot_data[match(phy_bb$tip.label, rownames(plot_data)),]

colMeans(log10(plot_data[,-1]))
plot_data_combined <- ((plot_data[,2:5]) + (plot_data[,6:9]))/2
plot_data_combined$Group <- rownames(plot_data_combined)
apply((plot_data_combined[,-5]), 2, median)

plot_data_long <- pivot_longer(plot_data_combined, cols = -Group) 

ggplot(plot_data_long, aes(x = name, y = log10(value), fill = name)) +
  geom_boxplot() +
  ggtitle("a) comparison of observed state transition rates") +
  theme_minimal() +
  coord_cartesian(ylim = c(-2, 2))
  

b <- ggplot(plot_data, aes(x = log10(trans_00), y = log10(turn_00))) +
  geom_point(aes(color = rownames(plot_data)), size = 3) +
  labs(x = expression(log10(transition~rates)),y = expression(log10(turnover~rates))) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_smooth(method = "lm") +
  ggtitle("b) regression of 00") +
  theme(legend.position = "none")

c <- ggplot(plot_data, aes(x = log10(trans_01), y = log10(turn_01))) +
  geom_point(aes(color = rownames(plot_data)), size = 3) +
  labs(x = expression(log10(transition~rates)),y = expression(log10(turnover~rates))) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_smooth(method = "lm") +
  ggtitle("c) regression of 01") +
  theme(legend.position = "none")

d <- ggplot(plot_data, aes(x = log10(trans_11), y = log10(turn_11))) +
  geom_point(aes(color = rownames(plot_data)), size = 3) +
  labs(x = expression(log10(transition~rates)),y = expression(log10(turnover~rates))) +
  theme_minimal() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_smooth(method = "lm") +
  ggtitle("d) regression of 11") +
  theme(legend.position = "none")

bcd <- grid.arrange(b,c,d, ncol = 1)
abcd <- grid.arrange(a, bcd, ncol = 2, widths = c(0.5, 1))

ggsave("plots/08_observed_state_plot.pdf", abcd)

lm_dat$trans_00 <- plot_data$trans_00
lm_dat$trans_01 <- plot_data$trans_01
lm_dat$trans_11 <- plot_data$trans_11

fit = phylolm(d_turns ~ d_trans + trans_00 + trans_01 + trans_11, data=lm_dat, phy=phy_bb)
summary(fit)

##############################
### TRANSITIONS VERSUS NETDIV 
##############################

##############################
### BY RATE CLASS
##############################
# look at the ones where our hypothesis is potentially supportable and see whether rates are higher in one than the other. need to consider how to measure "rates"
both_differ_model_list <- all_model_list[result_matrix$both_differ > 0.5]
both_differ_clade_names <- clade_names[result_matrix$both_differ > 0.5]
both_differ_model_table <- lapply(both_differ_model_list, get_model_table)

# plot data for heterogenous clades
plot_data_het <- cbind(Group = both_differ_clade_names, as.data.frame(do.call(rbind, lapply(both_differ_model_list, evaluate_difference))))
plot_data_het <- as.data.frame(plot_data_het)

# plot data for all clades
plot_data_all <- cbind(Group = clade_names, as.data.frame(do.call(rbind, lapply(all_model_list, evaluate_difference))))
plot_data_all <- as.data.frame(plot_data_all)

##############################
### set which to analyze
##############################
library(gridExtra)
library(aplot)
library(ggplotify)
plot_data <- plot_data_all
plot_data_long <- pivot_longer(plot_data, cols = -Group)

# backbone phylog10eny
phy_bb <- read.tree("backbone_tree.tre")
dat_names <- unlist(lapply(strsplit(paste(plot_data$Group), " "), function(x) x[[1]]))
phy_bb$tip.label <- unlist(lapply(strsplit(phy_bb$tip.label, "-"), function(x) x[[1]]))
phy_bb$tip.label <- unlist(lapply(strsplit(phy_bb$tip.label, "_"), function(x) x[[1]]))
to_drop <- phy_bb$tip.label[!phy_bb$tip.label %in% dat_names]
phy_bb <- drop.tip(phy_bb, to_drop)
plot(phy_bb)
tip_names <- phy_bb$tip.label

ggplot(plot_data_long, aes(x = log10(value), fill = name)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5, bins = 30) +
  geom_density(aes(color = name), size = 1, alpha = 0.5) +
  labs(title = "Overlaid Histogram", x = "Value", y = "Density") +
  theme_classic() +
  facet_wrap(~name)

# scatter plot of diff
lm_dat <- data.frame(row.names = gsub(" .*", "", plot_data[,1]),
                     d_trans = (plot_data$trans_b) - (plot_data$trans_a),
                     d_divs = (plot_data$net.div_b) - (plot_data$net.div_a))

lm_dat <- log(lm_dat - 2*min(lm_dat))

lm_dat <- lm_dat[match(phy_bb$tip.label, rownames(lm_dat)),]
lm_dat$trans_00
fit = phylolm(d_divs ~ d_trans, data=lm_dat, phy=phy_bb, boot = 1000)
summary(fit)
lm_dat$predicted <- predict(fit, lm_dat)

a <- ggplot(lm_dat, aes(x = lm_dat$d_trans, y = lm_dat$d_divs)) +
  geom_point(aes(color = rownames(lm_dat)), size = 3) +
  labs(
    x = expression(log_10(Delta * transition~rates)),
    y = expression(log_10(Delta * turnover~rates))) +
  theme_minimal() +
  ggtitle("a) Regression of differences") +
  theme(legend.position = "none") +
  geom_text(aes(x = lm_dat$d_trans, y = lm_dat$d_divs, label = rownames(lm_dat)), hjust = 0.4, vjust = 1, size = 3)
  

b1 <- ggtree(phy_bb) +
  geom_tiplab(size = 2) +
  coord_cartesian(xlim = c(0, 180)) +
  ggtitle("b) Rate differences by clade")

lm_dat$id <- rownames(lm_dat)
lm_dat_long <- pivot_longer(lm_dat, cols = -id)
filtered_data <- lm_dat_long[lm_dat_long$name %in% c("d_trans", "d_turns"), ]

b2 <- ggplot(filtered_data, aes(x=id, y=value, fill=name)) + 
  geom_bar(stat="identity", position="dodge") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),      # Remove axis titles
    axis.text.y = element_blank(),       # Remove axis text
    plot.title = element_blank(),      # Remove plot title
    legend.position = "none"
  ) +
  coord_flip()

# ggplot(result_matrix_long, aes(x = Variable, y = id, fill = Value)) +
#   ggtitle("b) Support for rate heterogeneity") +
#   geom_tile() +
#   scale_fill_viridis_c(option = "E", name = expression(sum(AIC[wt]))) +
#   labs(x = NULL, y = NULL)

b <- b2 %>% insert_left(b1, width = 0.75)

c <- ggplot(filtered_data, aes(y = (value), x = name, fill=name)) +
  ggtitle("c) Distribution of rate differences") + 
  geom_violin() +
  theme_minimal()

bc <- grid.arrange(as.grob(b), as.grob(c), ncol = 2, widths = c(1,0.5))

abc <- grid.arrange(a, bc, ncol = 1)

##############################
### BY OBSERVED STATE
##############################

# plot data for heterogenous clades
plot_data_het <- cbind(Group = both_differ_clade_names, as.data.frame(do.call(rbind, lapply(both_differ_model_list, function(x) evaluate_difference(x, type = "obs")))))
plot_data_het <- as.data.frame(plot_data_het)

# plot data for all clades
plot_data_all <- cbind(Group = clade_names, as.data.frame(do.call(rbind, lapply(all_model_list, function(x) evaluate_difference(x, type = "obs")))))
plot_data_all <- as.data.frame(plot_data_all)

##############################
### set which to analyze
##############################

plot_data <- plot_data_all

# backbone phylog10eny
phy_bb <- read.tree("backbone_tree.tre")
dat_names <- unlist(lapply(strsplit(paste(plot_data$Group), " "), function(x) x[[1]]))
phy_bb$tip.label <- unlist(lapply(strsplit(phy_bb$tip.label, "-"), function(x) x[[1]]))
phy_bb$tip.label <- unlist(lapply(strsplit(phy_bb$tip.label, "_"), function(x) x[[1]]))
to_drop <- phy_bb$tip.label[!phy_bb$tip.label %in% dat_names]
phy_bb <- drop.tip(phy_bb, to_drop)
plot(phy_bb)
tip_names <- phy_bb$tip.label

boxplot(log10(plot_data[,-1]))

# scatter plot of diff

rownames(plot_data) <- gsub(" .*", "", plot_data$Group)
plot_data$Group <- gsub(" .*", "", plot_data$Group)
plot_data <- plot_data[match(phy_bb$tip.label, rownames(plot_data)),]

plot_data_long <- pivot_longer(plot_data, cols = -Group) 

fit_00 = phylolm(log10(net.div_00) ~ log10(trans_00), data=plot_data, phy=phy_bb)
fit_01 = phylolm(log10(net.div_01) ~ log10(trans_01), data=plot_data, phy=phy_bb)
fit_11 = phylolm(log10(net.div_11) ~ log10(trans_11), data=plot_data, phy=phy_bb)

summary(fit_00)
summary(fit_01)
summary(fit_11)

subset_data <- plot_data_long[grep("net.div", plot_data_long$name),]

ggplot(subset_data, aes(x = name, y = (value), fill = name)) +
  geom_boxplot() +
  ggtitle("a) comparison of observed state diversification rates") +
  theme_minimal() + 
  coord_cartesian(ylim = c(-2,2))


##############################
### EXAMINING ALL parameters
##############################

# plot data for heterogenous clades
plot_data_het <- cbind(Group = both_differ_clade_names, as.data.frame(do.call(rbind, lapply(both_differ_model_list, function(x) evaluate_difference(x, type = "all")))))
plot_data_het <- as.data.frame(plot_data_het)

# plot data for all clades
plot_data_all <- cbind(Group = clade_names, as.data.frame(do.call(rbind, lapply(all_model_list, function(x) evaluate_difference(x, type = "all")))))
plot_data_all <- as.data.frame(plot_data_all)

plot_data <- plot_data_all
colnames(plot_data) <- gsub("\\..*", "", colnames(plot_data))
head(plot_data)
plot_data <- cbind(plot_data,result_matrix[,-5])

write.csv(plot_data, file = "all_par_table.csv", row.names = FALSE)
colMeans(plot_data[,-1])