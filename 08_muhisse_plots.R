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
library(ggsignif)
library(ggridges)

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
plot(phy_bb)
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
                     d_trans = log(rowMeans(rate_class_b_rate)) - log(rowMeans(rate_class_a_rate)),
                     d_turns = log(rowMeans(rate_class_b_turn)) - log(rowMeans(rate_class_a_turn)))

lm_dat <- lm_dat[match(phy_bb$tip.label, rownames(lm_dat)),]
fit = phylolm(d_turns ~ d_trans, data=lm_dat, phy=phy_bb, boot = 1000)
summary(fit)
lm_dat$predicted <- predict(fit, lm_dat)
lm_dat$ntaxa <-plot_data$ntip

a <- ggplot(lm_dat, aes(x = lm_dat$d_trans, y = lm_dat$d_turns)) +
  geom_point(aes(col = "red", alpha = 0.5)) +
  labs(
    x = expression(log(Delta * transition~rates)),
    y = expression(log(Delta * turnover~rates))) +
  theme_minimal() +
  ggtitle("a) Regression of differences") +
  geom_line(aes(x = lm_dat$d_trans, y = lm_dat$predicted), color = "blue") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme(legend.position = "none")
# geom_text(aes(x = lm_dat$d_trans, y = lm_dat$d_turns, label = rownames(lm_dat)), hjust = 0.4, vjust = 1, size = 3)

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

b <- b2 %>% insert_left(b1, width = 0.75)

c <- ggplot(filtered_data, aes(y = (value), x = name, fill=name)) +
  ggtitle("c) Distribution of rate differences") + 
  geom_violin() +
  theme_minimal()

bc <- grid.arrange(as.grob(b), as.grob(c), ncol = 2, widths = c(1,0.5))

abc <- grid.arrange(a, bc, ncol = 1, heights = c(0.5,1))
ggsave("plots/08_hidden_state_plot.pdf", abc)

##############################
### boxplot of observed state against turnover
##############################

# all clades
tn_00 <- cbind(plot_data$tn_00_a,
               plot_data$tn_00_b)
tn_01 <- cbind(plot_data$tn_01_a,
               plot_data$tn_01_b)
tn_11 <- cbind(plot_data$tn_11_a,
               plot_data$tn_11_b)

# phylo ttest
ttest_01_00 <- phyl.pairedttest(phy_bb, 
                                x1 = setNames(log(rowMeans(tn_01)), plot_data$Group),
                                x2 = setNames(log(rowMeans(tn_00)), plot_data$Group))
ttest_01_11 <- phyl.pairedttest(phy_bb, 
                                x1 = setNames(log(rowMeans(tn_01)), plot_data$Group),
                                x2 = setNames(log(rowMeans(tn_11)), plot_data$Group))
ttest_00_11 <- phyl.pairedttest(phy_bb, 
                                x1 = setNames(log(rowMeans(tn_00)), plot_data$Group),
                                x2 = setNames(log(rowMeans(tn_11)), plot_data$Group))

# Adding a group column to each dataset
tn_data <- data.frame(tn_00 = log(rowMeans(tn_00)),
                      tn_01 = log(rowMeans(tn_01)),
                      tn_11 = log(rowMeans(tn_11)))

long_data <- tn_data %>%
  pivot_longer(cols = everything(), names_to = "group", values_to = "values")

long_data$group <- factor(long_data$group, labels = c("C", "W", "O"))

max_value <- max(long_data$values, na.rm = TRUE)
buffer <- (max_value - min(long_data$values, na.rm = TRUE)) * 0.1  # 10% buffer

pairwise_comparisons <- data.frame(
  group1 = factor(c("W", "W", "C"),
                  levels = c("C", "W", "O")),
  group2 = factor(c("C", "O", "O"),
                  levels = c("C", "W", "O")),
  p.value = c(ttest_01_00$P.dbar, ttest_01_11$P.dbar, ttest_00_11$P.dbar),
  signif_label = ifelse(c(ttest_01_00$P.dbar, ttest_01_11$P.dbar, ttest_00_11$P.dbar) < 0.05, "*", ""))

pairwise_comparisons <- pairwise_comparisons %>%
  mutate(
    y_position = max_value + seq(1, n(), by = 1) * (buffer/2),  
    midpoint = ((as.numeric(group1) + as.numeric(group2)) / 2)
  )

ggplot() +
  geom_flat_violin(data = long_data, 
                   aes(x = group, y = values, fill = group, alpha = 0.5), 
                   position = position_nudge(x = -0.175)) +
  geom_point(data = long_data, 
             aes(x = group, y = values, color = group, size = 1),
             position = position_jitter(width = 0.05), alpha = 0.5) +
  geom_boxplot(data = long_data, 
               aes(x = group, y = values, fill = group, alpha = 0.5), 
               outlier.shape = NA, width = 0.25) +
  geom_segment(data = pairwise_comparisons,
               aes(x = group1, xend = group2, y = y_position, yend = y_position),
               linetype = "solid") +
  geom_text(data = pairwise_comparisons,
            aes(x = midpoint, y = y_position, label = signif_label), size = 10) +
  ggtitle("a)") +
  xlab("Observed State") +
  ylab(expression(log(Turnover))) +
  coord_cartesian(xlim = c(0.9, 3.1)) +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 18))

ggsave("plots/turnover_vs_obs_state.pdf")

##############################
### boxplot 
##############################

# all clades
f00_t01 <- cbind(plot_data$f00_t01_a,
                 plot_data$f00_t01_b)
f01_t00 <- cbind(plot_data$f01_t00_a,
                 plot_data$f01_t00_b)
f01_t11 <- cbind(plot_data$f01_t11_a,
                 plot_data$f01_t11_b)
f11_t01 <- cbind(plot_data$f11_t01_a,
                 plot_data$f11_t01_b)

# Adding a group column to each dataset
trans_data <- data.frame(f00_t01 = (rowMeans(f00_t01)),
                         f01_t00 = (rowMeans(f01_t00)),
                         f01_t11 = (rowMeans(f01_t11)),
                         f11_t01 = (rowMeans(f11_t01)))

ttest_a_b <- phyl.pairedttest(phy_bb, 
                              x1 = setNames(log(trans_data$f00_t01), plot_data$Group),
                              x2 = setNames(log(trans_data$f01_t00), plot_data$Group))
ttest_a_c <- phyl.pairedttest(phy_bb, 
                              x1 = setNames(log(trans_data$f00_t01), plot_data$Group),
                              x2 = setNames(log(trans_data$f01_t11), plot_data$Group))
ttest_a_d <- phyl.pairedttest(phy_bb, 
                              x1 = setNames(log(trans_data$f00_t01), plot_data$Group),
                              x2 = setNames(log(trans_data$f11_t01), plot_data$Group))

ttest_b_c <- phyl.pairedttest(phy_bb, 
                              x1 = setNames(log(trans_data$f01_t00), plot_data$Group),
                              x2 = setNames(log(trans_data$f01_t11), plot_data$Group))
ttest_b_d <- phyl.pairedttest(phy_bb, 
                              x1 = setNames(log(trans_data$f01_t00), plot_data$Group),
                              x2 = setNames(log(trans_data$f11_t01), plot_data$Group))

ttest_c_d <- phyl.pairedttest(phy_bb, 
                              x1 = setNames(log(trans_data$f01_t11), plot_data$Group),
                              x2 = setNames(log(trans_data$f11_t01), plot_data$Group))


long_data <- log(trans_data) %>%
  pivot_longer(cols = everything(), names_to = "group", values_to = "values")

long_data$group <- factor(long_data$group, labels = c("C to W", "W to C", "W to O", "O to W"))

max_value <- max(long_data$values, na.rm = TRUE)
buffer <- (max_value - min(long_data$values, na.rm = TRUE)) * 0.1  # 10% buffer

p_values <- c(ttest_a_b$P.dbar, ttest_a_c$P.dbar, ttest_a_d$P.dbar, ttest_b_c$P.dbar, ttest_b_d$P.dbar, ttest_c_d$P.dbar)
pairwise_comparisons <- data.frame(
  group1 = factor(
    c("C to W", "C to W", "C to W", "W to C", "W to C", "W to O"),
    levels = c("C to W", "W to C", "W to O", "O to W")),
  group2 = factor(
    c("W to C", "W to O", "O to W", "W to O", "O to W", "O to W"),
    levels = c("C to W", "W to C", "W to O", "O to W")),
  p.value = p_values,
  signif_label = ifelse(p_values < 0.05, "*", ""))

pairwise_comparisons <- pairwise_comparisons %>%
  mutate(
    y_position = max_value + seq(1, n(), by = 1) * (buffer/2),  
    midpoint = ((as.numeric(group1) + as.numeric(group2)) / 2)
  )

ggplot() +
  geom_flat_violin(data = long_data, 
                   aes(x = group, y = values, fill = group, alpha = 0.5), 
                   position = position_nudge(x = -0.175)) +
  geom_point(data = long_data, 
             aes(x = group, y = values, color = group, size = 1),
             position = position_jitter(width = 0.05), alpha = 0.5) +
  geom_boxplot(data = long_data, 
               aes(x = group, y = values, fill = group, alpha = 0.5), 
               outlier.shape = NA, width = 0.25) +
  geom_segment(data = pairwise_comparisons,
               aes(x = group1, xend = group2, y = y_position, yend = y_position),
               linetype = "solid") +
  geom_text(data = pairwise_comparisons,
            aes(x = midpoint, y = y_position, label = signif_label), size = 10) +
  ggtitle("b)") +
  xlab("Transition") +
  ylab(expression(log(Transition_Rate))) +
  coord_cartesian(xlim = c(1, 4)) +
  theme_minimal() +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 18))

ggsave("plots/trans_rate_vs_transition.pdf")

##############################
### examining extinction fraction dynamics
##############################

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
rate_class_a_ef <- cbind(plot_data$ef_00_a,
                         plot_data$ef_01_a,
                         plot_data$ef_11_a)
rate_class_b_ef <- cbind(plot_data$ef_00_b,
                         plot_data$ef_01_b,
                         plot_data$ef_11_b)

# scatter plot of diff
lm_dat <- data.frame(row.names = gsub(" .*", "", plot_data[,1]),
                     d_trans = log(rowMeans(rate_class_b_rate)) - log(rowMeans(rate_class_a_rate)),
                     d_ef = log(rowMeans(rate_class_b_ef)) - log(rowMeans(rate_class_a_ef)))

lm_dat <- lm_dat[match(phy_bb$tip.label, rownames(lm_dat)),]
fit = phylolm(d_ef ~ d_trans, data=lm_dat, phy=phy_bb, boot = 1000)
summary(fit)
lm_dat$predicted <- predict(fit, lm_dat)
lm_dat$ntaxa <- n_tips
lm_dat$aicwt <- result_matrix$both_differ

a <- ggplot(lm_dat, aes(x = lm_dat$d_trans, y = lm_dat$d_ef)) +
  geom_point(aes(col = "red", size = lm_dat$aicwt * 7, alpha = 0.5)) +
  labs(
    x = expression(log(Delta * transition~rates)),
    y = expression(log(Delta * ef~rates))) +
  theme_minimal() +
  ggtitle("a) Regression of differences") +
  geom_line(aes(x = lm_dat$d_trans, y = lm_dat$predicted), color = "blue") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)
# geom_text(aes(x = lm_dat$d_trans, y = lm_dat$d_turns, label = rownames(lm_dat)), hjust = 0.4, vjust = 1, size = 3)

b1 <- ggtree(phy_bb) +
  geom_tiplab(size = 2) +
  coord_cartesian(xlim = c(0, 180)) +
  ggtitle("b) Rate differences by clade")

lm_dat$id <- rownames(lm_dat)
lm_dat_long <- pivot_longer(lm_dat, cols = -id)
filtered_data <- lm_dat_long[lm_dat_long$name %in% c("d_trans", "d_ef"), ]

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

b <- b2 %>% insert_left(b1, width = 0.75)

c <- ggplot(filtered_data, aes(y = (value), x = name, fill=name)) +
  ggtitle("c) Distribution of rate differences") + 
  geom_violin() +
  theme_minimal()

bc <- grid.arrange(as.grob(b), as.grob(c), ncol = 2, widths = c(1,0.5))

abc <- grid.arrange(a, bc, ncol = 1)
ggsave("plots/08_hidden_state_plot.pdf", abc)


##############################
### regression of mean observed rate against turnover
##############################

# all clades
tn_00 <- cbind(plot_data$tn_00_a,
               plot_data$tn_00_b)
tn_01 <- cbind(plot_data$tn_01_a,
               plot_data$tn_01_b)
tn_11 <- cbind(plot_data$tn_11_a,
               plot_data$tn_11_b)

# phylo ttest
ttest_01_00 <- phyl.pairedttest(phy_bb, 
                                x1 = setNames(log(rowMeans(tn_01)), plot_data$Group),
                                x2 = setNames(log(rowMeans(tn_00)), plot_data$Group))
ttest_01_11 <- phyl.pairedttest(phy_bb, 
                                x1 = setNames(log(rowMeans(tn_01)), plot_data$Group),
                                x2 = setNames(log(rowMeans(tn_11)), plot_data$Group))
ttest_00_11 <- phyl.pairedttest(phy_bb, 
                                x1 = setNames(log(rowMeans(tn_00)), plot_data$Group),
                                x2 = setNames(log(rowMeans(tn_11)), plot_data$Group))

# Adding a group column to each dataset
tn_data <- data.frame(tn_00 = log(rowMeans(tn_00)),
                      tn_01 = log(rowMeans(tn_01)),
                      tn_11 = log(rowMeans(tn_11)))

long_data <- tn_data %>%
  pivot_longer(cols = everything(), names_to = "group", values_to = "values")

max_value <- max(long_data$values, na.rm = TRUE)
buffer <- (max_value - min(long_data$values, na.rm = TRUE)) * 0.1  # 10% buffer

pairwise_comparisons <- data.frame(
  group1 = factor(c('tn_00', 'tn_01', 'tn_00')),
  group2 = factor(c('tn_01', 'tn_11', 'tn_11')),
  p.value = c(ttest_01_00$P.dbar, ttest_01_11$P.dbar, ttest_00_11$P.dbar),
  signif_label = ifelse(c(ttest_01_00$P.dbar, ttest_01_11$P.dbar, ttest_00_11$P.dbar) < 0.05, "*", ""))

pairwise_comparisons <- pairwise_comparisons %>%
  mutate(
    y_position = max_value + seq(1, n(), by = 1) * buffer,  
    midpoint = ((as.numeric(group1) + as.numeric(group2)) / 2) + 0.5
  )

a <- ggplot(long_data, aes(x = group, y = values)) + 
  geom_boxplot() +
  geom_segment(data = pairwise_comparisons, 
               aes(x = group1, xend = group2, y = y_position, yend = y_position), 
               linetype = "solid") +
  geom_text(data = pairwise_comparisons, 
            aes(x = midpoint, y = y_position, label = signif_label)) +
  ggtitle("a) All clades included")

# hetero clades
plot_data <- plot_data[result_matrix$both_differ > 0.5,]
tn_00 <- cbind(plot_data$tn_00_a,
               plot_data$tn_00_b)
tn_01 <- cbind(plot_data$tn_01_a,
               plot_data$tn_01_b)
tn_11 <- cbind(plot_data$tn_11_a,
               plot_data$tn_11_b)

# phylo ttest
phy_bb <- keep.tip(phy_bb, plot_data$Group)
ttest_01_00 <- phyl.pairedttest(phy_bb, 
                                x1 = setNames(log(rowMeans(tn_01)), plot_data$Group),
                                x2 = setNames(log(rowMeans(tn_00)), plot_data$Group))
ttest_01_11 <- phyl.pairedttest(phy_bb, 
                                x1 = setNames(log(rowMeans(tn_01)), plot_data$Group),
                                x2 = setNames(log(rowMeans(tn_11)), plot_data$Group))
ttest_00_11 <- phyl.pairedttest(phy_bb, 
                                x1 = setNames(log(rowMeans(tn_00)), plot_data$Group),
                                x2 = setNames(log(rowMeans(tn_11)), plot_data$Group))

# Adding a group column to each dataset
tn_data <- data.frame(tn_00 = log(rowMeans(tn_00)),
                      tn_01 = log(rowMeans(tn_01)),
                      tn_11 = log(rowMeans(tn_11)))

long_data <- tn_data %>%
  pivot_longer(cols = everything(), names_to = "group", values_to = "values")

max_value <- max(long_data$values, na.rm = TRUE)
buffer <- (max_value - min(long_data$values, na.rm = TRUE)) * 0.1  # 10% buffer

pairwise_comparisons <- data.frame(
  group1 = factor(c('tn_00', 'tn_01', 'tn_00')),
  group2 = factor(c('tn_01', 'tn_11', 'tn_11')),
  p.value = c(ttest_01_00$P.dbar, ttest_01_11$P.dbar, ttest_00_11$P.dbar),
  signif_label = ifelse(c(ttest_01_00$P.dbar, ttest_01_11$P.dbar, ttest_00_11$P.dbar) < 0.05, "*", ""))

pairwise_comparisons <- pairwise_comparisons %>%
  mutate(
    y_position = max_value + seq(1, n(), by = 1) * buffer,  
    midpoint = ((as.numeric(group1) + as.numeric(group2)) / 2) + 0.5
  )

b <- ggplot(long_data, aes(x = group, y = values)) + 
  geom_boxplot() +
  geom_segment(data = pairwise_comparisons, 
               aes(x = group1, xend = group2, y = y_position, yend = y_position), 
               linetype = "solid") +
  geom_text(data = pairwise_comparisons, 
            aes(x = midpoint, y = y_position, label = signif_label)) +
  ggtitle("b) Hidden state support in turn and rates")

grid.arrange(a, b, nrow = 1)

