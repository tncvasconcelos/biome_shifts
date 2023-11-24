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
### Creating a master table
##############################
# look at the ones where our hypothesis is potentially supportable and see whether rates are higher in one than the other. need to consider how to measure "rates"
both_differ_model_list <- all_model_list[result_matrix$both_differ > 0.5]
both_differ_clade_names <- clade_names[result_matrix$both_differ > 0.5]
both_differ_model_table <- lapply(both_differ_model_list, get_model_table)

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
rownames(plot_data) <- gsub(" .*", "", plot_data$Group)
plot_data$Group <- gsub(" .*", "", plot_data$Group)

write.csv(plot_data, file = "all_par_table.csv", row.names = FALSE)
colMeans(plot_data[,-1])

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
lm_dat$ntaxa <- n_tips
lm_dat$aicwt <- result_matrix$both_differ

a <- ggplot(lm_dat, aes(x = lm_dat$d_trans, y = lm_dat$d_turns)) +
  geom_point(aes(col = "red", size = lm_dat$aicwt * 7, alpha = 0.5)) +
  labs(
    x = expression(log(Delta * transition~rates)),
    y = expression(log(Delta * turnover~rates))) +
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

abc <- grid.arrange(a, bc, ncol = 1)
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

ggplot(long_data, aes(x = group, y = values)) + 
  geom_boxplot() +
  geom_segment(data = pairwise_comparisons, 
               aes(x = group1, xend = group2, y = y_position, yend = y_position), 
               linetype = "solid") +
  geom_text(data = pairwise_comparisons, 
            aes(x = midpoint, y = y_position, label = signif_label)) +
  ggtitle("a) All clades included")

# hetero clades
# plot_data <- plot_data[result_matrix$both_differ > 0.5,]
# tn_00 <- cbind(plot_data$tn_00_a,
#                plot_data$tn_00_b)
# tn_01 <- cbind(plot_data$tn_01_a,
#                plot_data$tn_01_b)
# tn_11 <- cbind(plot_data$tn_11_a,
#                plot_data$tn_11_b)
# 
# # phylo ttest
# phy_bb <- keep.tip(phy_bb, plot_data$Group)
# ttest_01_00 <- phyl.pairedttest(phy_bb, 
#                                 x1 = setNames(log(rowMeans(tn_01)), plot_data$Group),
#                                 x2 = setNames(log(rowMeans(tn_00)), plot_data$Group))
# ttest_01_11 <- phyl.pairedttest(phy_bb, 
#                                 x1 = setNames(log(rowMeans(tn_01)), plot_data$Group),
#                                 x2 = setNames(log(rowMeans(tn_11)), plot_data$Group))
# ttest_00_11 <- phyl.pairedttest(phy_bb, 
#                                 x1 = setNames(log(rowMeans(tn_00)), plot_data$Group),
#                                 x2 = setNames(log(rowMeans(tn_11)), plot_data$Group))
# 
# # Adding a group column to each dataset
# tn_data <- data.frame(tn_00 = log(rowMeans(tn_00)),
#                       tn_01 = log(rowMeans(tn_01)),
#                       tn_11 = log(rowMeans(tn_11)))
# apply(exp(tn_data), 2, median)
# 
# long_data <- tn_data %>%
#   pivot_longer(cols = everything(), names_to = "group", values_to = "values")
# 
# max_value <- max(long_data$values, na.rm = TRUE)
# buffer <- (max_value - min(long_data$values, na.rm = TRUE)) * 0.1  # 10% buffer
# 
# pairwise_comparisons <- data.frame(
#   group1 = factor(c('tn_00', 'tn_01', 'tn_00')),
#   group2 = factor(c('tn_01', 'tn_11', 'tn_11')),
#   p.value = c(ttest_01_00$P.dbar, ttest_01_11$P.dbar, ttest_00_11$P.dbar),
#   signif_label = ifelse(c(ttest_01_00$P.dbar, ttest_01_11$P.dbar, ttest_00_11$P.dbar) < 0.05, "*", ""))
# 
# pairwise_comparisons <- pairwise_comparisons %>%
#   mutate(
#     y_position = max_value + seq(1, n(), by = 1) * buffer,  
#     midpoint = ((as.numeric(group1) + as.numeric(group2)) / 2) + 0.5
#   )
# 
# b <- ggplot(long_data, aes(x = group, y = values)) + 
#   geom_boxplot() +
#   geom_segment(data = pairwise_comparisons, 
#                aes(x = group1, xend = group2, y = y_position, yend = y_position), 
#                linetype = "solid") +
#   geom_text(data = pairwise_comparisons, 
#             aes(x = midpoint, y = y_position, label = signif_label)) +
#   ggtitle("b) Hidden state support in turn and rates")
# 
# grid.arrange(a, b, nrow = 1)


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

apply(trans_data, 2, median)

long_data <- trans_data %>%
  pivot_longer(cols = everything(), names_to = "group", values_to = "values")

ggplot(long_data, aes(x = group, y = log(values))) + 
  geom_boxplot() +
  ggtitle("b) Observed transition rates") +
  coord_cartesian(ylim=c(-5, 5))


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


