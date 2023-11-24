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

