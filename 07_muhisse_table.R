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
plot_data$ntip <- n_tips

write.csv(plot_data, file = "all_par_table.csv", row.names = FALSE)
colMeans(plot_data[,-1])



