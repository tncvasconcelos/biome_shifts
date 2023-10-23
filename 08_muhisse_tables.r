# rm(list=ls())
setwd("~/biome_shifts/")

library(ape)
library(hisse)
library(parallel)
library(phytools)
library(phylolm)

get_model_table <- function (model_list, type = "AIC") {
  nTip <- length(model_list[[1]]$phy$tip.label)
  AIC <- simplify2array(lapply(model_list, "[[", type))
  dAIC <- AIC - min(AIC)
  AICwt <- exp(-0.5 * dAIC)/sum(exp(-0.5 * dAIC))
  LogLik <- simplify2array(lapply(model_list, "[[", "loglik"))
  out <- data.frame(lnLik = LogLik, AIC = AIC, 
                    dAIC = dAIC, AICwt = AICwt)
  colnames(out) <- gsub("AIC", type, colnames(out))
  return(out)
}


convert_muhisse_pars <- function(model_res){
  est_pars <- model_res$solution[model_res$index.par < max(model_res$index.par)]
  trans_matrix <- model_res$trans.matrix
  hidden_states <- dim(trans_matrix)[1]/4
  trans_rates <- est_pars[grep("q", names(est_pars))]
  turnover <- eps <- setNames(rep(0, 4 * hidden_states), colnames(trans_matrix))
  turnover_tmp <- est_pars[grep("turnover", names(est_pars))]
  eps_tmp <- est_pars[grep("eps", names(est_pars))]
  if(hidden_states == 1){
    names(trans_rates) <- gsub("A", "", names(trans_rates))
    names(turnover_tmp) <- gsub("A", "", names(turnover_tmp))
    names(eps_tmp) <- gsub("A", "", names(eps_tmp))
    names(trans_rates) <- gsub("q", "", names(trans_rates))
    names(turnover_tmp) <- gsub("turnover", "", names(turnover_tmp))
    names(eps_tmp) <- gsub("eps", "", names(eps_tmp))
    from_to <- do.call(rbind, strsplit(names(trans_rates), "_"))
    for(i in 1:length(eps_tmp)){
      eps[grep(names(eps_tmp[i]), names(eps))] <- eps_tmp[i]
      turnover[grep(names(turnover_tmp[i]), names(eps))] <- turnover_tmp[i]
    }
  }else{
    names(trans_rates) <- gsub("q", "", names(trans_rates))
    names(turnover_tmp) <- gsub("turnover", "", names(turnover_tmp))
    names(eps_tmp) <- gsub("eps", "", names(eps_tmp))
    from_to <- do.call(rbind, strsplit(names(trans_rates), "_"))
    for(i in 1:length(eps_tmp)){
      eps[grep(names(eps_tmp[i]), names(eps))] <- eps_tmp[i]
      turnover[grep(names(turnover_tmp[i]), names(eps))] <- turnover_tmp[i]
    }
  }
  for(i in 1:dim(from_to)[1]){
    focal_trans <- from_to[i,]
    from_index <- grep(focal_trans[1], rownames(trans_matrix))
    to_index <- grep(focal_trans[2], rownames(trans_matrix))
    trans_matrix[from_index, to_index] <- trans_rates[i]
  }
  trans_matrix[trans_matrix == 0] <- NA
  return(list(turnover = turnover, eps = eps, trans_matrix = trans_matrix))
}

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


test_hypotheses <- function(model_res){
  converted_model <- convert_muhisse_pars(model_res)
  n_hidden <- dim(model_res$trans.matrix)[1]/4
  both_differ <- eps_differ <- turn_differ <-trans_differ <- hidden_incl <- FALSE
  if(n_hidden > 1){
    hidden_incl <- TRUE
    trans_ul <- model_res$trans.matrix[1:4, 1:4]
    trans_lr <- model_res$trans.matrix[5:8, 5:8]
    trans_differ <- !all(trans_ul == trans_lr, na.rm = TRUE)
    turn_differ <- !all(converted_model$turnover[1:4] == converted_model$turnover[5:8])
    eps_differ <- !all(converted_model$eps[1:4] == converted_model$eps[5:8])
  }
  both_differ <- trans_differ & turn_differ
  return(c(hidden_incl = hidden_incl,  eps_differ=eps_differ, turn_differ=turn_differ, trans_differ=trans_differ, both_differ=both_differ))
}

# Load the necessary packages
library(ggplot2)
library(dplyr)
library(tidyr)

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
result_matrix$RowNames <- rownames(result_matrix)
# Reshape the data from wide to long format
result_matrix_long <- result_matrix %>%
  pivot_longer(cols = -RowNames, names_to = "Variable", values_to = "Value")

result_matrix_long$Variable <- factor(result_matrix_long$Variable, unique(result_matrix_long$Variable))
# Create the heatmap
#D01B1B, #FF4242, #FFFFFF, #e7f9ff, #95D2EC and #47abd8. 
p <- ggplot(result_matrix_long, aes(x = Variable, y = RowNames, fill = Value)) +
  geom_tile() +
  scale_fill_viridis_c(option = "E", name = expression(sum(AIC[wt]))) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
print(p)
ggsave("plots/08_aic_support_plot.pdf", p, height = 8, width = 4)

# look at the ones where our hypothesis is potentially supportable and see whether rates are higher in one than the other. need to consider how to measure "rates"
both_differ_model_list <- all_model_list[result_matrix$both_differ > 0.5]
both_differ_clade_names <- clade_names[result_matrix$both_differ > 0.5]
both_differ_model_table <- lapply(both_differ_model_list, get_model_table)

evaluate_difference <- function(model_list, return="mean"){
  model_table <- get_model_table(model_list)
  par_list <- lapply(model_list, convert_muhisse_pars)
  out <- c()
  for(i in 1:dim(model_table)[1]){
    if(length(par_list[[i]]$turnover) == 4){ #single rate class
      rate_class_a <- par_list[[i]]$trans_matrix[c(1,2,4), c(1,2,4)]
      rate_class_b <- par_list[[i]]$trans_matrix[c(1,2,4), c(1,2,4)] # same rate class
      turns_rate_a <- par_list[[i]]$turnover[c(1,2,4)]
      turns_rate_b <- par_list[[i]]$turnover[c(1,2,4)]
    }else{ # hidden states
      rate_class_a <- par_list[[i]]$trans_matrix[c(1,2,4), c(1,2,4)]
      rate_class_b <- par_list[[i]]$trans_matrix[c(5,6,8), c(5,6,8)] # same rate class
      turns_rate_a <- par_list[[i]]$turnover[c(1,2,4)]
      turns_rate_b <- par_list[[i]]$turnover[c(5,6,8)]
    }
    out <- rbind(out, c(trans_a = mean(rate_class_a, na.rm = TRUE), 
                    trans_b = mean(rate_class_b, na.rm = TRUE),
                    turns_a = mean(turns_rate_a),
                    turns_b = mean(turns_rate_b)))
  }
  return(colSums(out * model_table$AICwt))
}

plot_data <- cbind(Group = clade_names, as.data.frame(do.call(rbind, lapply(all_model_list, evaluate_difference))))
plot_data <- as.data.frame(plot_data)


ggplot(plot_data, aes(x = log(trans_a), y = log(turns_a))) +
  geom_point(aes(color = Group, fill = Group), size = 3) +
  geom_point(data = plot_data, aes(x = log(trans_b), y = log(turns_b), color = Group, fill = Group), size = 3) +
  labs(x = "log(transition rates)", y = "log(turnover rates)", color = "Group") +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_segment(aes(xend = log(trans_b), yend = log(turns_b)), linetype = "dotted", color = "black") +
  geom_text(aes(x = log(trans_b), y = log(turns_b), label = Group), hjust = -0.2, vjust = 0.5, size = 3)


phy_bb <- read.tree("backbone_tree.tre")
dat_names <- unlist(lapply(strsplit(paste(plot_data$Group), " "), function(x) x[[1]]))
phy_bb$tip.label <- unlist(lapply(strsplit(phy_bb$tip.label, "-"), function(x) x[[1]]))
phy_bb$tip.label <- unlist(lapply(strsplit(phy_bb$tip.label, "_"), function(x) x[[1]]))
to_drop <- phy_bb$tip.label[!phy_bb$tip.label %in% dat_names]
phy_bb <- drop.tip(phy_bb, to_drop)
plot(phy_bb)
tip_names <- phy_bb$tip.label
for(focal_tip in tip_names){
  index <- which(phy_bb$tip.label == focal_tip)
  phy_bb$tip.label[index] <- paste0(focal_tip, "_A")
  phy_bb <- bind.tip(phy_bb, tip.label = paste0(focal_tip, "_B"), edge.length = 1e-2, where = index)
}
plot(phy_bb, type = "fan")

lm_dat <- data.frame(row.names = c(paste0(dat_names, "_A"), paste0(dat_names, "_B")),
                     trans_rt = c(plot_data$trans_a, plot_data$trans_b),
                     turn_rt = c(plot_data$turns_a, plot_data$turns_b))
lm_dat <- lm_dat[match(phy_bb$tip.label, rownames(lm_dat)),]

fit = phylolm(trans_rt~turn_rt, data=lm_dat, phy=phy_bb)
summary(fit)

