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

quick_check <- function(model_path){
  load(model_path)
  return(all(unlist(lapply(res, function(x) class(x) == "try-error"))))
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
