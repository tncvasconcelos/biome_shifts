
library(ape)
# Load the package (after installation, see above).
library(optimx)   # optimx seems better than R's default optim()
library(GenSA)    # GenSA seems better than optimx (but slower) on 5+ parameters, 
# seems to sometimes fail on simple problems (2-3 parameters)
library(FD)       # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)
library(rexpokit)
library(cladoRcpp)
# library(BioGeoBEARS)

resolve.names <- function(names_to_solve) {
  gnr_resolve_x <- function(x) {
    sources <- taxize::gnr_datasources()
    tmp.name <- suppressWarnings(taxize::gnr_resolve(names=x, data_source_ids=sources$id[sources$title == "GBIF Backbone Taxonomy"], best_match_only=TRUE)$matched_name)
    if(is.null(tmp.name)) {
      tmp.name <- paste0(x,"_UNMATCHED")
    }
    return(tmp.name)
  }
  all_names <- pbapply::pblapply(names_to_solve, gnr_resolve_x, cl=6)
  return(as.character(all_names))
}

load.trees <- function(tree.dir) {
  tree_files <- list.files(tree.dir, full.names = T)
  all_trees <- list()
  for(i in 1:length(tree_files)) {
    load(tree_files[i])
    if(exists("one_tree")) {
      all_trees[[i]] <- one_tree
      names(all_trees)[i] <- gsub(paste0(c(paste0(tree.dir,"/"), ".Rsave"), collapse="|"),"", tree_files[i])
      rm("one_tree")
    }
  }
  return(all_trees)
}


simplify.names.taxize <- function(names) {
  results <- c()
  for(name_index in 1:length(names)){
    one_tmp_string <- names[name_index]
    splitted_names <- strsplit(one_tmp_string," ")[[1]]
    genus <- splitted_names[1]
    epiphet <- splitted_names[2]
    if(any(grepl("indet_sp",splitted_names))) {
      full_name <- "tip_to_drop" # indet species
    } else {
      if(stringr::str_detect(epiphet,"[[:upper:]]")) {
        full_name <- "tip_to_drop" # indet species
      } else {
        if(length(splitted_names) == 2) { 
          full_name <- paste(c(genus, epiphet), collapse = " ")
        } else {
          if(length(splitted_names) > 2) {
            complement <- splitted_names[3:length(splitted_names)]
            if(grepl("[()]", complement[1])) {
              full_name <- paste(c(genus, epiphet), collapse = " ")
            } else {
              if(stringr::str_detect(complement[1],"[[:upper:]]")) {
                full_name <- paste(c(genus, epiphet), collapse = " ")
              } else {
                complement <- subset(complement, !stringr::str_detect(complement,"[[:upper:]]"))
                complement <- subset(complement, !grepl(paste(c("[()]","&","([0-9]+).*$","^ex$"), collapse="|"), complement))
                if(length(complement)==0){
                  full_name <- paste(c(genus, epiphet), collapse = " ")
                } else {
                  full_name <- paste(c(genus, epiphet, complement), collapse = " ")
                }
              }
            } 
          }
        }  
      }
    }
    results[name_index] <- full_name
  }
  return(results)
}


fix.names.taxize <- function(focal_species_trees) {
  for(name_index in 1:length(focal_species_trees)){
    one_tmp_string <- focal_species_trees[name_index]
    if(any(grepl("[()]", one_tmp_string))){
      splitted_names <- strsplit(one_tmp_string," ")[[1]]
      begin_author <- which(grepl("[()]", splitted_names))[1]
      species_name <- paste0(splitted_names[1:(begin_author-1)], collapse=" ")
      author <- splitted_names[begin_author:length(splitted_names)]
      old_authors <- author[grep("[()]", author)]
      end_first_half <- floor(length(old_authors)/2)
      before <- old_authors[1:end_first_half]
      after <- old_authors[(end_first_half+1):(length(old_authors))]
      if(paste(before,collapse = " ") == paste(after, collapse=" ")) {
        author <- paste(author[1:(length(author)/2)], collapse=" ")
        focal_species_trees[name_index] <- paste0(species_name, " ", author, collapse=" ")
      } else {
        author <- paste(author, collapse=" ")
        focal_species_trees[name_index] <- paste0(species_name, " ", author, collapse=" ")
      }
    }
  }
  return(focal_species_trees)
}

convert2Lambda <- function(pars){
  if(is.na(pars[1])){
    focal_pars <- sample(which(!is.na(pars)), size = 2, replace = FALSE)
    if(2 %in% focal_pars & 3 %in% focal_pars){
      # mu and div
      lambda <- pars[2] + pars[3]
    }
    if(2 %in% focal_pars & 4 %in% focal_pars){
      # mu and turn
      lambda <- pars[4] - pars[2]
    }
    if(2 %in% focal_pars & 5 %in% focal_pars){
      # mu and ef
      lambda <- pars[2]/pars[5]
    }
    if(3 %in% focal_pars & 4 %in% focal_pars){
      # div and turn
      lambda <- (pars[3]+pars[4])/2
    }
    if(3 %in% focal_pars & 5 %in% focal_pars){
      # div and ef
      lambda <- pars[3]/(1-pars[5])
    }
    if(4 %in% focal_pars & 5 %in% focal_pars){
      # turn and ef
      lambda <- pars[4]/(1+pars[5])
    }
  }else{
    lambda <- pars[1]
  }
  return(lambda)
}

convert2Mu <- function(pars){
  if(is.na(pars[2])){
    focal_pars <- sample(which(!is.na(pars)), size = 2, replace = FALSE)
    if(1 %in% focal_pars & 3 %in% focal_pars){
      # lambda and div
      mu <- pars[1] - pars[3]
    }
    if(1 %in% focal_pars & 4 %in% focal_pars){
      # lambda and turn
      mu <- pars[4] - pars[1]
    }
    if(1 %in% focal_pars & 5 %in% focal_pars){
      # lambda and ef
      mu <- pars[1]*pars[5]
    }
    if(3 %in% focal_pars & 4 %in% focal_pars){
      # div and turn
      mu <- (pars[4] - pars[3])/2
    }
    if(3 %in% focal_pars & 5 %in% focal_pars){
      # div and ef
      mu <- (pars[3]*pars[5])/(1-pars[5])
    }
    if(4 %in% focal_pars & 5 %in% focal_pars){
      # turn and ef
      mu <- (pars[4]*pars[5])/(1+pars[5])
    }
  }else{
    mu <- pars[2]
  }
  return(mu)
  
}

convertBetweenPars <- function(pars){
  # pars <- c("lambda", "mu", "net.div", "turn", "ef")
  if(length(which(!is.na(pars))) >= 3){
    warning("More than 2 paramaters are specified. Randomly choosing 2 for the calculations.")
  }
  if(is.na(pars[1])){
    lambda <- convert2Lambda(pars)
  }else{
    lambda <- pars[1]
  }
  if(is.na(pars[2])){
    mu <- convert2Mu(pars)
  }else{
    mu <- pars[2]
  }
  net.div <- lambda - mu
  turn <- lambda + mu
  ef <- mu/lambda
  out <- c(lambda=lambda, mu=mu, net.div=net.div, turn=turn, ef=ef)
  if(!setequal(round(out[which(!is.na(pars))], 5), round(pars[which(!is.na(pars))], 5))){
    stop("An error occured because the calculated output doesn't match the input. Please check that your input parameters can be combined in a way that is possible.")
  }
  return(out)
}

quickConvert <- function(par, par.class){
  base_p <- c(lambda=NA, mu=NA, net.div=NA, turn=NA, ef=NA)
  base_p[match(par.class, names(base_p))] <- par
  p <- convertBetweenPars(base_p)
  names(p) <- c("lambda", "mu", "net.div", "turn", "ef")
  return(p)
}

get_model_table <- function (model_list, type = "AIC") {
  nTip <- length(model_list[[1]]$phy$tip.label)
  AIC <- simplify2array(lapply(model_list, "[[", type))
  dAIC <- AIC - min(AIC)
  AICwt <- exp(-0.5 * dAIC)/sum(exp(-0.5 * dAIC))
  log10Lik <- simplify2array(lapply(model_list, "[[", "loglik"))
  out <- data.frame(lnLik = log10Lik, AIC = AIC, 
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


evaluate_difference <- function(model_list, type = "rate_class"){
  model_table <- get_model_table(model_list)
  par_list <- lapply(model_list, convert_muhisse_pars)
  out <- c()
  if(type == "rate_class"){
    for(i in 1:dim(model_table)[1]){
      if(length(par_list[[i]]$turnover) == 4){ #single rate class
        rate_class_a <- par_list[[i]]$trans_matrix[c(1,2,4), c(1,2,4)]
        rate_class_b <- par_list[[i]]$trans_matrix[c(1,2,4), c(1,2,4)] # same rate class
        turns_rate_a <- par_list[[i]]$turnover[c(1,2,4)]
        turns_rate_b <- par_list[[i]]$turnover[c(1,2,4)]
        ef_rate_a <- par_list[[i]]$eps[c(1,2,4)]
        ef_rate_b <- par_list[[i]]$eps[c(1,2,4)]
      }else{ # hidden states
        rate_class_a <- par_list[[i]]$trans_matrix[c(1,2,4), c(1,2,4)]
        rate_class_b <- par_list[[i]]$trans_matrix[c(5,6,8), c(5,6,8)] # same rate class
        turns_rate_a <- par_list[[i]]$turnover[c(1,2,4)]
        turns_rate_b <- par_list[[i]]$turnover[c(5,6,8)]
        ef_rate_a <- par_list[[i]]$eps[c(1,2,4)]
        ef_rate_b <- par_list[[i]]$eps[c(5,6,8)]
      }
      muhisse_pars <- data.frame(turn = c(turns_rate_a, turns_rate_b), ef = c(ef_rate_a, ef_rate_b))
      all_pars <- apply(muhisse_pars, 1, function(x) quickConvert(x, c("turn", "ef"))) |> t()
      all_pars_a <- all_pars[1:3,]
      all_pars_b <- all_pars[4:6,]
      mean_pars_a <- colMeans(all_pars_a, na.rm = TRUE)
      mean_pars_b <- colMeans(all_pars_b, na.rm = TRUE)
      names(mean_pars_a) <- paste0(names(mean_pars_a), "_a")
      names(mean_pars_b) <- paste0(names(mean_pars_b), "_b")
      out <- rbind(out, c(mean_pars_a,
                          trans_a = mean(rate_class_a, na.rm = TRUE),
                          mean_pars_b,
                          trans_b = mean(rate_class_b, na.rm = TRUE)))
    }
  }
  if(type == "obs"){
    for(i in 1:dim(model_table)[1]){
      if(length(par_list[[i]]$turnover) == 4){ #single rate class
        focal_mat <- par_list[[i]]$trans_matrix
        focal_mat[is.na(focal_mat)] <- 0
        diag(focal_mat) <- -rowSums(focal_mat)
        trans_00 <- -diag(focal_mat)[1]
        trans_01 <- -diag(focal_mat)[2]
        trans_11 <- -diag(focal_mat)[4]
        turns_00 <- par_list[[i]]$turnover[1]
        turns_01 <- par_list[[i]]$turnover[2]
        turns_11 <- par_list[[i]]$turnover[4]
        ef_00 <- par_list[[i]]$eps[1]
        ef_01 <- par_list[[i]]$eps[2]
        ef_11 <- par_list[[i]]$eps[4]
      }else{ # hidden states
        focal_mat <- par_list[[i]]$trans_matrix
        focal_mat[is.na(focal_mat)] <- 0
        diag(focal_mat) <- -rowSums(focal_mat)
        trans_00 <- -mean(diag(focal_mat)[c(1,5)])
        trans_01 <- -mean(diag(focal_mat)[c(2,6)])
        trans_11 <- -mean(diag(focal_mat)[c(4,8)])
        turns_00 <- mean(par_list[[i]]$turnover[c(1,5)])
        turns_01 <- mean(par_list[[i]]$turnover[c(2,6)])
        turns_11 <- mean(par_list[[i]]$turnover[c(4,8)])
        ef_00 <- mean(par_list[[i]]$eps[c(1,5)])
        ef_01 <- mean(par_list[[i]]$eps[c(2,6)])
        ef_11 <- mean(par_list[[i]]$eps[c(4,8)])
      }
      muhisse_pars <- data.frame(turn = c(turns_00, turns_01, turns_11), ef = c(ef_00, ef_01, ef_11))
      all_pars <- apply(muhisse_pars, 1, function(x) quickConvert(x, c("turn", "ef"))) |> t()
      all_pars_00 <- all_pars[1,]
      all_pars_01 <- all_pars[2,]
      all_pars_11 <- all_pars[3,]
      names(all_pars_00) <- paste0(names(all_pars_00), "_00")
      names(all_pars_01) <- paste0(names(all_pars_01), "_01")
      names(all_pars_11) <- paste0(names(all_pars_11), "_11")
      out <- rbind(out, c(all_pars_00, 
                          trans_00 = as.numeric(trans_00),
                          all_pars_01,
                          trans_01 = as.numeric(trans_01),
                          all_pars_11,
                          trans_11 = as.numeric(trans_11)))
    }
  }
  if(type == "trans"){
    for(i in 1:dim(model_table)[1]){
      if(length(par_list[[i]]$turnover) == 4){ #single rate class
        focal_mat <- par_list[[i]]$trans_matrix
        trans_rates <- c(f00_t01_a = focal_mat[1,2], 
                         f01_t00_a = focal_mat[2,1],
                         f01_t11_a = focal_mat[2,4],
                         f11_t01_a = focal_mat[4,2],
                         f00_t01_b = focal_mat[1,2], 
                         f01_t00_b = focal_mat[2,1],
                         f01_t11_b = focal_mat[2,4],
                         f11_t01_b = focal_mat[4,2])
      }else{ # hidden states
        focal_mat <- par_list[[i]]$trans_matrix
        trans_rates <- c(f00_t01_a = focal_mat[1,2], 
                         f01_t00_a = focal_mat[2,1],
                         f01_t11_a = focal_mat[2,4],
                         f11_t01_a = focal_mat[4,2],
                         f00_t01_b = focal_mat[5,6], 
                         f01_t00_b = focal_mat[6,5],
                         f01_t11_b = focal_mat[6,8],
                         f11_t01_b = focal_mat[8,6])
      }
      out <- rbind(out, trans_rates)
    }
  }
  if(type == "all"){
    for(i in 1:dim(model_table)[1]){
      if(length(par_list[[i]]$turnover) == 4){ #single rate class
        focal_mat <- par_list[[i]]$trans_matrix
        trans_rates <- c(f00_t01_a = focal_mat[1,2], 
                         f01_t00_a = focal_mat[2,1],
                         f01_t11_a = focal_mat[2,4],
                         f11_t01_a = focal_mat[4,2],
                         f00_t01_b = focal_mat[1,2], 
                         f01_t00_b = focal_mat[2,1],
                         f01_t11_b = focal_mat[2,4],
                         f11_t01_b = focal_mat[4,2])
        turns <- c(tn_00_a = par_list[[i]]$turnover[1],
                   tn_01_a = par_list[[i]]$turnover[2],
                   tn_11_a = par_list[[i]]$turnover[4],
                   tn_00_b = par_list[[i]]$turnover[1],
                   tn_01_b = par_list[[i]]$turnover[2],
                   tn_11_b = par_list[[i]]$turnover[4])
        ef <- c(ef_00_a = par_list[[i]]$eps[1], 
                ef_01_a = par_list[[i]]$eps[2],
                ef_11_a = par_list[[i]]$eps[4],
                ef_00_b = par_list[[i]]$eps[1], 
                ef_01_b = par_list[[i]]$eps[2],
                ef_11_b = par_list[[i]]$eps[4])
      }else{ # hidden states
        focal_mat <- par_list[[i]]$trans_matrix
        trans_rates <- c(f00_t01_a = focal_mat[1,2], 
                         f01_t00_a = focal_mat[2,1],
                         f01_t11_a = focal_mat[2,4],
                         f11_t01_a = focal_mat[4,2],
                         f00_t01_b = focal_mat[5,6], 
                         f01_t00_b = focal_mat[6,5],
                         f01_t11_b = focal_mat[6,8],
                         f11_t01_b = focal_mat[8,6])
        turns <- c(tn_00_a = par_list[[i]]$turnover[1],
                   tn_01_a = par_list[[i]]$turnover[2],
                   tn_11_a = par_list[[i]]$turnover[4],
                   tn_00_b = par_list[[i]]$turnover[5],
                   tn_01_b = par_list[[i]]$turnover[6],
                   tn_11_b = par_list[[i]]$turnover[8])
        ef <- c(ef_00_a = par_list[[i]]$eps[1], 
                ef_01_a = par_list[[i]]$eps[2],
                ef_11_a = par_list[[i]]$eps[4],
                ef_00_b = par_list[[i]]$eps[5], 
                ef_01_b = par_list[[i]]$eps[6],
                ef_11_b = par_list[[i]]$eps[8])
      }
      out <- rbind(out, c(trans_rates,turns,ef))
    }
  }
  return(colSums(out * model_table$AICwt))
}


quick_check <- function(model_path){
  load(model_path)
  return(all(unlist(lapply(res, function(x) class(x) == "try-error"))))
}



"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(ymin = min(y),
                     ymax = max(y),
                     xmin = x - width / 2,
                     xmax = x)
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data, 
                              xmaxv = x,
                              xminv = x + violinwidth * (xmin - x))
            
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                             plyr::arrange(transform(data, x = xmaxv), -y))
            
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1,])
            
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                            alpha = NA, linetype = "solid"),
          
          required_aes = c("x", "y")
  )
