
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

# WWFload is taken from speciesgeocodeR; all credit goes to the original authors
WWFload <- function(x = NULL) {
  if (missing(x)) {
    x <- getwd()
  }
  download.file("http://assets.worldwildlife.org/publications/15/files/original/official_teow.zip",
                destfile = file.path(x, "wwf_ecoregions.zip"), quiet=TRUE)
  unzip(file.path(x, "wwf_ecoregions.zip"), exdir = file.path(x, "WWF_ecoregions"))
  file.remove(file.path(x, "wwf_ecoregions.zip"))
  wwf <- maptools::readShapeSpatial(file.path(x, "WWF_ecoregions", "official",
                                              "wwf_terr_ecos.shp"))
  return(wwf)
}


localityToBiome <- function (points, lat="lat",lon="lon") {
  #colnames(points) <- c("acceptedScientificName","key","decimalLatitude","decimalLongitude","basisOfRecord","issues")
  cat("Getting biome from locality data...")
  points[,lat] <-  as.numeric(points[,lat])
  points[,lon] <-  as.numeric(points[,lon])
  locations.spatial <- sp::SpatialPointsDataFrame(coords=points[,c(which(colnames(points)==lon), which(colnames(points)==lat))], data=points)
  wwf <- WWFload(tempdir())
  mappedregions <- sp::over(locations.spatial, wwf)
  biomes <- c("Tropical & Subtropical Moist Broadleaf Forests", "Tropical & Subtropical Dry Broadleaf Forests", "Tropical & Subtropical Coniferous Forests", "Temperate Broadleaf & Mixed Forests", "Temperate Conifer Forests", "Boreal Forests/Taiga", "Tropical & Subtropical Grasslands, Savannas & Shrubland", "Temperate Grasslands, Savannas & Shrublands", "Flooded Grasslands & Savannas", "Montane Grasslands & Shrublands", "Tundra", "Mediterranean Forests, Woodlands & Scrub", "Deserts & Xeric Shrublands", "Mangroves")
  points$eco_name <- mappedregions$ECO_NAME
  points$biome <- biomes[mappedregions$BIOME]
  return(points)
}


# getting biomes for each species
getBiomes <- function (points, species="species") {
  cat("Summarizing biome from locality data...")
  points <- as.data.frame(points) # not sure how to do it without transforming back to data.frame
  points <- subset(points, !is.na(points[,"biome"]))
  categories <- unique(points[,"biome"])
  taxa <- as.character(unique(points[,species]))
  result <- matrix(0, nrow=length(taxa), ncol=length(categories))
  rownames(result) <- taxa
  colnames(result) <- categories
  cat("\n")
  for (taxon_index in seq_along(taxa)) {
    for (category_index in seq_along(categories)) {
      x0 <- points[,species]==taxa[taxon_index]
      x1 <- points[,"biome"]==categories[category_index]
      result[taxon_index, category_index] <- length(which(x0 & x1))
    }
    cat(taxon_index, "\r")
  }
  return(result)
}

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


evaluate_difference <- function(model_list, type = "rate_class", mod_avg = TRUE){
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
  if(mod_avg){
    return(colSums(out * model_table$AICwt))
  }else{
    tmp <- vector("numeric", length(model_table$AICwt))
    tmp[which.max(model_table$AICwt)] <- 1
    return(colSums(out * tmp))
  }
}


quick_check <- function(model_path){
  load(model_path)
  return(all(unlist(lapply(res, function(x) class(x) == "try-error"))))
}


# 
# "%||%" <- function(a, b) {
#   if (!is.null(a)) a else b
# }
# 
# geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
#                              position = "dodge", trim = TRUE, scale = "area",
#                              show.legend = NA, inherit.aes = TRUE, ...) {
#   layer(
#     data = data,
#     mapping = mapping,
#     stat = stat,
#     geom = GeomFlatViolin,
#     position = position,
#     show.legend = show.legend,
#     inherit.aes = inherit.aes,
#     params = list(
#       trim = trim,
#       scale = scale,
#       ...
#     )
#   )
# }
# 
# GeomFlatViolin <-
#   ggproto("GeomFlatViolin", Geom,
#           setup_data = function(data, params) {
#             data$width <- data$width %||%
#               params$width %||% (resolution(data$x, FALSE) * 0.9)
#             
#             # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
#             data %>%
#               group_by(group) %>%
#               mutate(ymin = min(y),
#                      ymax = max(y),
#                      xmin = x - width / 2,
#                      xmax = x)
#           },
#           
#           draw_group = function(data, panel_scales, coord) {
#             # Find the points for the line to go all the way around
#             data <- transform(data, 
#                               xmaxv = x,
#                               xminv = x + violinwidth * (xmin - x))
#             
#             # Make sure it's sorted properly to draw the outline
#             newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
#                              plyr::arrange(transform(data, x = xmaxv), -y))
#             
#             # Close the polygon: set first and last point the same
#             # Needed for coord_polar and such
#             newdata <- rbind(newdata, newdata[1,])
#             
#             ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
#           },
#           
#           draw_key = draw_key_polygon,
#           
#           default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
#                             alpha = NA, linetype = "solid"),
#           
#           required_aes = c("x", "y")
#   )


FilterWCVP_genus <- function(points, all_vars, twgd_data, lon="decimalLongitude", lat="decimalLatitude") {
  npoints_start <- nrow(points)
  tmp_points = as.data.frame(points)
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  tmp_points = subset(tmp_points, !is.na(tmp_points$x))
  tmp_points = subset(tmp_points, !is.na(tmp_points$y))
  # Load shape files and make sure they have the same name as the WCVP column with the TDWG areas
  #twgd_data <- suppressWarnings(maptools::readShapeSpatial(path))
  dubiousGBIF_ids <- c()
  all_genera <- unique(all_vars$genus)
  all_vars_genus_level <- all_vars[all_vars$taxon_rank=="Genus",]
  for(genus_index in 1:length(all_genera)) {
    gbif_subset <- subset(tmp_points, tmp_points$genus == all_genera[genus_index])
    if(nrow(gbif_subset)!=0) {
      wcvp_subset <- subset(all_vars_genus_level, all_vars_genus_level$genus == all_genera[genus_index])
      wcvp_subset <- subset(wcvp_subset, wcvp_subset$introduced==0)
      wcvp_subset <- subset(wcvp_subset, wcvp_subset$extinct==0)
      wcvp_subset <- subset(wcvp_subset, wcvp_subset$location_doubtful==0)
      occ_areas <- wcvp_subset$area_code_l3
      area_plus_buffer <- twgd_data[which(as.character(twgd_data$LEVEL3_COD) %in% occ_areas),]
      if(nrow(area_plus_buffer)>0) {
        coords <- gbif_subset[,c("x","y")]
        sp::coordinates(coords) <- ~ x + y
        answer <- which(is.na(sp::over(coords, area_plus_buffer)[,3]))
        if(length(answer) != 0) {
          dubiousGBIF_ids <- c(dubiousGBIF_ids, as.character(gbif_subset$gbifID[answer]))
        }
      }
    }
    cat(genus_index, "\r")
  }
  cleaned_points <- subset(points, !as.character(points$gbifID) %in% dubiousGBIF_ids)
  npoints_end <- nrow(cleaned_points)
  print(paste0(npoints_start - npoints_end, " points removed."))
  return(cleaned_points)
}


FilterWCVP <- function(points, all_vars, reference_table, twgd_data, species= "scientificName", lon="decimalLongitude", lat="decimalLatitude") {
  npoints_start <- nrow(points)
  tmp_points = as.data.frame(points)
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  tmp_points = subset(tmp_points, !is.na(tmp_points$x))
  tmp_points = subset(tmp_points, !is.na(tmp_points$y))
  # Load shape files and make sure they have the same name as the WCVP column with the TDWG areas
  #twgd_data <- suppressWarnings(maptools::readShapeSpatial(path))
  dubiousGBIF_ids <- c()
  for(species_index in 1:nrow(reference_table)) {
    gbif_subset <- subset(tmp_points, tmp_points$scientificName == reference_table$gbif_name[species_index])
    if(nrow(gbif_subset)!=0) {
      wcvp_subset <- subset(all_vars, all_vars$taxon_name == reference_table$wcvp_name[species_index])
      wcvp_subset <- subset(wcvp_subset, wcvp_subset$introduced==0)
      wcvp_subset <- subset(wcvp_subset, wcvp_subset$extinct==0)
      wcvp_subset <- subset(wcvp_subset, wcvp_subset$location_doubtful==0)
      occ_areas <- wcvp_subset$area_code_l3
      area_plus_buffer <- twgd_data[which(as.character(twgd_data$LEVEL3_COD) %in% occ_areas),]
      if(nrow(area_plus_buffer)>0) {
        coords <- gbif_subset[,c("x","y")]
        sp::coordinates(coords) <- ~ x + y
        answer <- which(is.na(sp::over(coords, area_plus_buffer)[,3]))
        if(length(answer) != 0) {
          dubiousGBIF_ids <- c(dubiousGBIF_ids, as.character(gbif_subset$gbifID[answer]))
        }
      }
    }
    cat(species_index, "\r")
  }
  cleaned_points <- subset(points, !as.character(points$gbifID) %in% dubiousGBIF_ids)
  npoints_end <- nrow(cleaned_points)
  print(paste0(npoints_start - npoints_end, " points removed."))
  return(cleaned_points)
}

# Getting climate means per TWGD region
MeanRasterWCVP <- function(path_raster="3_Landscape_instability/bio_1_instability.tif", path_tdwg="wgsrpd-master/level3/level3.shp") {
  twgd_data <- suppressWarnings(maptools::readShapeSpatial(path_tdwg))
  raster_example <- raster(path_raster)
  template_map <- raster_example
  template_map[!is.na(template_map[])] <- 0
  template_map <- aggregate(template_map, fact=25)
  raster_list <- list()
  for(area_index in 1:length(twgd_data)){
    one_area <- twgd_data[area_index,]
    cropped_raster <- NULL
    try(cropped_raster <- mask(crop(raster_example, one_area), one_area))
    if(!is.null(cropped_raster)) {
      mean_value <- mean(subset(cropped_raster[], !is.na(cropped_raster[])))
      template_raster <- cropped_raster
      template_raster[!is.na(template_raster[])] <- mean_value
      template_raster <- aggregate(template_raster, fact=25)
      template_raster <- raster::resample(template_raster, template_map)
      raster_list[[area_index]] <- raster::mask(template_raster, template_map) 
      print(area_index)
    } else { next }
  }
  raster_list[which(unlist(lapply(raster_list, is.null)))] <- NULL
  mm <- do.call(merge, raster_list)
  return(mm)
}

############################
#' Removes points in the sea
#' @param points A data.frame of distribution points with at least three columns where one column represents species names and the other two decimal coordinates
#' @param lat The name of the column in the data.frame with latitudes
#' @param lon The name of the column in the data.frame with longitudes
#' @param buffer A number of degrees around continental areas where points are still kept after filtering
RemoveSeaPoints <- function(points, lon="decimalLongitude", lat="decimalLatitude", buffer=0) {
  npoints_start <- nrow(points)
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  wrld_map <- rworldmap::getMap(resolution="low") # leaving both maps in the code for now, should probably drop one of them later
  coords <- tmp_points[,c("x","y")]
  sp::coordinates(coords) <- ~ x + y
  sp::proj4string(coords) <- sp::proj4string(wrld_map)
  country_plus_buffer <- raster::buffer(wrld_map, buffer) # adding buffer around polygons
  answer <- which(is.na(sp::over(coords, country_plus_buffer)))
  if(length(answer) > 0) {
    points <- points[-answer,]
    npoints_end <- nrow(points)
    print(paste0(npoints_start - npoints_end, " points removed."))
    return(points)
  } else {
    print("no points removed")
    return(points) }
}

#' Removes points that have 0 for both latitude and longitude
#' @param points A data.frame of distribution points with at least three columns where one column represents species names and the other two decimal coordinates
#' @param lat The name of the column in the data.frame with latitudes
#' @param lon The name of the column in the data.frame with longitudes
RemoveZeros <- function(points, lon="decimalLongitude", lat="decimalLatitude") {
  if(!inherits(points, "data.frame")) {
    stop("Argument points is not a data.frame.")
  }
  npoints_start <- nrow(points)
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  if(any(tmp_points$x==0 & tmp_points$y==0)) {
    points <- points[-which(tmp_points$x==0 & tmp_points$y==0),]
  }
  npoints_end <- nrow(points)
  print(paste0(npoints_start - npoints_end, " points removed."))
  return(points)
}

#' Removes outliers
#' @param points A data.frame of distribution points with at least three columns where one column represents species names and the other two decimal coordinates
#' @param species The name of the column in the data.frame with the names of species
#' @param lat The name of the column in the data.frame with latitudes
#' @param lon The name of the column in the data.frame with longitudes
RemoveOutliers <- function(points, species="scientificName", lon="decimalLongitude", lat="decimalLatitude") {
  npoints_start <- nrow(points)
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  colnames(tmp_points)[colnames(tmp_points)==species] <- "tmp_species"
  spp <- unique(tmp_points$species)
  all_points <- list()
  for(species_index in 1:length(spp)){
    sp0 <- tmp_points[tmp_points$species==spp[species_index],]
    out_lat <- grDevices::boxplot.stats(sp0$y)$out
    out_lon <- grDevices::boxplot.stats(sp0$x)$out
    sp0 <- sp0[!sp0$y %in% out_lat, ]
    sp0 <- sp0[!sp0$x %in% out_lon, ]
    all_points[[species_index]] <- sp0
  }
  points <- do.call(rbind, all_points)
  colnames(points)[colnames(points)=="x"] <- lon
  colnames(points)[colnames(points)=="y"] <- lat
  colnames(points)[colnames(points)=="tmp_species"] <- species
  npoints_end <- nrow(points)
  print(paste0(npoints_start - npoints_end, " points removed."))
  return(points)
}

#' Removes points that are located in the wrong country according to their GBIF labels
#' That will remove points that are not located in the countries where their labels say they were collected
#' @param points A data.frame of distribution points with at least five columns where one column represents species names and other two decimal coordinates.
#' @param lat The name of the column in the data.frame with latitudes
#' @param lon The name of the column in the data.frame with longitudes
#' @param buffer A number of degrees around each country where points are still considered part of that country
#' @details The input data.frame must have a column named countryCode and one named gbifID, as the .csv files downloaded directly from GBIF.
RemoveWrongCountries <- function(points, lon="decimalLongitude", lat="decimalLatitude", buffer=5, wrld_simpl="") {
  npoints_start <- nrow(points)
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  countries <- as.character(wrld_simpl[2]$ISO2)
  dubiousGBIF_ids <- c()
  for(country_index in 1:length(countries)) {
    tmp_country <- countries[country_index]
    if(length(which(tmp_points$countryCode %in% tmp_country)) > 0) {
      tmp_subset <- tmp_points[tmp_points$countryCode==tmp_country,]
      coords <- stats::na.omit(tmp_subset[,c("x","y")])
      sp::coordinates(coords) <- ~ x + y
      sp::proj4string(coords) <- sp::proj4string(wrld_simpl) <- "+proj=longlat +ellps=WGS84 +no_defs" 
      country_plus_buffer <- raster::buffer(wrld_simpl[country_index,], buffer) # adding buffer around country
      answer <- which(is.na(sp::over(coords, country_plus_buffer)))
      dubiousGBIF_ids <- c(dubiousGBIF_ids, tmp_subset$gbifID[answer])
    }
  }
  if(!is.null(dubiousGBIF_ids)) {
    points_cleaned <- points[-which(points$gbifID %in% dubiousGBIF_ids),]
  }
  npoints_end <- nrow(points_cleaned)
  print(paste0(npoints_start - npoints_end, " points removed."))
  return(points_cleaned)
}

#' Removes points that are located in country centroids
#' @param points A data.frame of distribution points with at least three columns where one column represents species names and the other two decimal coordinates
#' @param lat The name of the column in the data.frame with latitudes
#' @param lon The name of the column in the data.frame with longitudes
#' @param buffer A number in meters around each country centroid for points to be removed
RemoveCentroids <- function(points, lon="decimalLongitude", lat="decimalLatitude", buffer=75000) {
  npoints_start <- nrow(points)
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  wrld_map <- rworldmap::getMap(resolution="low")
  # here the buffer is in meters
  # Note: probably should do something about small countries where the buffer of 75km may be too broad
  coords <- tmp_points[,c("x","y")]
  sp::coordinates(coords) <- ~ x + y
  sp::proj4string(coords) <- sp::proj4string(wrld_map)
  centroids <- rgeos::gCentroid(wrld_map, byid=TRUE)
  centroids_plus_buffer <- raster::buffer(centroids, buffer) # adding buffer around polygons
  answer <- which(!is.na(sp::over(coords, centroids_plus_buffer)))
  if(length(answer) > 0) {
    points <- points[-answer,]
    npoints_end <- nrow(points)
    print(paste0(npoints_start - npoints_end, " points removed."))
    return(points)
  } else {
    print("no points removed")
    return(points) }
}

#' Removes duplicated latitudes and longitudes for the same species
#' @param points A data.frame of distribution points 
#' @param species The name of the column in the data.frame with the names of species
#' @param lat The name of the column in the data.frame with latitudes
#' @param lon The name of the column in the data.frame with longitudes
RemoveDuplicates <- function(points, species="scientificName", lon="decimalLongitude", lat="decimalLatitude") {
  npoints_start <- nrow(points)
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  colnames(tmp_points)[colnames(tmp_points)==species] <- "tmp_species"
  spp <- unique(tmp_points$species)
  all_points <- list()
  for(species_index in 1:length(spp)){
    tmp_subset <- as.data.frame(tmp_points[tmp_points$species==spp[species_index],])
    all_points[[species_index]] <- tmp_subset[-which(duplicated(tmp_subset[,c("x","y")])),]
  }
  points <- do.call(rbind, all_points)
  colnames(points)[colnames(points)=="x"] <- lon
  colnames(points)[colnames(points)=="y"] <- lat
  colnames(points)[colnames(points)=="tmp_species"] <- species
  npoints_end <- nrow(points)
  print(paste0(npoints_start - npoints_end, " points removed."))
  return(points)
}

#' Removes points with coordinates without decimal cases (probably innacurate)
#' @param points A data.frame of distribution points 
#' @param lat The name of the column in the data.frame with latitudes
#' @param lon The name of the column in the data.frame with longitudes
RemoveNoDecimal <- function(points, lon="decimalLongitude", lat="decimalLatitude") {
  npoints_start <- nrow(points)
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  length_decimal_lat <- nchar(sub("^[^.]*", "", tmp_points$y))
  length_decimal_lon <- nchar(sub("^[^.]*", "", tmp_points$x))
  points <- points[which(length_decimal_lat>=1 & length_decimal_lon>=1),]
  npoints_end <- nrow(points)
  print(paste0(npoints_start - npoints_end, " points removed."))
  return(points)
}


#' Thinning distribution data to smooth sampling bias
#' @param points A data.frame of three columns containing species and coordinates
#' @param species A character string indicating name of column with species names
#' @param lat A character string indicating name of column with latitudes
#' @param lon character string indicating name of column with longitudes
#' @param n A number indicating how many points to keep in each cell after thinning
Thinning <- function(points, species="species", lat = "decimalLatitude", lon="decimalLongitude", n = 3) {
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  spp <- unique(tmp_points[,species])
  results <- list()
  for(species_index in 1:length(spp)) {
    coords <- tmp_points[tmp_points[,species]==spp[species_index],c("y","x")]
    coords <- coords[!duplicated(coords[,"x"]) & !duplicated(coords[,"y"]),]
    if(nrow(coords) > 1) {
      sp::coordinates(coords) <- ~ y + x
      raster::crs(coords) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
      r0 <- raster::raster(coords)
      raster::res(r0) <- 1 # cell resolution
      r0 <- raster::extend(r0, raster::extent(r0) + 5) 
      res <- cbind(spp[species_index], as.data.frame(dismo::gridSample(coords, r0, n))) # n = maximum number of points per cell
      colnames(res) <- c("species", "lat","lon")
      results[[species_index]] <- res
    } else {
      res <- cbind(spp[species_index],coords)
      colnames(res) <- c("species", "lat","lon")
      results[[species_index]] <- res
    }
  }
  results <- do.call(rbind, results)
  return(results)
}


placeImageAtPoint <- function(img, points, pointName, imgWidth = 1, imgHeight = 1) {
  # Find the point by name
  pointRow <- which(rownames(points) == pointName)
  if(length(pointRow) == 0) {
    stop("Point name not found")
  }
  
  # Original point coordinates
  origX <- points[pointRow, 1]
  origY <- points[pointRow, 2]
  
  # Print the coordinates
  cat("Coordinates of", pointName, ": (", origX, ",", origY, ")\n")
  
  # Prompt for a click
  cat("Please click where you want to place the image for", pointName, "\n")
  loc <- locator(1)
  
  # Connect the original point to the image center with a dashed black line
  segments(origX, origY, loc$x, loc$y, col = "black", lty = 2)
  
  # Calculate corners to center the image on the click
  left <- loc$x - imgWidth / 2
  bottom <- loc$y - imgHeight / 2
  right <- loc$x + imgWidth / 2
  top <- loc$y + imgHeight / 2
  
  # Place the image
  rasterImage(img, left, bottom, right, top)
  
  # Place the point's name above the image
  text(loc$x, top + 0.05 * (top - bottom), labels = pointName, pos = 3)
}

getCoordsAtClick <- function(){
  coords <- locator(1)
  return(coords)
}

node_labels_custom <- function (text, node, adj = c(0.5, 0.5), frame = "rect", pch = NULL, thermo = NULL, pie = NULL, piecol = NULL, col = "black", bg = "lightblue", horiz = FALSE, width = NULL, height = NULL, border=NULL, offset=0, ...) 
{
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  if (missing(node)) 
    node <- (lastPP$Ntip + 1):length(lastPP$xx)
  XX <- lastPP$xx[node] + offset
  YY <- lastPP$yy[node]
  BOTHlabels_custom(text, node, XX, YY, adj, frame, pch, thermo, pie, 
             piecol, col, bg, horiz, width, height, border, ...)
}

BOTHlabels_custom <- function (text, sel, XX, YY, adj, frame, pch, thermo, pie, piecol, col, bg, horiz, width, height, border,...){
  if (missing(text)) 
    text <- NULL
  if (length(adj) == 1) 
    adj <- c(adj, 0.5)
  if (is.null(text) && is.null(pch) && is.null(thermo) && is.null(pie)) 
    text <- as.character(sel)
  frame <- match.arg(frame, c("rect", "circle", "none"))
  args <- list(...)
  CEX <- if ("cex" %in% names(args)) 
    args$cex
  else par("cex")
  if (frame != "none" && !is.null(text)) {
    if (frame == "rect") {
      width <- strwidth(text, units = "inches", cex = CEX)
      height <- strheight(text, units = "inches", cex = CEX)
      if ("srt" %in% names(args)) {
        args$srt <- args$srt%%360
        if (args$srt == 90 || args$srt == 270) {
          tmp <- width
          width <- height
          height <- tmp
        }
        else if (args$srt != 0) 
          warning("only right angle rotation of frame is supported;\n         try  `frame = \"n\"' instead.\n")
      }
      width <- xinch(width)
      height <- yinch(height)
      xl <- XX - width * adj[1] - xinch(0.03)
      xr <- xl + width + xinch(0.03)
      yb <- YY - height * adj[2] - yinch(0.02)
      yt <- yb + height + yinch(0.05)
      rect(xl, yb, xr, yt, col = bg)
    }
    if (frame == "circle") {
      radii <- 0.8 * apply(cbind(strheight(text, units = "inches", 
                                           cex = CEX), strwidth(text, units = "inches", 
                                                                cex = CEX)), 1, max)
      symbols(XX, YY, circles = radii, inches = max(radii), 
              add = TRUE, bg = bg)
    }
  }
  if (!is.null(thermo)) {
    parusr <- par("usr")
    if (is.null(width)) {
      width <- CEX * (parusr[2] - parusr[1])
      width <- if (horiz) 
        width/15
      else width/40
    }
    if (is.null(height)) {
      height <- CEX * (parusr[4] - parusr[3])
      height <- if (horiz) 
        height/40
      else height/15
    }
    if (is.vector(thermo) || ncol(thermo) == 1) 
      thermo <- cbind(thermo, 1 - thermo)
    thermo <- if (horiz) 
      width * thermo
    else height * thermo
    if (is.null(piecol)) 
      piecol <- rainbow(ncol(thermo))
    xl <- XX - width/2 + adj[1] - 0.5
    xr <- xl + width
    yb <- YY - height/2 + adj[2] - 0.5
    yt <- yb + height
    if (horiz) {
      rect(xl, yb, xl + thermo[, 1], yt, border = NA, col = piecol[1])
      for (i in 2:ncol(thermo)) rect(xl + rowSums(thermo[, 1:(i - 1), drop = FALSE]), yb, xl + rowSums(thermo[, 1:i]), yt, border = NA, col = piecol[i])
    }
    else {
      rect(xl, yb, xr, yb + thermo[, 1], border = NA, col = piecol[1])
      for (i in 2:ncol(thermo)) rect(xl, yb + rowSums(thermo[, 1:(i - 1), drop = FALSE]), xr, yb + rowSums(thermo[, 1:i]), border = NA, col = piecol[i])
    }
    s <- apply(thermo, 1, function(xx) any(is.na(xx)))
    xl[s] <- xr[s] <- NA
    rect(xl, yb, xr, yt, border = "black")
    if (!horiz) {
      segments(xl, YY, xl - width/5, YY)
      segments(xr, YY, xr + width/5, YY)
    }
  }
  if (!is.null(pie)) {
    if (is.data.frame(pie)) 
      pie <- as.matrix(pie)
    if (is.vector(pie) || ncol(pie) == 1) 
      pie <- cbind(pie, 1 - pie)
    xrad <- CEX * diff(par("usr")[1:2])/50
    xrad <- rep(xrad, length(sel))
    XX <- XX + adj[1] - 0.5
    YY <- YY + adj[2] - 0.5
    for (i in seq_along(sel)) {
      if (any(is.na(pie[i, ]))) 
        next
      floating.pie.asp_custom(XX[i], YY[i], pie[i, ], radius = xrad[i], 
                       col = piecol, border = border)
    }
  }
  if (!is.null(text)) 
    text(XX, YY, text, adj = adj, col = col, ...)
  if (!is.null(pch)) 
    points(XX + adj[1] - 0.5, YY + adj[2] - 0.5, pch = pch, 
           col = col, bg = bg, ...)
}

floating.pie.asp_custom <- function (xpos, ypos, x, edges = 200, radius = 1, col = NULL, startpos = 0, border, ...){
  u <- par("usr")
  user.asp <- diff(u[3:4])/diff(u[1:2])
  p <- par("pin")
  inches.asp <- p[2]/p[1]
  asp <- user.asp/inches.asp
  if (!is.numeric(x) || any(is.na(x) | x < 0)) 
    stop("floating.pie: x values must be non-negative")
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  nx <- length(dx)
  col <- if (is.null(col)) 
    rainbow(nx)
  else rep_len(col, nx)
  if (length(i <- which(dx == 1))) {
    symbols(xpos, ypos, circles = radius, inches = FALSE, 
            add = TRUE, fg = par("fg"), bg = col[i])
  }
  else {
    bc <- 2 * pi * (x[1:nx] + dx/2) + startpos
    for (i in seq_len(nx)) {
      n <- max(2, floor(edges * dx[i]))
      t2p <- 2 * pi * seq(x[i], x[i + 1], length = n) + 
        startpos
      xc <- c(cos(t2p) * radius + xpos, xpos)
      yc <- c(sin(t2p) * radius * asp + ypos, ypos)
      polygon(xc, yc, col = col[i], border=border, ...)
    }
  }
}

quickCol <- function(col, alpha){
  tmp <- c(col2rgb(col, alpha))/255
  out <- rgb(tmp[1], tmp[2], tmp[3], tmp[4])
  return(out)
}

axisPhylo_custom <- function(H){
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  
  lastPP$xx[lastPP$Ntip+1]
  lastPP$yy[lastPP$Ntip+1]
  segments(xa, ya, x, y)
  text(xa, ya, r0 - r, srt = srt, adj = c(0.5,1.1), ...)
  
}


circular.plot_custom <- function (edge, Ntip, Nnode, xx, yy, theta, r, edge.color, edge.width, edge.lty) {
  r0 <- r[edge[, 1]]
  r1 <- r[edge[, 2]]
  theta0 <- theta[edge[, 2]]
  costheta0 <- cos(theta0)
  sintheta0 <- sin(theta0)
  x0 <- r0 * costheta0
  y0 <- r0 * sintheta0
  x1 <- r1 * costheta0
  y1 <- r1 * sintheta0
  segments(x0, y0, x1, y1, col = edge.color, lwd = edge.width, 
           lty = edge.lty)
  tmp <- which(diff(edge[, 1]) != 0)
  start <- c(1, tmp + 1)
  Nedge <- dim(edge)[1]
  end <- c(tmp, Nedge)
  foo <- function(edge.feat, default) {
    if (length(edge.feat) == 1) 
      return(as.list(rep(edge.feat, Nnode)))
    edge.feat <- rep(edge.feat, length.out = Nedge)
    feat.arc <- as.list(rep(default, Nnode))
    for (k in 1:Nnode) {
      tmp <- edge.feat[start[k]]
      if (tmp == edge.feat[end[k]]) {
        feat.arc[[k]] <- tmp
      }
      else {
        if (nodedegree[k] == 2) 
          feat.arc[[k]] <- rep(c(tmp, edge.feat[end[k]]), 
                               each = 50)
      }
    }
    feat.arc
  }
  nodedegree <- tabulate(edge[, 1L])[-seq_len(Ntip)]
  co <- foo(edge.color, par("fg"))
  lw <- foo(edge.width, par("lwd"))
  ly <- foo(edge.lty, par("lty"))
  for (k in 1:Nnode) {
    i <- start[k]
    j <- end[k]
    X <- rep(r[edge[i, 1]], 100)
    Y <- seq(theta[edge[i, 2]], theta[edge[j, 2]], length.out = 100)
    x <- X * cos(Y)
    y <- X * sin(Y)
    x0 <- x[-100]
    y0 <- y[-100]
    x1 <- x[-1]
    y1 <- y[-1]
    if(any(x0 == 0) | any(x1 == 0)){
      print(k)
    }
    if(ly[[k]] != 1){
      focal_co <- rep(c("white", "black", "white"), 33)
    }else{
      focal_co <- co[[k]]
    }
    segments(x0, y0, x1, y1, col = focal_co, lwd = lw[[k]], 
             lty = ly[[k]])
  }
}

plot.phylo_custom <- function (x, type = "phylogram", use.edge.length = TRUE, node.pos = NULL, 
                               show.tip.label = TRUE, show.node.label = FALSE, edge.color = NULL, 
                               edge.width = NULL, edge.lty = NULL, node.color = NULL, node.width = NULL, 
                               node.lty = NULL, font = 3, cex = par("cex"), adj = NULL, 
                               srt = 0, no.margin = FALSE, root.edge = FALSE, label.offset = 0, 
                               underscore = FALSE, x.lim = NULL, y.lim = NULL, direction = "rightwards", 
                               lab4ut = NULL, tip.color = par("col"), plot = TRUE, rotate.tree = 0, 
                               open.angle = 0, node.depth = 1, align.tip.label = FALSE, add=FALSE,
                               ...) 
{
  Ntip <- length(x$tip.label)
  if (Ntip < 2) {
    warning("found fewer than 2 tips in the tree")
    return(NULL)
  }
  .nodeHeight <- function(edge, Nedge, yy) .C(node_height, 
                                              as.integer(edge[, 1]), as.integer(edge[, 2]), as.integer(Nedge), 
                                              as.double(yy))[[4]]
  .nodeDepth <- function(Ntip, Nnode, edge, Nedge, node.depth) .C(node_depth, 
                                                                  as.integer(Ntip), as.integer(edge[, 1]), as.integer(edge[, 
                                                                                                                           2]), as.integer(Nedge), double(Ntip + Nnode), as.integer(node.depth))[[5]]
  .nodeDepthEdgelength <- function(Ntip, Nnode, edge, Nedge, 
                                   edge.length) .C(node_depth_edgelength, as.integer(edge[, 
                                                                                          1]), as.integer(edge[, 2]), as.integer(Nedge), as.double(edge.length), 
                                                   double(Ntip + Nnode))[[5]]
  Nedge <- dim(x$edge)[1]
  Nnode <- x$Nnode
  if (any(x$edge < 1) || any(x$edge > Ntip + Nnode)) 
    stop("tree badly conformed; cannot plot. Check the edge matrix.")
  ROOT <- Ntip + 1
  type <- match.arg(type, c("phylogram", "cladogram", "fan", 
                            "unrooted", "radial", "tidy"))
  direction <- match.arg(direction, c("rightwards", "leftwards", 
                                      "upwards", "downwards"))
  if (is.null(x$edge.length)) {
    use.edge.length <- FALSE
  }
  else {
    if (use.edge.length && type != "radial") {
      tmp <- sum(is.na(x$edge.length))
      if (tmp) {
        warning(paste(tmp, "branch length(s) NA(s): branch lengths ignored in the plot"))
        use.edge.length <- FALSE
      }
    }
  }
  if (is.numeric(align.tip.label)) {
    align.tip.label.lty <- align.tip.label
    align.tip.label <- TRUE
  }
  else {
    if (align.tip.label) 
      align.tip.label.lty <- 3
  }
  if (align.tip.label) {
    if (type %in% c("unrooted", "radial") || !use.edge.length || 
        is.ultrametric(x)) 
      align.tip.label <- FALSE
  }
  if (type %in% c("unrooted", "radial") || !use.edge.length || 
      is.null(x$root.edge) || !x$root.edge) 
    root.edge <- FALSE
  phyloORclado <- type %in% c("phylogram", "cladogram", "tidy")
  horizontal <- direction %in% c("rightwards", "leftwards")
  if (type == "tidy" && any(x$edge.length < 0)) 
    stop("cannot plot in tidy mode with negative branch lengths. Check 'edge.length' vector.")
  xe <- x$edge
  if (phyloORclado) {
    phyOrder <- attr(x, "order")
    if (is.null(phyOrder) || phyOrder != "cladewise") {
      x <- reorder(x)
      if (!identical(x$edge, xe)) {
        ereorder <- match(x$edge[, 2], xe[, 2])
        if (length(edge.color) > 1) {
          edge.color <- rep(edge.color, length.out = Nedge)
          edge.color <- edge.color[ereorder]
        }
        if (length(edge.width) > 1) {
          edge.width <- rep(edge.width, length.out = Nedge)
          edge.width <- edge.width[ereorder]
        }
        if (length(edge.lty) > 1) {
          edge.lty <- rep(edge.lty, length.out = Nedge)
          edge.lty <- edge.lty[ereorder]
        }
      }
    }
    yy <- numeric(Ntip + Nnode)
    TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
    yy[TIPS] <- 1:Ntip
  }
  getStringLengthbyTip <- function(x, lab, sin, cex) {
    s <- strwidth(lab, "inches", cex = cex)
    lim <- getLimit(x, lab, sin, cex)
    alp <- lim/sin
    s * alp
  }
  getLimit <- function(x, lab, sin, cex) {
    s <- strwidth(lab, "inches", cex = cex)
    if (any(s > sin)) 
      return(1.5 * max(x))
    Limit <- 0
    while (any(x > Limit)) {
      i <- which.max(x)
      alp <- x[i]/(sin - s[i])
      Limit <- x[i] + alp * s[i]
      x <- x + alp * s
    }
    Limit
  }
  z <- reorder(x, order = "postorder")
  if (phyloORclado) {
    if (is.null(node.pos)) 
      node.pos <- if (type == "cladogram" && !use.edge.length) 
        2
    else 1
    if (node.pos == 1) {
      yy <- .nodeHeight(z$edge, Nedge, yy)
    }
    else {
      ans <- .C(node_height_clado, as.integer(Ntip), as.integer(z$edge[, 
                                                                       1]), as.integer(z$edge[, 2]), as.integer(Nedge), 
                double(Ntip + Nnode), as.double(yy))
      xx <- ans[[5]] - 1
      yy <- ans[[6]]
    }
    if (!use.edge.length) {
      if (node.pos != 2) 
        xx <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, 
                         node.depth) - 1
      xx <- max(xx) - xx
    }
    else {
      xx <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, Nedge, 
                                 z$edge.length)
    }
    if (type == "tidy") {
      if (!show.tip.label) {
        yy <- tidy.xy(z$edge, Ntip, Nnode, xx, yy)
      }
      else {
        xx.tips <- xx[1:Ntip]
        pin1 <- par("pin")[1]
        lab.strlength <- getStringLengthbyTip(xx.tips, 
                                              x$tip.label, pin1, cex)
        xx2 <- xx
        xx2[1:Ntip] <- xx2[1:Ntip] + lab.strlength
        yy <- tidy.xy(z$edge, Ntip, Nnode, xx2, yy)
      }
    }
  }
  else {
    twopi <- 2 * pi
    rotate.tree <- twopi * rotate.tree/360
    if (type != "unrooted") {
      TIPS <- x$edge[which(x$edge[, 2] <= Ntip), 2]
      xx <- seq(0, twopi * (1 - 1/Ntip) - twopi * open.angle/360, 
                length.out = Ntip)
      theta <- double(Ntip)
      theta[TIPS] <- xx
      theta <- c(theta, numeric(Nnode))
    }
    switch(type, fan = {
      theta <- .nodeHeight(z$edge, Nedge, theta)
      if (use.edge.length) {
        r <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, 
                                  Nedge, z$edge.length)
      } else {
        r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
        max_r <- max(r)
        r <- (max_r - r + 1)/max_r
      }
      theta <- theta + rotate.tree
      if (root.edge) r <- r + x$root.edge
      xx <- r * cos(theta)
      yy <- r * sin(theta)
    }, unrooted = {
      nb.sp <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
      XY <- if (use.edge.length) unrooted.xy(Ntip, Nnode, 
                                             z$edge, z$edge.length, nb.sp, rotate.tree) else unrooted.xy(Ntip, 
                                                                                                         Nnode, z$edge, rep(1, Nedge), nb.sp, rotate.tree)
      xx <- XY$M[, 1] - min(XY$M[, 1])
      yy <- XY$M[, 2] - min(XY$M[, 2])
    }, radial = {
      r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge, node.depth)
      r[r == 1] <- 0
      r <- 1 - r/Ntip
      theta <- .nodeHeight(z$edge, Nedge, theta) + rotate.tree
      xx <- r * cos(theta)
      yy <- r * sin(theta)
    })
  }
  if (phyloORclado) {
    if (!horizontal) {
      tmp <- yy
      yy <- xx
      xx <- tmp - min(tmp) + 1
    }
    if (root.edge) {
      if (direction == "rightwards") 
        xx <- xx + x$root.edge
      if (direction == "upwards") 
        yy <- yy + x$root.edge
    }
  }
  if (no.margin) 
    par(mai = rep(0, 4))
  if (show.tip.label) 
    nchar.tip.label <- nchar(x$tip.label)
  max.yy <- max(yy)
  if (is.null(x.lim)) {
    if (phyloORclado) {
      if (horizontal) {
        xx.tips <- xx[1:Ntip]
        if (show.tip.label) {
          pin1 <- par("pin")[1]
          tmp <- getLimit(xx.tips, x$tip.label, pin1, 
                          cex)
          tmp <- tmp + label.offset
        }
        else tmp <- max(xx.tips)
        x.lim <- c(0, tmp)
      }
      else {
        x.lim <- c(1, max(xx[1:Ntip]))
      }
    }
    else switch(type, fan = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
                        cex)
        x.lim <- range(xx) + c(-offset, offset)
      } else x.lim <- range(xx)
    }, unrooted = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
                        cex)
        x.lim <- c(0 - offset, max(xx) + offset)
      } else x.lim <- c(0, max(xx))
    }, radial = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.03 * cex)
        x.lim <- c(-1 - offset, 1 + offset)
      } else x.lim <- c(-1, 1)
    })
  }
  else if (length(x.lim) == 1) {
    x.lim <- c(0, x.lim)
    if (phyloORclado && !horizontal) 
      x.lim[1] <- 1
    if (type %in% c("fan", "unrooted") && show.tip.label) 
      x.lim[1] <- -max(nchar.tip.label * 0.018 * max.yy * 
                         cex)
    if (type == "radial") 
      x.lim[1] <- if (show.tip.label) 
        -1 - max(nchar.tip.label * 0.03 * cex)
    else -1
  }
  if (phyloORclado && direction == "leftwards") 
    xx <- x.lim[2] - xx
  if (is.null(y.lim)) {
    if (phyloORclado) {
      if (horizontal) {
        y.lim <- c(1, max(yy[1:Ntip]))
      }
      else {
        pin2 <- par("pin")[2]
        yy.tips <- yy[1:Ntip]
        if (show.tip.label) {
          tmp <- getLimit(yy.tips, x$tip.label, pin2, 
                          cex)
          tmp <- tmp + label.offset
        }
        else tmp <- max(yy.tips)
        y.lim <- c(0, tmp)
      }
    }
    else switch(type, fan = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
                        cex)
        y.lim <- c(min(yy) - offset, max.yy + offset)
      } else y.lim <- c(min(yy), max.yy)
    }, unrooted = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.018 * max.yy * 
                        cex)
        y.lim <- c(0 - offset, max.yy + offset)
      } else y.lim <- c(0, max.yy)
    }, radial = {
      if (show.tip.label) {
        offset <- max(nchar.tip.label * 0.03 * cex)
        y.lim <- c(-1 - offset, 1 + offset)
      } else y.lim <- c(-1, 1)
    })
  }
  else if (length(y.lim) == 1) {
    y.lim <- c(0, y.lim)
    if (phyloORclado && horizontal) 
      y.lim[1] <- 1
    if (type %in% c("fan", "unrooted") && show.tip.label) 
      y.lim[1] <- -max(nchar.tip.label * 0.018 * max.yy * 
                         cex)
    if (type == "radial") 
      y.lim[1] <- if (show.tip.label) 
        -1 - max(nchar.tip.label * 0.018 * max.yy * cex)
    else -1
  }
  if (phyloORclado && direction == "downwards") 
    yy <- y.lim[2] - yy
  if (phyloORclado && root.edge) {
    if (direction == "leftwards") 
      x.lim[2] <- x.lim[2] + x$root.edge
    if (direction == "downwards") 
      y.lim[2] <- y.lim[2] + x$root.edge
  }
  asp <- if (type %in% c("fan", "radial", "unrooted")) 
    1
  else NA
  if(!add){
    plot.default(0, type = "n", xlim = x.lim, ylim = y.lim, xlab = "", 
                 ylab = "", axes = FALSE, asp = asp, ...)
  }
  if (plot) {
    if (is.null(adj)) 
      adj <- if (phyloORclado && direction == "leftwards") 
        1
    else 0
    if (phyloORclado && show.tip.label) {
      MAXSTRING <- max(strwidth(x$tip.label, cex = cex))
      loy <- 0
      if (direction == "rightwards") {
        lox <- label.offset + MAXSTRING * 1.05 * adj
      }
      if (direction == "leftwards") {
        lox <- -label.offset - MAXSTRING * 1.05 * (1 - 
                                                     adj)
      }
      if (!horizontal) {
        psr <- par("usr")
        MAXSTRING <- MAXSTRING * 1.09 * (psr[4] - psr[3])/(psr[2] - 
                                                             psr[1])
        loy <- label.offset + MAXSTRING * 1.05 * adj
        lox <- 0
        srt <- 90 + srt
        if (direction == "downwards") {
          loy <- -loy
          srt <- 180 + srt
        }
      }
    }
    if (type %in% c("phylogram", "tidy")) {
      phylogram.plot(x$edge, Ntip, Nnode, xx, yy, horizontal, 
                     edge.color, edge.width, edge.lty, node.color, 
                     node.width, node.lty)
    }
    else {
      if (is.null(edge.color)) {
        edge.color <- par("fg")
      }
      if (is.null(edge.width)) {
        edge.width <- par("lwd")
      }
      if (is.null(edge.lty)) {
        edge.lty <- par("lty")
      }
      if (type == "fan") {
        ereorder <- match(z$edge[, 2], x$edge[, 2])
        if (length(edge.color) > 1) {
          edge.color <- rep_len(edge.color, Nedge)
          edge.color <- edge.color[ereorder]
        }
        if (length(edge.width) > 1) {
          edge.width <- rep_len(edge.width, Nedge)
          edge.width <- edge.width[ereorder]
        }
        if (length(edge.lty) > 1) {
          edge.lty <- rep_len(edge.lty, Nedge)
          edge.lty <- edge.lty[ereorder]
        }
        circular.plot_custom(z$edge, Ntip, Nnode, xx, yy, theta, 
                      r, edge.color, edge.width, edge.lty)
      }
      else cladogram.plot(x$edge, xx, yy, edge.color, edge.width, 
                          edge.lty)
    }
    if (root.edge) {
      rootcol <- if (length(edge.color) == 1) 
        edge.color
      else par("fg")
      rootw <- if (length(edge.width) == 1) 
        edge.width
      else par("lwd")
      rootlty <- if (length(edge.lty) == 1) 
        edge.lty
      else par("lty")
      if (type == "fan") {
        # rootlty <- 3
        rootcol <- "white"
        tmp <- polar2rect(x$root.edge, theta[ROOT])
        segments(0, 0, tmp$x, tmp$y, col = rootcol, lwd = rootw, 
                 lty = rootlty)
      }
      else {
        switch(direction, rightwards = segments(0, yy[ROOT], 
                                                x$root.edge, yy[ROOT], col = rootcol, lwd = rootw, 
                                                lty = rootlty), leftwards = segments(xx[ROOT], 
                                                                                     yy[ROOT], xx[ROOT] + x$root.edge, yy[ROOT], 
                                                                                     col = rootcol, lwd = rootw, lty = rootlty), 
               upwards = segments(xx[ROOT], 0, xx[ROOT], x$root.edge, 
                                  col = rootcol, lwd = rootw, lty = rootlty), 
               downwards = segments(xx[ROOT], yy[ROOT], xx[ROOT], 
                                    yy[ROOT] + x$root.edge, col = rootcol, lwd = rootw, 
                                    lty = rootlty))
      }
    }
    if (show.tip.label) {
      if (is.expression(x$tip.label)) 
        underscore <- TRUE
      if (!underscore) 
        x$tip.label <- gsub("_", " ", x$tip.label)
      if (phyloORclado) {
        if (align.tip.label) {
          xx.tmp <- switch(direction, rightwards = max(xx[1:Ntip]), 
                           leftwards = min(xx[1:Ntip]), upwards = xx[1:Ntip], 
                           downwards = xx[1:Ntip])
          yy.tmp <- switch(direction, rightwards = yy[1:Ntip], 
                           leftwards = yy[1:Ntip], upwards = max(yy[1:Ntip]), 
                           downwards = min(yy[1:Ntip]))
          segments(xx[1:Ntip], yy[1:Ntip], xx.tmp, yy.tmp, 
                   lty = align.tip.label.lty)
        }
        else {
          xx.tmp <- xx[1:Ntip]
          yy.tmp <- yy[1:Ntip]
        }
        text(xx.tmp + lox, yy.tmp + loy, x$tip.label, 
             adj = adj, font = font, srt = srt, cex = cex, 
             col = tip.color)
      }
      else {
        angle <- if (type == "unrooted") 
          XY$axe
        else atan2(yy[1:Ntip], xx[1:Ntip])
        lab4ut <- if (is.null(lab4ut)) {
          if (type == "unrooted") 
            "horizontal"
          else "axial"
        }
        else match.arg(lab4ut, c("horizontal", "axial"))
        xx.tips <- xx[1:Ntip]
        yy.tips <- yy[1:Ntip]
        if (label.offset) {
          xx.tips <- xx.tips + label.offset * cos(angle)
          yy.tips <- yy.tips + label.offset * sin(angle)
        }
        if (lab4ut == "horizontal") {
          y.adj <- x.adj <- numeric(Ntip)
          sel <- abs(angle) > 0.75 * pi
          x.adj[sel] <- -strwidth(x$tip.label)[sel] * 
            1.05
          sel <- abs(angle) > pi/4 & abs(angle) < 0.75 * 
            pi
          x.adj[sel] <- -strwidth(x$tip.label)[sel] * 
            (2 * abs(angle)[sel]/pi - 0.5)
          sel <- angle > pi/4 & angle < 0.75 * pi
          y.adj[sel] <- strheight(x$tip.label)[sel]/2
          sel <- angle < -pi/4 & angle > -0.75 * pi
          y.adj[sel] <- -strheight(x$tip.label)[sel] * 
            0.75
          text(xx.tips + x.adj * cex, yy.tips + y.adj * 
                 cex, x$tip.label, adj = c(adj, 0), font = font, 
               srt = srt, cex = cex, col = tip.color)
        }
        else {
          if (align.tip.label) {
            POL <- rect2polar(xx.tips, yy.tips)
            POL$r[] <- max(POL$r)
            REC <- polar2rect(POL$r, POL$angle)
            xx.tips <- REC$x
            yy.tips <- REC$y
            segments(xx[1:Ntip], yy[1:Ntip], xx.tips, 
                     yy.tips, lty = align.tip.label.lty)
          }
          if (type == "unrooted") {
            adj <- abs(angle) > pi/2
            angle <- angle * 180/pi
            angle[adj] <- angle[adj] - 180
            adj <- as.numeric(adj)
          }
          else {
            s <- xx.tips < 0
            angle <- angle * 180/pi
            angle[s] <- angle[s] + 180
            adj <- as.numeric(s)
          }
          font <- rep(font, length.out = Ntip)
          tip.color <- rep(tip.color, length.out = Ntip)
          cex <- rep(cex, length.out = Ntip)
          for (i in 1:Ntip) text(xx.tips[i], yy.tips[i], 
                                 x$tip.label[i], font = font[i], cex = cex[i], 
                                 srt = angle[i], adj = adj[i], col = tip.color[i])
        }
      }
    }
    if (show.node.label) 
      text(xx[ROOT:length(xx)] + label.offset, yy[ROOT:length(yy)], 
           x$node.label, adj = adj, font = font, srt = srt, 
           cex = cex)
  }
  L <- list(type = type, use.edge.length = use.edge.length, 
            node.pos = node.pos, node.depth = node.depth, show.tip.label = show.tip.label, 
            show.node.label = show.node.label, font = font, cex = cex, 
            adj = adj, srt = srt, no.margin = no.margin, label.offset = label.offset, 
            x.lim = x.lim, y.lim = y.lim, direction = direction, 
            tip.color = tip.color, Ntip = Ntip, Nnode = Nnode, root.time = x$root.time, 
            align.tip.label = align.tip.label)
  assign("last_plot.phylo", c(L, list(edge = xe, xx = xx, yy = yy)), 
         envir = .PlotPhyloEnv)
  invisible(L)
}


draw_circle <- function(radius, col) {
  theta <- seq(0, 2*pi, length.out=100) # Create 100 points along the circle
  x <- radius * cos(theta)
  y <- radius * sin(theta)
  polygon(x, y, border=NA, col=col) # Draw filled circle with no border
}

getPolygonCoords <- function(tip_coord, inner_ratio = 1.01, outer_ratio = 1.05, index=NULL, return.df=FALSE, order=FALSE){
  if(!is.null(index)){
    tip_coord <- tip_coord[index,]
  }
  center_x <- 0
  center_y <- 0
  r <- sqrt((tip_coord$x - center_x)^2 + (tip_coord$y - center_y)^2)
  theta <- atan2(tip_coord$y - center_y, tip_coord$x - center_x)
  inner_r <- r * inner_ratio
  outer_r <- r * outer_ratio
  extended_x_inner <- center_x + inner_r * cos(theta)
  extended_x_outer <- center_x + outer_r * cos(theta)
  extended_y_inner <- center_y + inner_r * sin(theta)
  extended_y_outer <- center_y + outer_r * sin(theta)
  if(return.df){
    xy_df <- data.frame(x0=extended_x_inner, 
                        y0=extended_y_inner, 
                        x1=extended_x_outer,
                        y1=extended_y_outer)
    return(xy_df)
  }
  if(order){
    ixi <- order(extended_x_inner)
    ixo <- order(extended_x_outer)
    vertices_x <- c(extended_x_inner[ixi], rev(extended_x_outer[ixo]))
    vertices_y <- c(extended_y_inner[ixi], rev(extended_y_outer[ixo]))
  }else{
    vertices_x <- c(extended_x_inner, rev(extended_x_outer))
    vertices_y <- c(extended_y_inner, rev(extended_y_outer))
  }
  poly_coords <- data.frame(vertices_x=vertices_x, vertices_y=vertices_y)
  # poly_coords <- rbind(poly_coords, poly_coords[1,])
  return(poly_coords)
}

find_centroid <- function(x, y){
  Cx <- Cy <- A <- c()
  for(i in 1:(length(x)-1)){
    Cx <- c(Cx, (x[i] + x[i+1]) * ((x[i]*y[i+1]) - (x[i+1]*y[i])))
    Cy <- c(Cy, (y[i] + y[i+1]) * ((x[i]*y[i+1]) - (x[i+1]*y[i])))
    A <- c(A, ((x[i]*y[i+1]) - (x[i+1]*y[i])))
  }
  A <- sum(A)/2
  Cx <- sum(Cx)/(6*A)
  Cy <- sum(Cy)/(6*A)
  return(c(Cx, Cy))
}

rescale_vector <- function(v, min_new, max_new) {
  min_old <- min(v)
  max_old <- max(v)
  v_rescaled <- min_new + (v - min_old) * (max_new - min_new) / (max_old - min_old)
  return(v_rescaled)
}


make_less_vibrant <- function(colors, saturation_decrease = 0.5, brightness_decrease = 0.9) {
  a <- col2rgb(colors) / 255
  b <- apply(a, 2, function(x) 
    rgb2hsv(x[1], x[2], x[3], 1))
  c <- apply(b, 2, function(x) 
    hsv(x[1], x[2]*saturation_decrease, x[3]*brightness_decrease))
  return(c)
}
