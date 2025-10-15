# rm(list=ls())
setwd("~/biome_shifts/")

library(ape)
library(hisse)
library(parallel)

# finished model sets for particular datsets
to_load <- dir("5_results/", full.names = TRUE)

quick_check <- function(model_path){
  load(model_path)
  return(all(unlist(lapply(res, function(x) class(x) == "try-error"))))
}

failed_load <- to_load[sapply(to_load, quick_check)]
to_load <- to_load[!sapply(to_load, quick_check)]

nm_idx1 <- to_load[grep("idx1.RData", to_load)]
all_res <- vector("list", length(nm_idx1))
model_table_list <- vector("list", length(nm_idx1))

for(i in seq_along(nm_idx1)){
  clade_name <- sub("_idx1\\.RData$","", basename(nm_idx1[i]))
  files <- to_load[grepl(paste0("^", clade_name, "_idx[0-9]+\\.RData$"), basename(to_load))]
  if(length(files) == 0) next
  idxs <- as.integer(sub(".*_idx([0-9]+)\\.RData$","\\1", files))
  ord <- order(idxs)
  files <- files[ord]; idxs <- idxs[ord]
  clade_res <- lapply(files, function(f){ load(f); res })
  names(clade_res) <- paste0("M", idxs)
  all_res[[i]] <- clade_res
  model_table_list[[i]] <- setNames(GetAICWeights(clade_res), names(clade_res))
}

# combine AIC weight vectors into a rectangular data.frame (NA where a model is missing)
all_models <- sort(unique(unlist(lapply(model_table_list, names))))
model_table <- t(sapply(model_table_list, function(x){
  v <- setNames(rep(NA_real_, length(all_models)), all_models)
  if(!is.null(x)) v[names(x)] <- x
  v
}))
rownames(model_table) <- sub("_idx1\\.RData$","", basename(nm_idx1))
rownames(model_table) <- gsub("res_state", "", rownames(model_table))
model_table <- as.data.frame(model_table)
model_table <- model_table[,paste0("M", 1:36)]

to_run <- c()
best_models <- apply(model_table, 1, which.max)
for(i in 1:length(best_models)){
  focal_id <- paste0("idx", best_models[i], ".RData")
  focal_clade <- names(best_models)[i]
  tmp_load <- to_load[grep(focal_clade, to_load)]
  final_load <- tmp_load[grep(focal_id, tmp_load)]
  to_run <- c(to_run, final_load)
}

individual_recon <- function(model_res){
  if(class(model_res)[1] == "try-error"){
    recon <- NA
  }else{
    hidden.states <- ncol(model_res$trans.matrix )/4
    recon <- MarginReconMuHiSSE(phy = model_res$phy, data = model_res$data, hidden.states = hidden.states, f = model_res$f, pars = model_res$solution, root.type = model_res$root.type, root.p = model_res$root.p, get.tips.only = FALSE)
  }
  return(recon)
}

run_recon_path <- function(to_load){
  load(to_load)
  file_name <- gsub("5_results/", "6_recons/", to_load)
  file_name <- gsub("res_state", "recon", file_name)
  recon <- try(individual_recon(res))
  save(recon, file=file_name)
  return(NULL)
}

num_cores <- 13
mclapply(to_run, run_recon_path, mc.cores = num_cores)


