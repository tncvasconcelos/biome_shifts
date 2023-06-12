# rm(list=ls())
setwd("~/biome_shifts/")

library(ape)
library(hisse)
library(parallel)
require(dplyr)

# finished model sets for particular datsets
clades <- dir("5_results/") %>% gsub("results_", "", .) %>% gsub(".RData", "", .) %>% unique(.)
to_load_recon <- dir("6_recons/", full.names = TRUE)
to_load_results <- dir("5_results/", full.names = TRUE)

# function for pulling out AIC and stuff like that
getModelRes <- function(model_res){
  if(class(model_res)[1] == "try-error"){
      out <- c(lnLik=NA, k=NA, AIC=NA)
  }else{
      out <- c(lnLik=model_res$loglik, k=max(model_res$index.par, na.rm=TRUE), AIC=model_res$AIC)
  }
  return(out)
}

getModelTable <- function(model_list){
  model_table <- as.data.frame(do.call(rbind, lapply(model_list, getModelRes)))
  model_table$dAIC <- model_table$AIC - min(model_table$AIC, na.rm=TRUE)
  model_table$AICwt <- exp(-0.5 * model_table$dAIC)/sum(exp(-0.5 * model_table$dAIC), na.rm=TRUE)
  return(model_table)
}

tables <- list()
for(i in 1:length(to_load_results)){
  load(to_load_results[i])
  tables[[i]] <- getModelTable(res)
  write.csv(tables[[i]], file=paste0("tables/", clades[i], "_model_table.csv"), row.names=FALSE)
}
names(tables) <- clades

