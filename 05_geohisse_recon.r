# rm(list=ls())
setwd("~/biome_shifts/")

library(ape)
library(hisse)
library(parallel)

# finished model sets for particular datsets
to_load <- dir("pilot_results/", full.names = TRUE)

load(to_load[1])
model_res <- res[[5]]

individual_recon <- function(model_res){
        if(class(model_res)[1] == "try-error"){
                recon <- NA
        }else{
                hidden.states <- ncol(model_res$trans.matrix )/3
                recon <- MarginReconGeoSSE(phy = model_res$phy, data = model_res$data, hidden.states = hidden.states, f = model_res$f, pars = model_res$solution, root.type = model_res$root.type, root.p = model_res$root.p, n.cores = 4, assume.cladogenetic = FALSE, get.tips.only = TRUE)
        }
        return(recon)
}

for(i in 1:length(to_load)){
        load(to_load[i])
        recon <- mclapply(res, individual_recon, mc.cores=26)
        file_name <- gsub("pilot_results_", "recon_", to_load[i])
        save(recon, file=file_name)
}