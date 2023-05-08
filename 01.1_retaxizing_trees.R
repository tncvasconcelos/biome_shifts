

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

source("00_utility_functions.R")

all_trees <- load.trees("2_trees")

for(tree_index in 1:length(all_trees)){
  all_trees[[tree_index]]$tip.label <- resolve.names(all_trees[[tree_index]]$tip.label)
  one_tree <- all_trees[[tree_index]]
  save(one_tree, file=paste0("2.1_trees/",names(all_trees)[tree_index],".Rsave"))
}

