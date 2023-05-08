
resolve.names <- function(names_to_solve) {
  gnr_resolve_x <- function(x) {
    sources <- taxize::gnr_datasources()
    tmp.name <- suppressWarnings(taxize::gnr_resolve(names=x, data_source_ids=sources$id[sources$title == "GBIF Backbone Taxonomy"], best_match_only=TRUE)$matched_name)
    if(is.null(tmp.name)) {
      tmp.name <- paste0(x,"_UNMATCHED")
    }
    return(tmp.name)
  }
  all_names <- parallel::mclapply(names_to_solve, gnr_resolve_x, mc.cores=6)
  return(as.character(all_names))
}

species_list <- readRDS("species_list.Rdata")

taxized_names <- resolve.names(species_list[1:1000]) # This function adjust the names to the GBIF taxonomic backbone
reference_table <- data.frame(wcvp_name = species_list[1:100000], gbif_name = taxized_names) # you will need this table later
write.csv(reference_table, file="reference_table1.csv", row.names = F) # saving table that you will need later

taxized_names <- resolve.names(species_list[100001:200000]) # This function adjust the names to the GBIF taxonomic backbone
reference_table <- data.frame(wcvp_name = species_list[100001:200000], gbif_name = taxized_names) # you will need this table later
write.csv(reference_table, file="reference_table2.csv", row.names = F) # saving table that you will need later

taxized_names <- resolve.names(species_list[200001:300000]) # This function adjust the names to the GBIF taxonomic backbone
reference_table <- data.frame(wcvp_name = species_list[200001:300000], gbif_name = taxized_names) # you will need this table later
write.csv(reference_table, file="reference_table3.csv", row.names = F) # saving table that you will need later

taxized_names <- resolve.names(species_list[300000:400000]) # This function adjust the names to the GBIF taxonomic backbone
reference_table <- data.frame(wcvp_name = species_list[300000:400000], gbif_name = taxized_names) # you will need this table later
write.csv(reference_table, file="reference_table4.csv", row.names = F) # saving table that you will need later

taxized_names <- resolve.names(species_list[400000:length(species_list)]) # This function adjust the names to the GBIF taxonomic backbone
reference_table <- data.frame(wcvp_name = species_list[400000:length(species_list)], gbif_name = taxized_names) # you will need this table later
write.csv(reference_table, file="reference_table5.csv", row.names = F) # saving table that you will need later
