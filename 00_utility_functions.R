
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
library(BioGeoBEARS)


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
