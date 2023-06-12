# rm(list=ls())
setwd("~/biome_shifts/")

library(ape)
library(hisse)
library(parallel)
require(dplyr)

# finished model sets for particular datsets
clades <- dir("5_results/") %>% gsub("results_", "", .) %>% gsub(".RData", "", .) %>% unique(.)
to_load_results <- dir("5_results/", full.names = TRUE)
to_load_model_tables <- dir("tables/", full.names = TRUE)

model_table <- read.csv(to_load_model_tables[1])
load(to_load_results[1])
