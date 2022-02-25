## ---------------------------
## Purpose of script: Subset the fusion data to relevant clusters and convert into binary matrix
## Author: Henri Chung
## Date Created: 2021-01-12
## ---------------------------
## Notes:
## ---------------------------
library(data.table)
library(tidyverse)
library(dtplyr)
library(Matrix)
library(parallel)
library(vegan)
#####
rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
N = 1L
k = 5

message(Sys.time(), ": Number of cores: ", N)
message(Sys.time(), ": Reading in data")

subset_data <- read_tsv("data/0_fusion_data_subset.tsv")
fusion_to_ec <- read_csv("outputs/0_fusion_to_ec.csv")
cast_data <- subset_data %>%
  select(assembly_accession, fusion_lvl_1) %>%
  filter(fusion_lvl_1 %in% fusion_to_ec$fusion_lvl_1) %>%
  group_by(fusion_lvl_1) %>%
  mutate(val = 1) %>%
  as.data.table()

# Split data
message(Sys.time(), ": Finding interval.")
cast_data[ , cast_cat := findInterval(fusion_lvl_1, seq(1, 800000, 10000))]
message(Sys.time(), ": Splitting by interval.")
w_list <- split(cast_data , by = 'cast_cat' )
w_list <- mclapply( w_list , function( x ) x[ , cast_cat := NULL ] , mc.cores = N)
message(Sys.time(), ": Casting each element into the list as wide.")
w_list <- mclapply( w_list , function( z ) data.table::dcast( z , assembly_accession ~ fusion_lvl_1 , value.var = 'val' , drop = FALSE, fill = 0)  , mc.cores = N)
message(Sys.time(), ": Attempting to join the list")
result <- Reduce( function( ... ) merge( ... , by = 'assembly_accession' , all = TRUE ) , w_list )

# Cast into wide format
message(Sys.time(), ": Casting into wide format.")
binDTM <- result %>%
  column_to_rownames("assembly_accession") %>%
  as.matrix() %>%
  t() %>%
  as.data.frame() %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  as.matrix()

# Binarize
binDTM[binDTM > 1] <- 1
# filter our rows with <= 5 samples.
binDTM <- binDTM[rowSums(binDTM) > k,]

# save output
saveRDS(binDTM, "outputs/1_subsample_binDTM.RDS")
