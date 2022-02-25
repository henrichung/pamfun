## ---------------------------
## Purpose of script: index of subset fusion dataset with top 50,000 fusion clusters and treemer 
## Author: Henri Chung
## Date Created: 2021-01-12
## ---------------------------
## Notes:
## ---------------------------

# Load necessary libraries
library(data.table)
library(tidyverse)
library(dtplyr)

# Read in fusion data
message(Sys.time(), ": Number of cores: ", N)
message(Sys.time(), ": Reading in data")
fusion_data <- data.table::fread("data/fusion_data.tsv")

# read in treemer trimmed representative assemblies
assemblies <- readr::read_tsv("data/balanced_organism_set.treemmer-RTL-0.9.gtdb.incl_taxonomy.tsv") %>%
  mutate(assembly_accession = paste("GCA_", assembly_id, sep = "")) %>%
  pull(assembly_accession)

# read in top 50,000 fusion clusters
fusion_clusters <- read.table("data/final_fusion_datafile.reduced_organism_set.treemmer.gtdb.rand_50k.txt")

# subset
mydata <- dtplyr::lazy_dt(fusion_data)

subset_data <- mydata %>%
  mutate(assembly_accession2 = gsub("\\.[0-9]", "", assembly_accession)) %>% # rename assembly accessions to remove version number
  filter(assembly_accession2 %in% assemblies) %>% # filter to assemblies in subset
  select(-assembly_accession2) %>% # remove renamed assembly column
  filter(fusion_lvl_1 %in% fusion_clusters$V1) %>% # filter to fusion clusters in subset
  as.data.table()
rm(fusion_data, mydata); gc()

# write to file.
write_tsv(subset_data, "data/0_fusion_data_subset.tsv")

