# Set working directory
#setwd("E:/Projects/pamfun")
# Load necessary Libraries
library(data.table)
library(tidyverse)
library(dtplyr)
library(Matrix)
library(parallel)
rm(list = ls())
#Subset larger data into working set (once)
# Sample 1mil rows
#ind <- sample(1:nrow(fusion_data), 1000000)
#working <-  fusion_data[ind,]
#readr::write_tsv(working, "fusion_data_working2.tsv")
#rm(fusion_data)
# Read in working data

#acc_to_pfam <- fread("data/fusion_pfam_scan_matches_domain_arrangement_preserved.tsv")
#seguid_to_ec <- fread("data/uniprot-pe-exp_seguid_to_ec_mapping.tsv")
#taxa_ref <- readr::read_csv("outputs/taxa_reference.csv")
#fusion_data <- data.table::fread("data/fusion_data.tsv")
#mydata = fusion_data
N = 16L
message("Number of cores: ", N)
message("Reading in data")
fusion_data <- data.table::fread("data/fusion_data.tsv")
taxa_ref <- readr::read_csv("data/taxa_reference.csv")
mydata = fusion_data

#mydata <- readr::read_tsv("data/fusion_data_working2.tsv")
mydata2 <- dtplyr::lazy_dt(mydata)

message("Removing Singletons")
# select assembly accession and fusion lvl 1 clusters
# remove singletons

d0 <- mydata2 %>%
	select(fusion_lvl_1, ncbi_accession) %>%
	group_by(fusion_lvl_1) %>%
	summarize(n = n()) %>%
	filter(n > 10) %>%
	pull(fusion_lvl_1)

# remove singletons and clusters <= 10
# filter to single phylum
# filter to just Proteobacteria

bact <- filter(taxa_ref, phylum == "Cyanobacteria") %>% pull(organism_taxid)
d1 <- mydata2 %>% 
  filter(organism_taxid %in% bact) %>%
  select(assembly_accession, fusion_lvl_1) %>%
  filter(fusion_lvl_1 %in% d0) %>%
  filter(!grepl("s|S", fusion_lvl_1)) %>%
  group_by(fusion_lvl_1) %>%
  mutate(val = 1) %>%
  as.data.table()


message("Dimensions of data:", print(dim(d1)))
 # three million x one thousand

message("Finding interval")
d1[ , cast_cat := findInterval(fusion_lvl_1, seq(10000 , 200000 , 10000) ) ]
message("Splitting by interval")
w_list <- split(d1 , by = 'cast_cat' )
w_list <- mclapply( w_list , function( x ) x[ , cast_cat := NULL ] , mc.cores = N)
message("casting each element into the list as wide")
w_list <- mclapply( w_list , function( z ) data.table::dcast( z , assembly_accession ~ fusion_lvl_1 , value.var = 'val' , drop = FALSE, fill = 0)  , mc.cores = N)
message("attempting to join the list")
result <- Reduce( function( ... ) merge( ... , by = 'assembly_accession' , all = TRUE ) , w_list )

message("casting into wide format")
# cast into wide format.
d2 <- result %>%
  column_to_rownames("assembly_accession") %>%
  as.matrix() %>%
  t()

saveRDS(d2, "binDTM_nosingletons_proteo.RDS")


