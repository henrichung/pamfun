# Set working directory
#setwd("E:/Projects/pamfun")
# Load necessary Libraries
library(data.table)
library(tidyverse)
library(dtplyr)
library(Matrix)
library(parallel)
library(vegan)
#####
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
args = commandArgs(trailingOnly=TRUE)

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

phylum_s = args[1]
phylum_s = "Firmicutes"
bact <- filter(taxa_ref, phylum == phylum_s) %>%
  group_by(species) %>%
  slice_sample(n = 5) %>%
  pull(organism_taxid)

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
binDTM <- result %>%
  column_to_rownames("assembly_accession") %>%
  as.matrix() %>%
  t() %>%
  as.data.frame() %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  as.matrix()


fusion_data2 <- mydata2 %>% 
  filter(organism_taxid %in% bact) %>%
  select(fusion_lvl_1, protein_name) %>%
  unique() %>%
  as.data.table()


message(Sys.time(), " SPLITTING DATAFRAME")
a <- binDTM %>%
  as.data.frame() %>%
  rownames_to_column("fusion_lvl_1") %>%
  group_by(across(c(-fusion_lvl_1))) %>%
  group_split() 

message(Sys.time(), " PULLING FUSION LABELS")
b <- lapply(a, function(.x){select(.x, "fusion_lvl_1") %>% pull("fusion_lvl_1")}) 


message(Sys.time(), " SUBSETTING FUSION DATA")
c <- lapply(b, function(spec){
  res <- fusion_data %>%
    filter(fusion_lvl_1 %in% spec) %>%
    select(fusion_lvl_1, protein_name) %>% 
    mutate(protein_name = tolower(protein_name)) %>%
    group_by(fusion_lvl_1, protein_name) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    arrange(fusion_lvl_1, protein_name) %>%
    as_tibble()
  return(res)
  })

message(Sys.time(), " WRITING TO FILE")
d <- bind_rows(c, .id = "group") %>%
  filter(protein_name != "") %>% 
  group_by(group) %>%
  arrange(desc(n), .by_group = TRUE) %>%
  ungroup() %>%
  mutate(protein_name = paste0(protein_name, " (", n, ")")) %>%
  group_by(group, fusion_lvl_1) %>%
  summarize(
    protein_names = paste0(protein_name, collapse = " | "),
    protein_names_n = sum(n)) %>%
  mutate(group = as.numeric(group)) %>%
  ungroup() %>% group_by(group) %>%
  mutate(fusion_lvl_1_n = n()) %>%
  ungroup() %>%
  arrange(group, desc(protein_names_n)) %>%
  select(group, fusion_lvl_1, fusion_lvl_1_n, protein_names, protein_names_n) %>%
  arrange(desc(fusion_lvl_1_n))

write_csv(d, paste(phylum_s, ".csv", sep = ""))
saveRDS(d, "proteo_all.RDS")
