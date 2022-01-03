library(data.table)
library(tidyverse)
library(dtplyr)
library(Matrix)
library(parallel)
library(vegan)
#####
rm(list = ls())
args = commandArgs(trailingOnly=TRUE)

N = 12L
message(Sys.time(), ": Number of cores: ", N)
message(Sys.time(), ": Reading in data")
fusion_data <- data.table::fread("data/fusion_data.tsv")
taxa_ref <- readr::read_csv("data/taxa_reference.csv")

#subset data to assembly accessions (species) so that there are only 5 per species
#subset to top 50,000 functions.
message(Sys.time(), " Subsetting Dataset")
subset <- readRDS("data/subsamples.RDS")

#where subset is a list object with function_ids and assembly_accessions
assembly <- paste0("GCA_", unlist(subset$assembly))

mydata <- dtplyr::lazy_dt(fusion_data)

subset_data <- mydata %>%
	mutate(assembly_accession2 = gsub("\\.[0-9]", "", assembly_accession)) %>%
	filter(assembly_accession2 %in% assembly) %>%
	select(-assembly_accession2) %>%
	filter(fusion_lvl_1 %in% unlist(subset$functions)) %>%
	as.data.table()
write_tsv(subset_data, "fusion_data_subset.tsv")
rm(fusion_data, mydata); gc()


cast_data <- subset_data %>%
  select(assembly_accession, fusion_lvl_1) %>%
  group_by(fusion_lvl_1) %>%
  mutate(val = 1) %>%
  as.data.table()

message(Sys.time(), ": Finding interval.")
cast_data[ , cast_cat := findInterval(fusion_lvl_1, seq(10000 , 200000 , 10000) ) ]
message(Sys.time(), ": Splitting by interval.")
w_list <- split(cast_data , by = 'cast_cat' )
w_list <- mclapply( w_list , function( x ) x[ , cast_cat := NULL ] , mc.cores = N)
message(Sys.time(), ": Casting each element into the list as wide.")
w_list <- mclapply( w_list , function( z ) data.table::dcast( z , assembly_accession ~ fusion_lvl_1 , value.var = 'val' , drop = FALSE, fill = 0)  , mc.cores = N)
message(Sys.time(), ": Attempting to join the list")
result <- Reduce( function( ... ) merge( ... , by = 'assembly_accession' , all = TRUE ) , w_list )

message(Sys.time(), ": Casting into wide format.")
# cast into wide format.
binDTM <- result %>%
  column_to_rownames("assembly_accession") %>%
  as.matrix() %>%
  t() %>%
  as.data.frame() %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  as.matrix()

binDTM[binDTM > 1] <- 1
saveRDS(binDTM, "subsample_binDTM.RDS")

binDTM <- readRDS("subsample_binDTM.RDS")

message(Sys.time(), " SPLITTING DATAFRAME")
split_df <- binDTM %>%
  as.data.frame() %>%
  rownames_to_column("fusion_lvl_1") %>%
  group_by(across(c(-fusion_lvl_1))) %>%
  group_split() 

message(Sys.time(), " PULLING FUSION LABELS")
labels <- mclapply(split_df, function(.x){select(.x, "fusion_lvl_1") %>% pull("fusion_lvl_1")}, mc.cores = N) 


message(Sys.time(), " SUBSETTING FUSION DATA")
split_fusion <- mclapply(labels, function(spec){
  res <- subset_data %>%
    filter(fusion_lvl_1 %in% spec) %>%
    select(fusion_lvl_1, protein_name) %>% 
    mutate(protein_name = tolower(protein_name)) %>%
    group_by(fusion_lvl_1, protein_name) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    arrange(fusion_lvl_1, protein_name) %>%
    as_tibble()
  return(res)
  }, mc.cores = N) 

message(Sys.time(), " WRITING TO FILE")
grouped_fusion <- bind_rows(split_fusion, .id = "group") %>%
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


#write_csv(grouped_fusion, "subset_groups.csv")
saveRDS(grouped_fusion, "grouped_fusion_all.RDS")


grouped_fusion <- readRDS("grouped_fusion_all.RDS")

grouped_fusion_over1 <- grouped_fusion %>%
	filter(fusion_lvl_1_n > 1)
write_csv(grouped_fusion_over1, "subset_groups.csv")

grouped_fusion_singletons <- grouped_fusion %>%
	filter(fusion_lvl_1_n <= 1)
write_csv(grouped_fusion_singletons, "subset_groups_singletons.csv")



