library(data.table)
library(tidyverse)
library(dtplyr)
library(Matrix)
library(parallel)
library(vegan)
#####
rm(list = ls())

binDTM <- readRDS("outputs/subsample_binDTM.RDS")
subset_data <- data.table::fread("data/0_fusion_data_subset.tsv")
# Parameters
N = 4

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

message(Sys.time(), " DESCRIBING FILE")
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

saveRDS(grouped_fusion, "outputs/1_grouped_fusion_all.RDS")

grouped_fusion_over1 <- grouped_fusion %>%
  filter(fusion_lvl_1_n > 1)
write_csv(grouped_fusion_over1, "outputs/1_subset_groups.csv")

grouped_fusion_singletons <- grouped_fusion %>%
  filter(fusion_lvl_1_n <= 1)
write_csv(grouped_fusion_singletons, "outputs/1_subset_groups_singletons.csv")



