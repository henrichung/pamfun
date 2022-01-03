library(tidyverse)
library(dtplyr)
library(data.table)
library(parallel)

binDTM <- readRDS("subsample_binDTM.RDS")

taxa_ref <- readr::read_csv("data/taxa_reference.csv")
subset_data <- read_tsv("fusion_data_subset.tsv")
subset_groups <- readr::read_csv("subset_groups.csv")
group_no = 11987


group_ind <- subset_groups %>%
	arrange(desc(fusion_lvl_1_n)) %>%
	filter(group == group_no) %>%
	pull(fusion_lvl_1) %>%
	as.character()

binDTM_int <- binDTM[group_ind, ]

cols <- colSums(binDTM_int)

orgs <- names(which(cols > 1))

group_int <- subset_groups %>% 
	filter(group == group_no)
write_csv(group_int, "group_11987.csv")

taxa_int <- filter(taxa_ref, assembly_accession %in% orgs)
write_csv(taxa_int, "taxa_11987.csv")

subset_int <- subset_data %>%
	filter(fusion_lvl_1 %in% group_ind)
write_csv(subset_int, "fusion_11987.csv")