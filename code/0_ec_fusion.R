## ---------------------------
## Purpose of script: Adjust the dataset to only fusion datasets with EC annotations using an adjustable HFSP parameter.
## Author: Henri Chung
## Date Created: 2021-01-12
## ---------------------------
## Notes:
## ---------------------------

# Load necessary libraries
library(data.table)
library(tidyverse)
rm(list = ls())
hfsp_threshold = 20
ec_threshold = 0.1

# read in fusion to hfsp mapping
seguid_to_ec <- read_tsv("data/uniprot-pe1-exp-ec.extended_by_hfsp.mapping") %>%
  filter(hfsp > hfsp_threshold) %>%
  mutate(ec_number3 = gsub("\\.[0-9]+\\Z", "", ec_number, perl = TRUE))
 head(seguid_to_ec)

# read in fusion data
subset_data <- data.table::fread("data/0_fusion_data_subset.tsv")
mydata <- dtplyr::lazy_dt(subset_data)

# Join fusion_to_ec and seguid_to_ec by seguid
fusion_to_ec <- mydata %>%
  select(fusion_lvl_1, seguid) %>%
  unique() %>% 
  left_join(seguid_to_ec) %>%
  filter(!is.na(ec_number)) %>% 
  select(seguid, fusion_lvl_1, hfsp, ec_number, ec_number3) %>%
  as.data.table() 

# Filter out EC to fusion annotations 
# Filter out EC annotations to fusion clusters where the EC annotations,
# is under a threshold of the total proportion of EC annotations in a cluster.

# Filtering proportions of EC numbers occurs first.
# then truncation of 4th EC number.
filtered_fusion_to_ec <- fusion_to_ec %>%
	group_by(fusion_lvl_1, ec_number) %>%
	summarize(n = n()) %>%
	ungroup() %>% 
	group_by(fusion_lvl_1) %>%
	summarize(prop = n / sum(n), ec_number = ec_number) %>%
	filter(prop > ec_threshold) %>%
	select(fusion_lvl_1, ec_number) %>%
	left_join(select(fusion_to_ec, fusion_lvl_1, ec_number3)) %>%
	unique() %>%
	as.data.table()
write_csv(filtered_fusion_to_ec, "outputs/0_fusion_to_ec.csv")

# Calculate how many fusion clusters have identical EC numbers.
ec_duplicates <- fusion_to_ec2 %>%
	group_by(ec_number) %>%
	filter(n() > 1) %>%
	as.data.table()
head(ec_duplicates)
length(unique(ec_duplicates$ec_number)) #291
length(unique(ec_duplicates$fusion_lvl_1)) #718
write_csv(ec_duplicates, "outputs/0_ec_duplicates.csv")

# Calculate how many EC numbers are assigned to multiple fusion clusters.
fusion_duplicates <- fusion_to_ec2 %>% 
	group_by(fusion_lvl_1) %>% 
	filter(n() > 1) %>% 
	as.data.table()
head(fusion_duplicates)
length(unique(fusion_duplicates$ec_number)) #731
length(unique(fusion_duplicates$fusion_lvl_1)) #280
write_csv(fusion_duplicates, "outputs/0_fusion_duplicates.csv")