library(parallel)
library(tidyverse)
library(vegan)
library(data.table)
#####

message(Sys.time(), " READING DATA")
binDTM <- readRDS("binDTM_nosingletons_proteo.RDS") 

message(Sys.time(), " READING TAXONOMIC DATA")
taxa_ref <- readr::read_csv("data/taxa_reference.csv")
bact <- filter(taxa_ref, phylum == "Proteobacteria") %>% pull(organism_taxid)

message(Sys.time(), " READING FUSION DATA")
fusion_data <- data.table::fread("data/fusion_data.tsv")
mydata = fusion_data
mydata2 <- dtplyr::lazy_dt(mydata)

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


message(Sys.time(), " SUBSETTING FUSION DATA DATA")
c <- lapply(b, function(.x){
  res <- fusion_data %>%
    filter(fusion_lvl_1 %in% .x) %>%
    select(fusion_lvl_1, protein_name) %>% 
    mutate(protein_name = tolower(protein_name)) %>%
    unique() %>%
    arrange(fusion_lvl_1, protein_name) %>%
    as_tibble()
  return(res)
  })

message(Sys.time(), " WRITING TO FILE")
d <- bind_rows(c, .id = "group")
saveRDS(d, "proteo_all.RDS")
#write_csv(d, "proteo_all.csv")