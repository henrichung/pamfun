library(data.table)
library(tidyverse)
library(dtplyr)
library(Matrix)
library(parallel)
library(vegan)
library(pvclust)
#####
rm(list = ls())
#mydata = fusion_data
args = commandArgs(trailingOnly=TRUE)

N = 36L
message(Sys.time(), " Number of cores: ", N)
message(Sys.time(), " Reading in data")
fusion_data <- data.table::fread("data/fusion_data.tsv")
taxa_ref <- readr::read_csv("data/taxa_reference.csv")
mydata = fusion_data

#mydata <- readr::read_tsv("data/fusion_data_working2.tsv")
mydata2 <- dtplyr::lazy_dt(mydata)

message(Sys.time(), " Removing Singletons")
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
phylum_s = "Proteobacteria"
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


message(Sys.time(), " Dimensions of data:", print(dim(d1)))
 # three million x one thousand

message(Sys.time(), " Finding interval")
d1[ , cast_cat := findInterval(fusion_lvl_1, seq(10000 , 200000 , 10000) ) ]
message(Sys.time(), " Splitting by interval")
w_list <- split(d1 , by = 'cast_cat' )
w_list <- mclapply( w_list , function( x ) x[ , cast_cat := NULL ] , mc.cores = N)
message(Sys.time(), " casting each element into the list as wide")
w_list <- mclapply( w_list , function( z ) data.table::dcast( z , assembly_accession ~ fusion_lvl_1 , value.var = 'val' , drop = FALSE, fill = 0)  , mc.cores = N)
message(Sys.time(), " attempting to join the list")
result <- Reduce( function( ... ) merge( ... , by = 'assembly_accession' , all = TRUE ) , w_list )

message(Sys.time(), " casting into wide format")
# cast into wide format.
binDTM <- result %>%
  column_to_rownames("assembly_accession") %>%
  as.matrix() %>%
  t() %>%
  as.data.frame() %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  as.matrix()

message(Sys.time(), " calculating distance matrix")
dist = parDist(binDTM, method = "manhattan", threads = N)  

message(Sys.time(), " calculating distance matrix")
clust = fastcluster::hclust(dist, method = "complete")


