## ---------------------------
## Purpose of script: Subset the fusion data to relevant clusters and convert into binary matrix
## Author: Henri Chung
## Date Created: 2021-01-12
## ---------------------------
## Notes:
## ---------------------------

# Package names
packages <- c("data.table", "tidyverse", "dtplyr", "Matrix", "parallel", "vegan", "argparse")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
lapply(packages, function(x){suppressPackageStartupMessages(library(x, character.only = TRUE))})

# clear working directory
rm(list = ls())

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument('fusion', type= "character", help='Fusion data to convert to binary.')
parser$add_argument("-n", "--cores", default = 1, help="Number of threads to use in parallel distance calculation.")
parser$add_argument("-k", "--min", default = 5, help="Minimum number of taxa to be considered in binary matrix.")
parser$add_argument("-f", "--ecmap", type="character", default= "fusion_to_ec.csv", help="Reference table mapping fusion #s to EC values.")
parser$add_argument("-o", "--output", type="character", default= "binDTM.RDS", help="Output file name.")

args <- parser$parse_args()

N = args$cores
k = args$min
file_out <- args$output

message(Sys.time(), " | Number of cores: ", N)
message(Sys.time(), " | Reading in data")

# subset_data <- read_tsv("data/fusion_data_subset.tsv")
subset_data <- read_tsv(args$fusion)
# fusion_to_ec <- read_tsv("data/fusion_to_ec.csv")
fusion_to_ec <- read_csv(args$ecmap)

# select assembly and fusion information from subset data.
cast_data <- subset_data %>%
  select(assembly_accession, fusion_lvl_1) %>%
  filter(fusion_lvl_1 %in% fusion_to_ec$fusion_lvl_1) %>%
  group_by(fusion_lvl_1) %>%
  mutate(val = 1) %>%
  as.data.table()

# Split data
message(Sys.time(), " | Finding interval.")
cast_data[ , cast_cat := findInterval(fusion_lvl_1, seq(1, nrow(cast_data), 10000))]
message(Sys.time(), " | Splitting by interval.")
w_list <- split(cast_data , by = 'cast_cat' )
w_list <- mclapply( w_list , function( x ) x[ , cast_cat := NULL ] , mc.cores = N)
message(Sys.time(), " | Casting each element into the list as wide.")
w_list <- mclapply( w_list , function( z ) data.table::dcast( z , assembly_accession ~ fusion_lvl_1 , value.var = 'val' , drop = FALSE, fill = 0, fun.aggregate = length)  , mc.cores = N)
message(Sys.time(), " | Attempting to join the list.")
result <- Reduce( function( ... ) merge( ... , by = 'assembly_accession' , all = TRUE ) , w_list )

# Cast into wide format
message(Sys.time(), " | Reshaping into Matrix.")
binDTM <- result %>%
  column_to_rownames("assembly_accession") %>%
  as.matrix() %>%
  t() %>%
  as.data.frame() %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  as.matrix()

# Binarize
binDTM[binDTM > 1] <- 1
message(Sys.time(), " | Filtering to samples greater than ", k)
# filter our rows with <= 5 samples.
binDTM <- binDTM[rowSums(binDTM) > k,]

# save output
saveRDS(binDTM, file_out)
