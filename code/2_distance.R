## ---------------------------
## Purpose of script: Rude potential fusion orthologs.
## Author: Henri Chung
## Date Created: 2021-01-12
## ---------------------------
## Notes:
## ---------------------------

## Usage

#Rscript 2_distance.R binDTM.RDS -n 1 -f data/fusion_to_ec.csv -o outputs/red_binDTM.RDS


# Package names
packages <- c("fastcluster", "parallelDist", "stringr", "tidyverse", "dtplyr", "data.table", "argparse")
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
parser$add_argument('matrix', type= "character", help='Binary fusion matrix.')
parser$add_argument("-n", "--threads", type="integer", default= 1, help="Number of threads to use in parallel distance calculation.")
parser$add_argument("-f", "--ecmap", type="character", default= "fusion_to_ec.csv", help="Reference table mapping fusion #s to EC values.")
parser$add_argument("-o", "--output", type="character", default= "reduced_binDTM.RDS", help="Output file name.")

args <- parser$parse_args()

N <- args$threads
file_in <- args$matrix
file_out <- args$output
# read in subsample binary matrix
binDTM <- as.matrix(readRDS(file_in))

fusion_to_ec <- read_csv(args$ecmap) %>%
	mutate(fusion_lvl_1 = as.character(fusion_lvl_1))

# given a pair of rows and matrix, return the 
# binary sum of those rows.
join_rows <- function(pair, matrix){
	pair <- as.character(pair)
	sum_row <- colSums(matrix[pair,])
	sum_row[sum_row > 1] <- 1
	return(sum_row)
}

reduce_dtm <- function(row_pairs, binDTM){
	row_pairs2 <- split(as.data.table(row_pairs), seq(nrow(as.data.table(row_pairs))))
	# get names of rows to be summed
	ind <- unlist(unique(row_pairs2))
	# separate rows that do not need to be summed
	binDTM_base <- binDTM[!(rownames(binDTM) %in% ind),]

	# go over each pair of rows and sum them, bind the new rows into a new matrix
	binDTM_joined <- as.matrix(bind_rows(lapply(row_pairs2, function(x){join_rows(pair = x, matrix = binDTM)})))
	# apply new rownames to the new summed rows
	rownames(binDTM_joined) <- unlist(lapply(row_pairs2, function(x){paste(as.character(x), collapse = "_")}))

	# return the new summed rows with the base rows
	res <- rbind(binDTM_base, binDTM_joined)
	return(res)
}

# given a key table of fusion clusters to EC numbers.
# function to calculate all two-way combinations of fustion clusters
# that point to the same EC
key_pairs <- function(key){
	self <- key %>%
		rename(fusion_lvl_1a = "fusion_lvl_1") %>%
		mutate(fusion_lvl_1b = fusion_lvl_1a) %>%
		group_by(ec) %>%
		complete(fusion_lvl_1a, fusion_lvl_1b) %>%
		filter(fusion_lvl_1a != fusion_lvl_1b) %>%
		rowwise() %>%
		mutate(fusion_lvl_1 = paste(fusion_lvl_1a, fusion_lvl_1b)) %>%
		select(ec, fusion_lvl_1) %>%
		mutate(fusion_lvl_1 = gsub(" ", "_", fusion_lvl_1)) %>%

	return(self)
}

collapse_orthologs <- function(binDTM, fusion_ec_key, dist_value = 0.80, threads = 1){

	# Calculate binary distance for DTM
	dist <- parDist(binDTM, method = "binary")

	# convert to matrix and reshape to long format
	dist_mat<- as.matrix(dist)
	dist_mat[lower.tri(dist_mat, diag = TRUE)] <- NA
	dist_molten <- reshape2::melt(dist_mat)
	
	# calculate fusion clusters which point to the same EC annotation
	fusion_self <- key_pairs(fusion_ec_key)

	# calculate pairs of fusion clusters that have 
	# large distances and point to the same EC enzyme number
	dist_dt <- dtplyr::lazy_dt(dist_molten)

	selfs <- dist_dt %>%
		filter(!is.na(value)) %>%
		filter(value > dist_value) %>%
		mutate(combo = paste(Var1, Var2, sep = "_")) %>%
		filter(combo %in% fusion_self$fusion_lvl_1) %>%
		group_by(Var1) %>% slice_head(n = 1) %>% ungroup() %>%
		as.data.table()
	if(nrow(selfs) == 0 ){message("No rows left to join"); return(list(binDTM, fusion_ec_key, dist_mat))}
	selfs2 <- selfs %>%
		filter(!(Var2 %in% selfs$Var1)) %>%
		group_by(Var2) %>% slice_head(n = 1) %>% ungroup() %>%
		mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
		as.data.table()

	# reduce binDTM by summing over othologous clusters
	new_binDTM <- reduce_dtm(select(selfs2, Var1, Var2), binDTM)

	# return new fusion key with joined clusters
	new_fusion_ec_key <- bind_rows(filter(fusion_self, fusion_lvl_1 %in% selfs2$combo), fusion_ec_key)
	return(list(new_binDTM, new_fusion_ec_key, dist_mat))
}


# collapse orthologs
fusion3 <- select(rename(fusion_to_ec, ec = "ec_number3"), c("fusion_lvl_1", "ec"))
temp <- collapse_orthologs(binDTM, fusion_ec_key = fusion3, threads = N)
i = 1

while(i > 0){
	res <- collapse_orthologs(temp[[1]], fusion_ec_key = temp[[2]], threads = N)
	if(!identical(dim(temp[[1]]), dim(res[[1]]))){message(Sys.time(), "NROW: ", nrow(temp[[1]]), " ITER: ", i); temp <- res}else{break}
	i = i + 1
}

# save results
names(res) <- c("binDTM", "fusion_ec_key", "dist_mat")
saveRDS(res, file_out)


# create random matrix by permuting columns of sample binDTM
rand_binDTM <- t(apply(res[[1]], 1, function(.x){ sample(size = length(.x), x = .x, replace = FALSE)}))

# Calculate distance distance matrix.
message(Sys.time(), " | Starting random distance calculation.")
rand_dist <- parDist(rand_binDTM, method = "binary", threads = N)
message(Sys.time(), " | Converting random distance to matrix.")
rand_mat <- as.matrix(rand_dist)
message(Sys.time(), " | Saving random distance matrix.")
saveRDS(rand_mat, paste("rand_distmat_",file_out,  sep = ""))
quit()

