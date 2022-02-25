## ---------------------------
## Purpose of script: Scrape EC pathways and protist information from KEGG website
## Author: Henri Chung
## Date Created: 2021-01-12
## ---------------------------
## Notes:
## ---------------------------

# Load necessary libraries
library(fastcluster)
library(parallelDist)
library(stringr)
rm(list = ls())
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = args[1]
}

file_in <- args[1]
file_out <- str_remove(args[2], ".RDS")
# read in subsample binary matrix
binDTM <- as.matrix(readRDS(file_in))
#outputs/4_reduced_binDTM1.RDS
# create random matrix by permuting columns of sample binDTM
rand_binDTM <- t(apply(binDTM, 1, function(.x){ sample(size = length(.x), x = .x, replace = FALSE)}))

# Calculate distance and cluster random distance matrix.
message(Sys.time(), " | RAND DISTANCE CALCULATION START")
rand_dist <- parDist(rand_binDTM, method = "binary")
message(Sys.time(), " | RAND DISTANCE CALCULATION END")
saveRDS(rand_dist, paste(file_out, "_rand_dist.RDS", sep = ""))

message(Sys.time(), " | RAND CLUSTERING START")
rand_clust <- fastcluster::hclust(rand_dist, method = "complete")
message(Sys.time(), " | RAND CLUSTERING END")
saveRDS(rand_clust, paste(file_out, "_rand_cluster.RDS", sep = ""))

message(Sys.time(), " | CONVERT DISTANCE TO MATRIX")
rand_mat <- as.matrix(rand_dist)
saveRDS(rand_mat, paste(file_out, "_rand_distmat.RDS", sep = ""))


# Calculate distance and cluster sample matrix.
message(Sys.time(), " | DISTANCE CALCULATION START")
dist <- parDist(binDTM, method = "binary")
message(Sys.time(), " | DISTANCE CALCULATION END")
saveRDS(dist, paste(file_out, "_dist.RDS", sep = ""))

message(Sys.time(), " | CLUSTERING START")
clust <- fastcluster::hclust(dist, method = "complete")
message(Sys.time(), " | CLUSTERING END")
saveRDS(clust, paste(file_out, "_cluster.RDS", sep = ""))

message(Sys.time(), " | CONVERT DISTANCE TO MATRIX")
dist_mat <- as.matrix(dist)
saveRDS(dist_mat, paste(file_out, "_distmat.RDS", sep = ""))
