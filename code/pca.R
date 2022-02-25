library(fastcluster)
library(parallelDist)

rand_binDTM <- readRDS("subsample_binDTM.RDS")

message(Sys.time(), " | RAND DISTANCE CALCULATION START")
rand_dist <- parDist(rand_binDTM, method = "binary")
message(Sys.time(), " | RAND DISTANCE CALCULATION END")
saveRDS(rand_dist, "rand_binDTM_dist.RDS")

message(Sys.time(), " | RAND CLUSTERING START")
rand_clust <- fastcluster::hclust(rand_dist, method = "complete")
message(Sys.time(), " | RAND CLUSTERING END")
saveRDS(rand_clust, "rand_cluster.RDS")


binDTM <- readRDS("subsample_binDTM.RDS")

message(Sys.time(), " | DISTANCE CALCULATION START")
dist <- parDist(binDTM, method = "binary")
message(Sys.time(), " | DISTANCE CALCULATION END")
saveRDS(dist, "binDTM_dist.RDS")

message(Sys.time(), " | CLUSTERING START")
clust <- fastcluster::hclust(dist, method = "complete")
message(Sys.time(), " | CLUSTERING END")
saveRDS(clust, "cluster.RDS")
