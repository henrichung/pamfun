library(fastcluster)
library(parallelDist)

binDTM <- readRDS("subsample_binDTM.RDS")

message(Sys.time(), " | DISTANCE CALCULATION START")
dist <- parDist(binDTM)
message(Sys.time(), " | DISTANCE CALCULATION END")
saveRDS(dist, "binDTM_dist.RDS")

message(Sys.time(), " | CLUSTERING START")
clust <- fastcluster::hclust(dist, method = "complete")
message(Sys.time(), " | CLUSTERING END")
saveRDS(clust, "cluster.RDS")
