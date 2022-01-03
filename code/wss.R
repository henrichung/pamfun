library(parallel)
library(tidyverse)
library(vegan)
#####

message(Sys.time(), " READING DATA")
binDTM <- readRDS("binDTM_nosingletons_cyano.RDS") 

# convert binDTM to appropriate format
cluster_by_sample <- binDTM %>%
	replace(is.na(.), 0) 

# calculate distance matrix using jaccard
message(Sys.time(), " STARTING DISTANCE CLUSTERING")
dist_mat <- vegan::vegdist(cluster_by_sample, method = "jaccard")
message(Sys.time(), " FINISHED DISTANCE CLUSTERING")


#clust <- hclust(dist_mat, method = "complete")
message(Sys.time(), " STARTING FAST CLUSTERING")
clust_agnes <- fastcluster::hclust(dist_mat, method = "complete")
message(Sys.time(), " FINISHED FAST CLUSTERING")

# determine optimal number of clusters
calc_ss <- function(df) sum(as.matrix(vegan::vegdist(df)^2)) / (2 * nrow(df))

custom_wss <- function(k){
  sub_grp <- cutree(clust_agnes, k = k)

  temp <- cluster_by_sample %>%
    as.data.frame() %>%
    rownames_to_column("fusion_lvl_1") %>%
    mutate(cluster = sub_grp) %>%
    select(fusion_lvl_1, cluster, everything()) %>%
    group_by(cluster) %>%
    nest() %>%
    mutate(data = map(data, function(.x){column_to_rownames(.x, "fusion_lvl_1")})) %>%
    mutate(wss = map(data, function(.x){calc_ss(.x)})) %>%
    unnest(wss)

  return(temp$wss)
}
k <- as.list(c(5, 10, 25, 50, 100, 150, 200, 300, 500, 1000, 2000))
message(Sys.time(), " STARTING WSS CALCULATION")
elbow_data <- mclapply(k, function(.x){custom_wss(.x)}, mc.cores = 6)
message(Sys.time(), " FINISHED WSS CALCULATION")

saveRDS(elbow_data, "elbow_cyano.RDS")
quit()