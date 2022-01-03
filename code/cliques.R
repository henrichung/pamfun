library(igraph)
library(data.table)
library(tidyverse)
library(dtplyr)
library(Matrix)
library(quanteda)

message("READING COOCCOUNTS")
coocCounts <- readRDS("bin_cooc_sp.RDS")
binDTM <- read_csv("binDTM_nosingletons_.csv")
colnames(coocCounts) <- binDTM$fusion_lvl_1
rownames(coocCounts) <- binDTM$fusion_lvl_1

message("MELTING")
cooc_long <- mefa4::Melt(coocCounts)
verts <- unique(cooc_long[,1])

message("CREATING IGRAPH")
graph <- graph_from_data_frame(cooc_long, directed=TRUE, vertices= verts)

message("IDENTIFYING CLIQUES")
cliques <- max_cliques(graph)
saveRDS(graph, "graph.RDS")
saveRDS(cliques, "cliques.RDS")