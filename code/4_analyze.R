## ---------------------------
## Purpose of script: Scrape EC pathways and protist information from KEGG website
## Author: Henri Chung
## Date Created: 2021-01-12
## ---------------------------
## Notes:
## ---------------------------

# Load necessary libraries
library(tidyverse)
library(dtplyr)
library(data.table)
library(parallel)
library(Matrix)
rm(list = ls())

fusion_to_ec2 <- readRDS("outputs/4_fusion_to_ec.RDS")
ec_list <- readRDS("outputs/4_ec_list.RDS")
dist_mat2 <- readRDS("outputs/4_reduced_binDTM1_distmat.RDS")

subset_melt <- function(mat, clusters){
	cols <- colnames(mat) %in% clusters
	rows <- rownames(mat) %in% clusters
	sub_dist <- mat[cols, rows]
	#sub_dist[lower.tri(sub_dist, diag = TRUE)] <- NA
	# reshape distance matrix to long format
	dist_molten <- reshape2::melt(sub_dist) %>%
		filter(!is.na(value)) %>%
		mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
		rowwise() %>%
		mutate(combo = paste(Var1, Var2)) 
}

ec_molten2 <- subset_melt(dist_mat2, fusion_to_ec2$fusion_lvl_1)

pre_df2 <- bind_rows(ec_list, .id = "ec_id") %>%
	rename(ec_number = Var1) %>% left_join(fusion_to_ec2) %>% rename(Var1 = fusion_lvl_1) %>% rename(ec1 = ec_number) %>%
	rename(ec_number = Var2) %>% left_join(fusion_to_ec2) %>% rename(Var2 = fusion_lvl_1) %>% rename(ec2 = ec_number) %>%
	rowwise() %>%
	filter(Var1 != "NA" & Var2 != "NA" & Var1 != Var2) %>%
	filter(Var1 %in% rownames(dist_mat2)) %>% filter(Var2 %in% rownames(dist_mat2)) %>%
	mutate(combo = paste(Var1, Var2)) %>%
	select(ec_id, ec1, ec2, combo) %>%
	unique()



self2 <- filter(pre_df2, ec1 == ec2)
ec_df2 <- filter(pre_df2, ec1 != ec2)

#self2 <- filter(self2, ec_id %in% most_common)
#ec_df2 <- filter(ec_df2, ec_id %in% most_common)

ec_molten_dt2 <- dtplyr::lazy_dt(ec_molten2)

self_combine2 <- ec_molten2 %>%
	filter(value > 0.80) %>%
	filter(combo %in% self2$combo) %>%
	group_by(Var1) %>%
	slice_head(n = 1) %>% 
	ungroup() %>%
	group_by(Var2) %>%
	slice_head(n = 1) %>%
	ungroup() %>%
	arrange(desc(value)) %>%
	as.data.table() 

# filter to just protist pathways 
protist_pathways <- readRDS("outputs/2_protist_pathways.RDS")

protist_df <- protist_pathways %>%
	bind_rows(.id = "abbr") %>%
	filter(num != "") %>%
	select(-abbr) %>% 
	unique() %>%
	rowwise() %>%
	mutate(num = paste("ec", num, sep = "")) %>%
	rename(ec_id = "num")

protist_ec <- protist_df %>%
	pull(ec_id)

#####

pathways_dist <- filter(ec_molten_dt2, combo %in% ec_df2$combo) %>%
	left_join(select(ec_df2, ec_id, combo)) %>% 	
	filter(Var1 != Var2) %>%
	distinct(combo, .keep_all = TRUE) %>%
	filter(ec_id %in% protist_ec) %>%
	as.data.table()

#test <- as.data.frame(filter(pathways_dist, ec_id == "ec00970"))
#test2 <- filter(ec_df2, ec_id == "ec00970")
dim(as.data.table(test))
#saveRDS(pathways_dist, "pathways_dist.RDS")
other_dist <- filter(ec_molten_dt2, !(combo %in% ec_df2$combo)) %>% 
	left_join(select(ec_df2, ec_id, combo)) %>%
	filter(Var1 != Var2) %>%
	distinct(combo, .keep_all = TRUE) %>%
	filter(!(ec_id %in% protist_ec)) %>%
	as.data.table()  

# Comparison to Random Distance Matrix 
rand_dist_mat <- readRDS("outputs/4_reduced_binDTM1_rand_distmat.RDS")
rand_cols <- colnames(rand_dist_mat) %in% fusion_to_ec2$fusion_lvl_1
rand_rows <- rownames(rand_dist_mat) %in% fusion_to_ec2$fusion_lvl_1
sub_dist_mat <- rand_dist_mat[rand_cols, rand_rows]
rand_dist <- reshape2::melt(sub_dist_mat) %>%
	mutate(ec_id = NA, combo = "NA", Var1 = as.character(1:nrow(.)), Var2 = as.character(1:nrow(.)), value = value) %>%
	sample_n(nrow(other_dist)) %>%
	as.data.table()

#colnames(rand_molten) <- c("x", "y", "values")
#rand_molten <- filter(rand_molten, values != 0)
#saveRDS(rand_molten, "rand_molten.RDS")

p1_list <- list(as.data.frame(other_dist), as.data.frame(pathways_dist), as.data.frame(rand_dist))
names(p1_list) = c("unshared", "shared", "random") 
p1_data <- bind_rows(p1_list, .id = "group") %>%
	rowwise() %>%
	mutate(value = as.numeric(value)) %>%
	ungroup() %>%
	mutate(group = factor(group, levels = c("unshared", "shared", "random")))

p1 <- p1_data %>%
	ggplot(aes(x = group, y = value, fill = group)) + 
	geom_violin() +
	ylab("Jaccard Distance") + xlab("") + 
	scale_fill_brewer(palette="Pastel1") + 
	labs(title = "Distribution of Vector Distance between functions.", caption = "N Shared = 4,657, N Unshared = 1,080,941, N Random = 1080941")
pdf("Distribution of Vector Distance between functions.pdf", width = 8, height = 6)
p1
dev.off()


ec_pathways <- readRDS("outputs/EC_pathways.RDS")
desc <- read_csv("outputs/EC_pathway_desc.csv")


rand_dist2 <- sample_n(rand_dist, nrow(other_dist)) %>% as.data.frame()
p2_data <- p1_data %>%
	right_join(desc, by = "ec_id") %>% 
	as.data.frame() 

enzyme_counts <- ec_pathways %>%
	filter(ind %in% unique(pathways_dist$ec_id)) %>%
	filter(values %in% unique(fusion_to_ec2$ec_number)) %>%
	group_by(ind) %>%
	summarize(enzymes_n = n()) %>%
	mutate(ind = as.character(ind)) %>%
	rename(ec_id = "ind")

p3_data <- p2_data %>%
	group_by(Subprocess) %>%
	summarize(n = n()) %>%
	left_join(p2_data) %>%
	left_join(enzyme_counts) %>%
	mutate(Subprocess = paste(Subprocess, " (", n, ")[", enzymes_n, "]", sep = "")) %>%
	plyr::rbind.fill(rand_dist2) %>%
	mutate(Subprocess = ifelse(is.na(Class), "Random", Subprocess), Process = ifelse(is.na(Class), "Random", Process)) %>%
	mutate(metabolism = ifelse(grepl("biosynthesis", Subprocess), "Biosynthesis", "Metabolism")) %>%
	filter(!is.na(value)) 

filter(p3_data, Subprocess == "Random") %>% head()
unique(p3_data$Process)

custom_plot <- function(x, y){
	res <- x %>%
		filter(Process == y) %>%
		ggplot(aes(x = Subprocess, y = value)) + 
		geom_violin(fill = "gray") +
		ylab("Jaccard Distance") + xlab("") + 
		coord_flip() +
		labs(title = y) +
		theme_bw() +
		facet_wrap(~metabolism, nrow = 2, scales = "free") +
		labs(caption = "(#) refers to the number of distance values in each subprocess. \n [#] refers to the number of possible fusion associated enzymes involved in each subprocess.")
  	return(res) 
}

p3a <- custom_plot(p3_data, "Carbohydrate metabolism")
p3b <- custom_plot(p3_data, "Amino acid metabolism")
p3c <- custom_plot(p3_data, "Metabolism of cofactors and vitamins")
p3d <- custom_plot(p3_data, "Glycan biosynthesis and metabolism")
p3e <- custom_plot(p3_data, "Metabolism of other amino acids")
p3f <- custom_plot(p3_data, "Nucleotide metabolism")
p3g <- custom_plot(p3_data, "Energy metabolism")
p3h <- custom_plot(p3_data, "Translation")
p3i <- custom_plot(p3_data, "Xenobiotics biodegradation and metabolism")
p3j <- custom_plot(p3_data, "Biosynthesis of other secondary metabolites")
p3k <- custom_plot(p3_data, "Metabolism of terpenoids and polyketides")
p3l <- custom_plot(p3_data, "Random")

pdf("Distributions by Subprocess.pdf", width = 10, height = 6) 
p3a
p3b
p3c
p3d
p3e
p3f
p3g
p3h
p3i
p3j
p3k
p3l
dev.off()

library(igraph)
graph_data <- p3b$data %>%
	filter(Subprocess == "Phenylalanine, tyrosine and tryptophan biosynthesis (63)[12]")

f1_data <- select(graph_data, c("Var1", "Var2", "value")) %>% mutate(value = as.numeric(value))
f1 <-  graph_from_data_frame(f1_data, directed = FALSE)
E(f1)$weight = f1_data$value
E(f1)$width <- (E(f1)$weight + min(E(f1)$weight)) * 10 # offset=1
pdf("test.pdf")
#plot(f1, edge.label=round(E(f1)$weight, 3))
plot(f1)
dev.off()
# Before
a <- filter(p3_data, Process == "Nucleotide metabolism")
c(unique(a$Var1, a$Var2)) %>% length()


b <- filter(p3_data, Subprocess == "O-Antigen nucleotide sugar biosynthesis (161)[63]")
c(unique(b$Var1, b$Var2)) %>% length()


quit()
