## ---------------------------
## Purpose of script: Use the results of first distance calculation to reduce potential orthologues.
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

# recalculate distances

fusion_to_ec <- read_csv("outputs/0_fusion_to_ec.csv") %>%
	select(ec_number, fusion_lvl_1) %>%
	filter(!is.na(fusion_lvl_1)) %>%
	unique() 

dist_mat<- readRDS("outputs/1_subsample_binDTM_distmat.RDS")

ec_molten <- reshape2::melt(dist_mat) %>%
		filter(!is.na(value)) %>%
		mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
		rowwise() %>%
		mutate(combo = paste(Var1, Var2)) 

pre_df <- bind_rows(ec_list, .id = "ec_id") %>%
	rename(ec_number = Var1) %>% left_join(fusion_to_ec) %>% rename(Var1 = fusion_lvl_1) %>% rename(ec1 = ec_number) %>%
	rename(ec_number = Var2) %>% left_join(fusion_to_ec) %>% rename(Var2 = fusion_lvl_1) %>% rename(ec2 = ec_number) %>%
	rowwise() %>%
	filter(Var1 != "NA" & Var2 != "NA" & Var1 != Var2) %>%
	filter(Var1 %in% rownames(dist_mat)) %>% filter(Var2 %in% rownames(dist_mat)) %>%
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


pathways_dist <- filter(ec_molten_dt, combo %in% ec_df$combo) %>%
	left_join(select(ec_df, ec_id, combo)) %>% 	
	filter(Var1 != Var2) %>%
	distinct(combo, .keep_all = TRUE) %>%
	filter(ec_id %in% protist_ec) %>%
	as.data.table()
#saveRDS(pathways_dist, "pathways_dist.RDS")
other_dist <- filter(ec_molten_dt, !(combo %in% ec_df$combo)) %>% 
	left_join(select(ec_df, ec_id, combo)) %>%
	filter(Var1 != Var2) %>%
	distinct(combo, .keep_all = TRUE) %>%
	filter(!(ec_id %in% protist_ec)) %>%
	as.data.table()  

pathways_dist2 <- filter(ec_molten_dt2, combo %in% ec_df2$combo) %>%
	left_join(select(ec_df2, ec_id, combo)) %>% 	
	filter(Var1 != Var2) %>%
	distinct(combo, .keep_all = TRUE) %>%
	filter(ec_id %in% protist_ec) %>%
	as.data.table()

#saveRDS(pathways_dist, "pathways_dist.RDS")
other_dist2 <- filter(ec_molten_dt2, !(combo %in% ec_df2$combo)) %>% 
	left_join(select(ec_df2, ec_id, combo)) %>%
	filter(Var1 != Var2) %>%
	distinct(combo, .keep_all = TRUE) %>%
	filter(!(ec_id %in% protist_ec)) %>%
	as.data.table()  
	
dim(pathways_dist2)
dim(other_dist2)

summary(pathways_dist2$value)
summary(other_dist2$value)

nrow(pathways_dist2)
nrow(other_dist2)

test <- pathways_dist2 %>%
	left_join(protist_df)

# quick ec lookup of pathways
library(rvest)
ecs <- unique(pathways_dist2$ec_id)
base_url <- "https://www.genome.jp/entry/"
desc_list <- list()
for(i in 1:length(ecs)){
	a <- read_html(paste(base_url,ecs[i], sep = "")) %>%
		html_nodes("tr") %>%
		html_text()
	b <- a[grepl("Class", a)]
	desc[[i]] <- b[nchar(b) < 100]
	message(i)
	}
names(desc_list) <- ecs

desc <- desc_list %>%
	stack() %>%
	as.data.frame() %>%
	rename(ec_id = "ind", desc = "values") %>%
	rowwise() %>%
	mutate(desc = gsub("BRITE hierarchy", "", desc)) %>%
	mutate(desc = gsub("Class\n", "", desc)) %>%
	separate(desc, into = c("Class", "Process"))
saveRDS(desc, "desc.RDS")

# Comparison to Random Distance Matrix 
rand_dist_mat <- readRDS("outputs/4_reduced_binDTM1_rand_distmat.RDS")
rand_cols <- colnames(rand_dist_mat) %in% fusion_to_ec2$fusion_lvl_1
rand_rows <- rownames(rand_dist_mat) %in% fusion_to_ec2$fusion_lvl_1
sub_dist_mat <- rand_dist_mat[rand_cols, rand_rows]
rand_molten <- reshape2::melt(sub_dist_mat) %>%
	mutate(group = "random", ec_id = NA, combo = "NA", Var1 = 1:nrow(.), Var2 = 1:nrow(.), value = value)


#colnames(rand_molten) <- c("x", "y", "values")
#rand_molten <- filter(rand_molten, values != 0)
#saveRDS(rand_molten, "rand_molten.RDS")

p1_list <- list(as.data.frame(other_dist2), as.data.frame(pathways_dist2),as.data.frame(other_dist), as.data.frame(pathways_dist))
names(p1_list) = c("unshared_post", "shared_post", "unshared_pre", "shared_pre") 
p1_data <- bind_rows(p1_list, .id = "group") %>%
	rowwise() %>%
	mutate(value = as.numeric(value)) %>%
	ungroup() %>%
	rbind(rand_molten) %>%
	mutate(group = factor(group, levels = c("random", "unshared_pre", "unshared_post", "shared_pre", "shared_post")))

p1 <- p1_data %>%
	ggplot(aes(x = group, y = value, fill = group)) + 
	geom_violin() +
	ylab("Jaccard Distance") + xlab("") + 
	scale_fill_brewer(palette="Pastel1") + 
	labs(title = "Distribution of Vector Distance between functions.", caption = "N Shared = 18,147, N Unshared = 2,598,132")
pdf("Distribution of Vector Distance between functions OP.pdf", width = 8, height = 6)
p1
dev.off()

p2 <- p1_data %>%
  ggplot(aes(x = value, color = group)) +
  geom_density(alpha = .2)  +
  xlab("Jaccard Distance") + ylab("Density") +
  scale_fill_brewer(palette="Pastel1") +
labs(title = "Density of Vector Distance between functions.", caption = "N Shared = 18,147, N Unshared = 2,598,132")
pdf("Density of Vector Distance between functions OP.pdf") 
p2
dev.off()

# Before

quit()
