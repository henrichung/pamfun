## ---------------------------
## Purpose of script: Scrape EC pathways and protist information from KEGG website
## Author: Henri Chung
## Date Created: 2021-01-12
## ---------------------------
## Notes:
## ---------------------------


## Usage
#Rscript code/4_analze.R red_binDTM.RDS -n 1 -e data/enzyme_interactions.csv 

#Rscript 2_distnace.R binDTM.RDS -n 1 -f data/fusion_to_ec.csv -o outputs/red_binDTM.RDS
packages <- c("Matrix", "parallel", "stringr", "tidyverse", "dtplyr", "data.table", "argparse")
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
parser$add_argument('dist_mat', type= "character", help='Binary fusion matrix.')
parser$add_argument("-o", "--output", type="character", default= "reduced_binDTM.RDS", help="Output file name.")
parser$add_argument("-n", "--threads", type="integer", default= 1, help="Number of threads to use in parallel distance calculation.")
parser$add_argument("-e", "--enzymes", type="character", default= "enzyme_interactions.csv", help="Reference table mapping fusion #s to EC values.")
parser$add_argument("-f", "--ecmap", type="character", default= "fusion_to_ec.csv", help="Reference table mapping fusion #s to EC values.")

args <- parser$parse_args()

dist_mat <- readRDS(args$dist_mat)
fusion_to_ec <- read_csv(args$ecmap)
enzyme_interactions <- read_csv(args$enzymes)

fusion_to_ec <- read_csv("data/fusion_to_ec.csv")
res <- readRDS("reduced_binDTM.RDS")
rand_dist_mat <- readRDS("rand_distmat_reduced_binDTM.RDS")

paths <- readRDS("data/protist_pathways.RDS") %>%
	bind_rows() %>%
	pull(num) %>%
	unique() 
bacteria_ecs <- paste("ec", paths, sep = "")

rand_dist_mat[lower.tri(rand_dist_mat, diag = TRUE)] <- NA
rand_distances <- rand_dist_mat %>%
	reshape2::melt() %>% 
	filter(!is.na(value)) %>%
	mutate(combo = NA, id = NA)
# dataframe of two way enzyme interactions in a pathway
modules_df <- read_csv("data/module_info.csv")
enzyme_interactions <- read_csv("data/enzyme_interactions.csv") 

ec_enzymes <- enzyme_interactions %>%
	filter(group == "pathway") %>%
	rename(id = "ec_id")
module_enzymes <- enzyme_interactions %>%
	filter(group == "module") %>%
	rename(id = "module_id") %>%
	filter(id %in% modules_df$module_id)

# dataframe of two way enzyme interactions in a pathway
ec_enzymes3 <- ec_enzymes %>%
	mutate(enzymes_a = gsub("\\.[0-9]+\\Z", ".-", enzymes_a, perl = TRUE)) %>%
	mutate(enzymes_b = gsub("\\.[0-9]+\\Z", ".-", enzymes_b, perl = TRUE)) %>%
	unique()
module_enzymes3 <- module_enzymes %>%
	mutate(enzymes_a = gsub("\\.[0-9]+\\Z", ".-", enzymes_a, perl = TRUE)) %>%
	mutate(enzymes_b = gsub("\\.[0-9]+\\Z", ".-", enzymes_b, perl = TRUE)) %>%
	filter(enzymes_a != enzymes_b) %>%
	unique()

# output data from distance calculations


#1.2.7.-
#2.7.1.-

custom_analyze <- function(data,enzymes){
	# extract data from output
	binDTM <- data[[1]]
	fusion_ec_key <- data[[2]] 
	dist_mat <-  data[[3]]

	# reshape to long format
	mod_fusion <- fusion_key(fusion_ec_key = fusion_ec_key, ec_enzymes =  enzymes)

	dist_molten <- reshape2::melt(dist_mat) %>%
		filter(!is.na(value)) %>%
		mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
		rowwise() %>%
		mutate(combo = paste(Var1, Var2)) %>%
		left_join(select(mod_fusion, id, combo)) %>% 	
		filter(Var1 != Var2) %>%
		distinct(combo, .keep_all = TRUE) %>%
		as.data.table()

	pathways_dist <- dist_molten %>% filter(combo %in% mod_fusion$combo) %>% as.data.frame()

	other_dist <- dist_molten %>% filter(!(combo %in% mod_fusion$combo)) %>% as.data.frame()

	res <- list(pathways_dist, other_dist)
	names(res) <- c("shared", "unshared")
	return(res)

}

fusion_key <- function(fusion_ec_key, ec_enzymes){
	ec_dt <- dtplyr::lazy_dt(ec_enzymes)
	res2 <- ec_dt %>%
		rename(ec = enzymes_a) %>% left_join(fusion_ec_key, by = "ec") %>% rename(fusion_a = fusion_lvl_1) %>% 
		rename(enzyme_a = ec) %>% 
		rename(ec = enzymes_b) %>% left_join(fusion_ec_key,  by = "ec") %>% rename(fusion_b = fusion_lvl_1) %>% 
		rename(enzyme_b = ec) %>% 
		mutate(combo = paste(fusion_a, fusion_b)) %>%
		filter(fusion_a != fusion_b)
	return(as.data.frame(res2))
}

ec3_module <- custom_analyze(res,  enzymes = module_enzymes3)
distances <- append(ec3_module,list(rand_distances))
names(distances) <- c(names(ec3_module), "random")


p_data <- distances %>%
	bind_rows(.id = "group") %>%
	rowwise() %>%
	mutate(value = as.numeric(value)) %>%
	ungroup() %>%
	mutate(group = factor(group, levels =  c("shared", "unshared", "random")))


p1 <- p_data %>%
	ggplot(aes(x = value, y = stat(density),  fill = group,  linetype = group, color = group)) + 
	geom_density(alpha = 0.1) +
	xlab("Jaccard Distance") + ylab("Density") +
	scale_fill_brewer(palette="Dark2") + scale_color_brewer(palette="Dark2") + 
	labs(title = "Jaccard Density.") +
	theme_bw()
pdf("Distribution of Jaccard Clusters.pdf", width = 8, height = 6) 
p1
dev.off()
#

p2_data <- p_data %>%
	filter(group == "shared") %>%
	rename(module_id = "id") %>%
	left_join(modules_df) 

carbs <- filter(p2_data, subclass == "Carbohydrate metabolism") %>% 
	pull(pathway) %>% 
	table() %>% 
	sort() %>% 
	stack() %>%
	filter(values > 500) %>% 
	pull(ind) %>%
	as.character()


p2 <- p2_data %>%
	filter(subclass == "Carbohydrate metabolism" & pathway %in% carbs) %>%
	ggplot(aes(x = value, y = stat(density),  fill = pathway, color = pathway)) + 
	geom_density(alpha = 0) +
	xlab("Jaccard Distance") + ylab("Density") +
	scale_fill_brewer(palette="Dark2") + scale_color_brewer(palette="Dark2") + 
	labs(title = "Carbohydrate Metabolism") +
	theme_bw()

pdf("Distribution of Jaccard Clusters - Carbohydrate Metabolism.pdf", width = 8, height = 6) 
p2
dev.off()

#
a <- filter(p2_data, grepl("Gluconeogenesis", desc)) %>% pull(combo) %>% unique()
a <- unlist(lapply(a, str_split, pattern = " "))

b <- filter(p2_data, grepl("Glycolysis", desc)) %>% pull(combo) %>% unique()
b <- unlist(lapply(b, str_split, pattern = " "))

c <- intersect(a, b)

p3_data <- p2_data %>%
	filter(pathway == "Glycolysis / Gluconeogenesis") %>%
	filter(!(Var1 %in% c) & !(Var2 %in% c))


p3 <- p3_data %>%
	ggplot(aes(x = value, y = stat(density),  fill = desc, color = desc)) + 
	geom_density(alpha = 0) +
	xlab("Jaccard Distance") + ylab("Density") +
	scale_fill_brewer(palette="Dark2", name = "module") + scale_color_brewer(palette="Dark2") + 
	labs(title = "Glycolysis/Gluconeogenesis.") +
	theme_bw()
pdf("Distribution of Jaccard Clusters - Glycolysis_Gluconeogenesis (no overlap).pdf", width = 8, height = 6) 
p3
dev.off()	


p5 <- p2_data %>%
	filter(pathway == "Glycolysis / Gluconeogenesis") %>%
	ggplot(aes(x = value, y = stat(density),  fill = desc, color = desc)) + 
	geom_density(alpha = 0) +
	xlab("Jaccard Distance") + ylab("Density") +
	scale_fill_brewer(palette="Dark2") + scale_color_brewer(palette="Dark2") + 
	labs(title = "Glycolysis/Gluconeogenesis.") +
	theme_bw()
pdf("Distribution of Jaccard Clusters - Glycolysis_Gluconeogenesis (overlap).pdf", width = 8, height = 6) 
p5
dev.off()	

p4 <- p3_data %>%
	ggplot(aes(x = desc, y =value,  fill = desc, color = desc)) + 
	geom_violin() +
	xlab("") + ylab("Density") +
	scale_fill_brewer(palette="Dark2") + scale_color_brewer(palette="Dark2") + 
	labs(title = "Glycolysis/Gluconeogenesis.") +
	theme_bw()
pdf("Distribution of Jaccard Clusters - Glycolysis_Gluconeogenesis_violin.pdf", width = 8, height = 6) 
p4
dev.off()	
#
#quit()

write_csv(p3_data, "outputs/glyco_gluco.csv")