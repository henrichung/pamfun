library(tidyverse)
library(dtplyr)
library(data.table)
library(parallel)
#distance <- readRDS("subsample_dist2.RDS")
#dist_mat <- as.matrix(distance)
#saveRDS(dist_mat, "dist_mat.RDS")
dist_mat <- readRDS("dist_mat.RDS")
seguid_to_fusion <- readRDS("seguid_to_fusion.RDS")
cols <- colnames(dist_mat) %in% seguid_to_fusion$fusion_lvl_1
rows <- rownames(dist_mat) %in% seguid_to_fusion$fusion_lvl_1
ec_dist <- dist_mat[cols, rows]
ec_molten <- reshape2::melt(ec_dist)

ec_pathways <- readRDS("EC_pathways.RDS") %>%
	group_by(ind) %>%
	group_split()

a <- lapply(ec_pathways, function(.x){
		mat = matrix(nrow = nrow(.x), ncol = nrow(.x))
		colnames(mat) = .x$values
		rownames(mat) = .x$values
		mat[upper.tri(mat, diag = FALSE)] = 1
		res <- reshape2::melt(mat) %>%
			filter(!is.na(value))
		return(res)
	})

b <- bind_rows(a) %>%
	rename(ec_number = Var1) %>%
	left_join(seguid_to_fusion) %>%
	rename(Var1 = fusion_lvl_1) %>% select(-c(ec_number, seguid)) %>%
	rename(ec_number = Var2) %>%
	left_join(seguid_to_fusion) %>%
	rename(Var2 = fusion_lvl_1) %>% select(-c(ec_number, seguid)) %>%
	rowwise() %>%
	mutate(combo = paste(Var1, Var2))

c <- ec_molten %>%
	rowwise() %>%
	mutate(combo = paste(Var1, Var2)) 

d <- dtplyr::lazy_dt(c)

e <- filter(d, combo %in% b$combo)

f <- as.data.table(e)


t1 <- stack(table(ec_molten$value)/2)
t2 <- stack(table(f$value))
t3 <- list(as.data.frame(t1), as.data.frame(t2))
names(t3) = c("unshared", "shared") 
t4 <- bind_rows(t3, .id = "group")
pdf("test.pdf")
p1 <- t4 %>%
	mutate(val = as.numeric(values) * as.numeric(ind)) %>%
	ggplot(aes(x = group, y = val)) + 
	geom_violin()
p1 <- ggplot(t4, aes(x = ind, y = values, color = group)) + geom_point(size = 0.5)
p1
dev.off()

pdf("test2.pdf")
p2 <- ggplot(t4, aes(x = ind, y = values, color = group)) + geom_point(size = 0.5)
p2
dev.off()
##############
dist_mat <- readRDS("dist_mat.RDS")
seguid_to_fusion <- readRDS("seguid_to_fusion.RDS")
ec_pathways <- readRDS("EC_pathways.RDS") %>%
	group_by(ind) %>%
	group_split()
cols <- colnames(dist_mat) %in% seguid_to_fusion$fusion_lvl_1
rows <- rownames(dist_mat) %in% seguid_to_fusion$fusion_lvl_1
ec_dist <- dist_mat[cols, rows]
ec_molten <- reshape2::melt(ec_dist)


list_group <- function(x1, x2, y){
	temp <- mclapply(y, FUN = function(.y){in_group(x1, x2, y = .y)}, mc.cores = 36L)
	res <- which(unlist(temp) == TRUE)
	if(length(res) == 0){res <- 0}
	return(res)
}

in_group <- function(x1, x2, y){
	if(x1 %in% y$values & x2 %in% y$values){
		return(TRUE)
	}else{
		return(FALSE)
	}
}
message(Sys.time(), " CALCULATING ROWS")
res <- ec_molten %>%
	mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
	rowwise() %>%
	mutate(group = list_group(Var1, Var2, ec_pathways))
saveRDS(res, "res.RDS")
quit()


for(i in 1:nrow(ec_molten)){

	x1 <- ec_molten$Var1[i]
	x2 <- ec_molten$Var2[i]
	temp <- lapply(ec_pathways, function(.y){in_group(x1, x2, y = .y)})
	res <- which(unlist(temp) == TRUE)
	if(length(res) == 0){res <- 0}
	message(i, " " , res)
	ec_molten$group[i] = res
}

####

rm(distance); gc()
dist_mat1 <- dist_mat[1:25000,]
dist_molten1 <- reshape2::melt(dist_mat1)
rm(dist_mat1);gc()
dist_mat2 <- dist_mat[25000:nrow(dist_mat),]
dist_molten2 <- reshape2::melt(dist_mat2)
rm(dist_mat2);gc()
seguid_to_fusion <- readRDS("seguid_to_fusion.RDS")
dist_ec1 <- dist_molten1 %>%
	filter(Var1 %in% seguid_to_fusion$fusion_lvl_1) %>%
	filter(Var2 %in% seguid_to_fusion$fusion_lvl_1)

dist_ec2 <- dist_molten2 %>%
	filter(Var1 %in% seguid_to_fusion$fusion_lvl_1) %>%
	filter(Var2 %in% seguid_to_fusion$fusion_lvl_1)


seguid_to_fusion <- readRDS("seguid_to_fusion.RDS")
identical = readRDS("identical.RDS")
offbyone <- readRDS("offbyone.RDS")
ec_pathways <- readRDS("EC_pathways.RDS")

test <- identical %>%
	filter(Var1 %in% seguid_to_fusion$fusion_lvl_1) %>%
	filter(Var2 %in% seguid_to_fusion$fusion_lvl_1)

test <- offbyone %>%
	filter(Var1 %in% seguid_to_fusion$fusion_lvl_1) %>%
	filter(Var2 %in% seguid_to_fusion$fusion_lvl_1)

#####
fusion_groups <- read_csv("subset_groups.csv")
EC_pathways <- readRDS("EC_pathways.RDS")
seguid_to_ec <- data.table::fread("data/uniprot-pe-exp_seguid_to_ec_mapping.tsv")
subset_data <- data.table::fread("data/fusion_data_subset.tsv")


mydata <- dtplyr::lazy_dt(subset_data)
seguid_to_fusion <- mydata %>%
	filter(seguid %in% seguid_to_ec$seguid) %>%
	select(seguid, fusion_lvl_1) %>%
	unique() %>%
	left_join(seguid_to_ec, by = "seguid") %>%
	as.data.table()
saveRDS(seguid_to_fusion,"seguid_to_fusion.RDS")

ec_groups <- fusion_groups %>% 
	mutate(fusion_lvl_1 = as.character(fusion_lvl_1)) %>%
	left_join(seguid_to_fusion, by = "fusion_lvl_1") %>%
	filter(!is.na(seguid))

###########

library(parallelDist)

binDTM <- readRDS("data/subsample_binDTM.RDS")
binDTM_matrix <- as.matrix(binDTM)

message(Sys.time(), " starting distance calculation")
dist <- parallelDist(binDTM_matrix, distance = "manhattan") 
message(Sys.time(), " stopping distance calculation")

saveRDS(dist, "subsample_dist.RDS")

quit()