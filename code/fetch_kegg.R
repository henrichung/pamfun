library(rvest)
library(tidyverse)
library(data.table)
rm(list = ls())

message(Sys.time(), ": Reading in Subset Dataset")
subset_data <- read_tsv("data/fusion_data_subset.tsv")

message(Sys.time(), ": seguid to EC")
seguid_to_ec <- fread("data/uniprot-pe-exp_seguid_to_ec_mapping.tsv")
fusion_to_ec <- subset_data %>%
  select(seguid, fusion_lvl_1) %>%
  filter(seguid %in% seguid_to_ec$seguid) %>%
  left_join(seguid_to_ec, by = "seguid") %>%
  as.data.table()

message(Sys.time(), ": Pulling EC queries.")
ec_queries <- unique(gsub("n", "", fusion_to_ec$ec_number))

fetch_term <- function(.x, base, term){
	url = paste0(base, .x, sep ="")
	html <- read_html(url) %>%
		html_nodes("tr") %>%
		html_text()
	search <- grepl(term, html)
	n <- nchar(html[search]) < 1000
	res <- html[search][n]
	if(identical(res, character(0))){
		return(NA)
	}else{
		return(res)
	}
}

parse_results <- function(x, string){
	if(is.na(x)){return(NA)}
	temp <- unlist(strsplit(unlist(x), split = "\n|  "))
	res <- temp[grepl(string, temp)]
	return(res)
}

res_list <- list()
res_names <- list()
while(length(ec_queries) != 0){
	message(ec_queries[1], " ", length(ec_queries))
	Sys.sleep(1)

	temp_pathway <- fetch_term(ec_queries[1],base = "https://www.kegg.jp/entry/", term = "Pathway")
	if(is.na(temp_pathway)){ec_queries <- ec_queries[-1]; next}

	pathway <- parse_results(temp_pathway, string = "^ec[0-9]+")

	temp_enzyme <- fetch_term(pathway,base = "https://www.kegg.jp/entry/pathway+", term = "Enzyme")
	if(is.na(temp_enzyme)){ec_queries <- ec_queries[-1]; next}

	enzyme <- gsub(" ", "", parse_results(temp_enzyme, string = "[0-9]+.[0-9]+.[0-9]+.[0-9]+"))

	ec_queries <- ec_queries[!(ec_queries %in% enzyme)]
	res_list <- append(res_list, list(enzyme))
	res_names <- append(res_names, pathway)
}
names(res_list) <- res_names
saveRDS(stack(res_list), "EC_pathways.RDS")
quit()

##########################


message(Sys.time(), ": Fetching Pathways Links.")
pathway_list <- list()

message(Sys.time(), ": Fetching Pathways Links.")
pathway_list <- list()
for(i in 1:length(ec_queries)){
	message(ec_queries[i])
	Sys.sleep(1)
	pathway_list[i] <- fetch_term(ec_queries[i],base = "https://www.kegg.jp/entry/", term = "Pathway")
}

pathways0 <- lapply(pathway_list, function(.x){parse_results(.x, string = "^ec[0-9]+")})
pathways <- unique(pathways0[!is.na(pathways0)])

message(Sys.time(), ": Fetch Enzymes in a Pathway.")
enzyme_list <- list()
for(i in 1:length(pathways)){
	message(pathways[i])
	Sys.sleep(1)
	enzyme_list[i] <- fetch_term(pathways[i],base = "https://www.kegg.jp/entry/pathway+", term = "Enzyme")
}

enzymes0 <- lapply(enzyme_list, function(.x) {parse_results(.x, string = "[0-9]+.[0-9]+.[0-9]+.[0-9]+")})
enzymes <- lapply(enzymes0, function(.x){gsub(" ", "", .x)})
names(enzymes) <- unlist(pathways)

saveRDS(enzymes, "data/enzyme_pathways.RDS")


