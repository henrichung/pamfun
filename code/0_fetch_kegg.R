## ---------------------------
## Purpose of script: Download relevant information from KEGG
## Author: Henri Chung
## Date Created: 2021-01-12
## ---------------------------
## Notes:
## ---------------------------

# Load necessary libraries
library(rjson)
library(rvest)
library(tidyverse)
library(data.table)
library(stringr)
rm(list = ls())

# Read in data
seguid_to_ec <- read_tsv("data/uniprot-pe1-exp-ec.extended_by_hfsp.mapping") 

ec_queries <- seguid_to_ec %>%
	pull(ec_number) %>%
	unique()

# read in Taxa file from KEGG
kegg_taxa <- fromJSON(file= "data/br08611.json")

# Identify what pathway each EC number participates in.
# --------------------------- 

## Custom functions to parse KEGG
fetch_term <- function(.x, base, term){
	url = paste0(base, .x, sep ="")
	html <- read_html(url) %>%
		html_nodes("tr") %>%
		html_text()
	search <- grepl(term, html)
	temp <- html[search]
	remove <- grepl("LinkDB", temp)
	res <- temp[!remove]
	Sys.sleep(1)
	if(identical(res, character(0))){
		return(NA)
	}else{
		return(res)
	}
}

parse_results <- function(x, string){
	if(is.na(x)){return(NA)}
	res <- unlist(str_extract_all(x, string))
	return(res)
}



# First fetch pathway information from list of matched ECs
pathways_list <- list()
pathways_names <- list()
while(length(ec_queries) != 0){
	
	message(Sys.time(), ec_queries[1], " ", length(ec_queries))
	# Fetch the pathway entry for an EC number
	temp_pathway <- fetch_term(ec_queries[1],base = "https://www.kegg.jp/entry/", term = "Pathway")
	# If there is no pathway information, skip.
	if(is.na(temp_pathway)){ec_queries <- ec_queries[-1]; next}
	
	# append the result to pathways_list
	pathways_list <- append(pathways_list, unique(parse_results(temp_pathway, string = "ec[0-9]{5}")))
	pathways_names <- append(pathways_names, ec_queries[1])
	# removed retrieved ec entry.
	ec_queries <- ec_queries[-1]
}
pathways <- unique(unlist(pathways_list))
saveRDS(pathways, "pathways.RDS")


# Look up each pathway and identify the modules. 
# --------------------------- 

# create results list to store pathway information
res_list <- list()
res_names <- list()

# loop through each pathway
for(i in 1:length(pathways)){
	message(Sys.time(), pathways[i], " ", length(pathways)-i)
	temp_info <- fetch_term(pathways[i],base = "https://www.kegg.jp/entry/pathway+", term = "Name\\n|Class\\n|Module\\n|Enzyme\\n")
	res_list <- append(res_list, list(temp_info))
	res_names <- append(res_names, pathways[i])
}
names(res_list) <- res_names
saveRDS(res_list, "res_list.RDS")

res_list <- readRDS("res_list.RDS")
# remove duplicates
res_list <- res_list[!duplicated(res_list)]
# remove NA
res_list <- res_list[!is.na(res_list)]
ind <- unlist(lapply(res_list, length)) == 4
res_list <- res_list[ind]

# Function to parse the name, class, module, and enzyme information.
custom_parse <- function(.x){
	# extract ec path name.
	ec <- as.character(str_match(.x[3], "ec[0-9]{5}"))
	# name of pathway
	name <- unlist(trimws(gsub("Name\\n", "", .x[1])))
	# identify class of pathway
	class <- unlist(trimws(unlist(str_split(gsub("Class\\n|BRITE hierarchy", "", .x[2]), pattern = ";"))))
	# identify modules of pathway
	module <- unlist(trimws(unlist(str_split(gsub("Module\\n", "", .x[3]), pattern = "\\[PATH:ec[0-9]{5}\\]"))))
	# identify enzymes in pathway
	enzymes <- unlist(trimws(unlist(str_split(gsub("Enzyme\\n", "", .x[4]), pattern = "   "))))
	temp_info <- list(name, class, module, enzymes)
	names(temp_info) <- c("name", "class", "module", "enzymes")
	return(list(list(ec), temp_info))
}

parsed_res <- lapply(res_list, custom_parse)

ecs <- lapply(parsed_res, function(x) x[[1]])
info <- lapply(parsed_res, function(x) x[[2]])
names(info) <- unlist(ecs)
saveRDS(info, "info.RDS")
info <- readRDS("info.RDS")


# Look up each module and identify KOs
# --------------------------- 

# extract module numbers from pathway information
modules <- unlist(lapply(info, function(x){x["module"]}))
modules <- unique(modules[modules != ""])
modules <- str_split(modules, pattern = "  ")
modules_df <- do.call(rbind.data.frame, modules)
colnames(modules_df) <- c("module_id", "desc")

module_ids <- modules_df$module_id
# Match KOs to Enzymes

# create modules to store modules
module_list <- list()
module_names <- list()

# loop over module ids to fetch what KOs participate in them
for(i in 1:length(module_ids)){
	message(Sys.time(), " ", module_ids[i], " ", length(module_ids)-i)
	temp_mod <- fetch_term(module_ids[i],base = "https://www.kegg.jp/module/", term = "Definition")
	module_list[[i]] <- temp_mod
	module_names[[i]] <- module_ids[i]
}
names(module_list) <- module_names


parse_module <- function(x){
	temp <- gsub("[\n()|Definition|Ortholog table|Taxonomy|Module]", "", x)
	temp <- gsub("\\+", ",", temp)
	temp <- gsub("K", ",K", temp)
	temp <- gsub(",,", ",", temp)
	res <- unlist(str_split(temp, pattern = ","))
	return(res[res != ""])
}

modules_ko <- lapply(module_list, parse_module) %>%
	stack() %>%
	mutate(values = gsub("[-]+", "", values)) %>%
	rename(ko_id = "values", module_id = "ind") %>%
	unique()
saveRDS(modules_ko, "modules_ko.RDS")

# Map KO information from modules to the EC number.
modules_ko <- readRDS("modules_ko.RDS")
modules_ko_df <- jsonlite::fromJSON("data/ko00001.json") %>%
	as.data.frame()
colnames(modules_ko_df ) <- c("ko_head", "class", "children")

modules_to_ko <- modules_ko_df  %>% 
	as_tibble() %>%
	unnest(children) %>%
	rename(subclass = "name") %>%
	unnest(children) %>%
	rename(module = "name") %>%
	unnest(children) %>%
	rename(ko_id = "name") %>%
	mutate(ec_id = str_match(ko_id, "\\[EC.*\\]")) %>%
	filter(!is.na(ec_id)) %>% 
	mutate(ec_id = gsub("\\[|\\]|EC:", "", ec_id)) %>% 
	separate(ec_id, into = c("a", "b"), sep = " ") %>%
	pivot_longer(cols = c("a", "b"), names_to = "temp", values_to = "ec_id") %>%
	select(-temp) %>% filter(!is.na(ec_id)) %>%
	unique()
write_csv(modules_to_ko, "modules_to_ko.csv")
# Identify EC Pathways specific to protists
# --------------------------- 
quit()
# Separate out Bacter and Archae sections
bacteria <- kegg_taxa[[2]][[5]]
archea <- kegg_taxa[[2]][[6]]


# function to extract information from list format to dataframes
custom_extract <- function(.x){
  res <- data.frame(matrix(unlist(.x), ncol=1, byrow=T))
  return(res)
  }

# list of dataframes
bact_list <- lapply(bacteria, custom_extract)[[2]]
archae_list <- lapply(archea, custom_extract)[[2]]

# name lists
colnames(bact_list) = "name"
colnames(archae_list) = "name"

# convert named lists to dataframes
bact_df <- bact_list %>%
 	separate(name, into = c("abr", "name"), sep = "  ") %>%
 	filter(!is.na(name))

arch_df <- archae_list %>%
 	separate(name, into = c("abr", "name"), sep = "  ") %>%
 	filter(!is.na(name))

# save bacteria and archae dataframes as single list
kegg_abb <- list(bact_df, arch_df) %>%
	bind_rows(.id = "Kingdom") 

abbr <- pull(kegg_abb, abr)


paths <- list()
assemblies <- list()
for(i in 1:length(abbr)){
		a <- read_html(paste("https://www.genome.jp/kegg-bin/show_organism?menu_type=pathway_maps&org=",abbr[i], sep = ""))%>%
			html_nodes("tr") %>%
			html_text()
		
		b <- str_split(a, "\n") %>%
			unlist() %>%
			as.data.frame()
		colnames(b) = "name"
		
		c <- filter(b, name != "" & name != " "  & name != "  ") %>%
			separate(name, into = c("num", "name"), sep = "  ") %>%
			filter(!is.na(name))
		message(Sys.time(), " ", abbr[i])

		d <- read_html(paste("https://www.genome.jp/kegg-bin/show_organism?menu_type=genome_info&org=",abbr[i], sep = ""))%>%
			html_nodes("tr") %>%
			html_text()
		
		e <- str_split(d, "\n") %>%
			unlist() %>%
			as.data.frame()
		
		assemblies[[i]] <- e
		paths[[i]] <- c
}
names(paths) = abbr
names(assemblies) = abbr

saveRDS(paths, "outputs/1_protist_pathways.RDS")
saveRDS(paths, "outputs/1_protist_geneinfo.RDS")


