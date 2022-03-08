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


# Identify what pathway each EC number participates in.
# --------------------------- 
# Read in data
seguid_to_ec <- read_tsv("data/uniprot-pe1-exp-ec.extended_by_hfsp.mapping") 

ec_queries <- seguid_to_ec %>%
	pull(ec_number) %>%
	unique()

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



# First list of all pathways known EC annotations participate in.
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

# Look up each pathway and identify the modules. 
# --------------------------- 

# create pathwayinfoults list to store pathway information
pathwayinfo_list <- list()
pathwayinfo_names <- list()

# loop through each pathway
for(i in 1:length(pathways)){
	message(Sys.time(), pathways[i], " ", length(pathways)-i)
	temp_info <- fetch_term(pathways[i],base = "https://www.kegg.jp/entry/pathway+", term = "Name\\n|Class\\n|Module\\n|Enzyme\\n")
	pathwayinfo_list <- append(pathwayinfo_list, list(temp_info))
	pathwayinfo_names <- append(pathwayinfo_names, pathways[i])
}
names(pathwayinfo_list) <- pathwayinfo_names
# remove duplicates
pathwayinfo_list <- pathwayinfo_list[!duplicated(pathwayinfo_list)]
# remove NA
pathwayinfo_list <- pathwayinfo_list[!is.na(pathwayinfo_list)]
ind <- unlist(lapply(pathwayinfo_list, length)) == 4
pathwayinfo_list <- pathwayinfo_list[ind]

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

parsed_res <- lapply(pathwayinfo_list, custom_parse)

ecs <- lapply(parsed_res, function(x) x[[1]])
pathway_info <- lapply(parsed_res, function(x) x[[2]])
names(pathway_info) <- unlist(ecs)
saveRDS(pathway_info, "pathway_info.RDS")



# Calculate the two way combinations of every enzyme for each pathway.
ec_enzymes <- lapply(info, function(x)x["enzymes"]) %>%
	bind_rows(.id = "ec_id") %>%
	filter(enzymes != "") %>%
	rename(enzymes_a = "enzymes") %>%
	mutate(enzymes_b = enzymes_a) %>%
	group_by(ec_id) %>%
	complete(enzymes_a, enzymes_b)

# optional bit to remove reverse order interactions Ex: A x B, and B x A 
ec_enzymes_one <- ec_enzymes %>%
 	group_by(grp = paste(pmax(enzymes_a, enzymes_b), pmin(enzymes_a, enzymes_b), sep = "_")) %>%
 	slice(1) %>%
 	ungroup() %>%
 	select(-grp) %>%
 	arrange(ec_id)

# Look up each module and identify KOs
# --------------------------- 

# filter pathway info to just bacterial pathways
paths <- readRDS("data/protist_pathways.RDS") %>%
	bind_rows() %>%
	pull(num) %>%
	unique() 
bacteria_ecs <- paste("ec", paths, sep = "")
ind <- names(pathway_info) %in% bacteria_ecs

pathway_info <- pathway_info[ind]
# extract module numbers from pathway information
modules <- lapply(pathway_info, function(x){as.data.frame(unlist(x["module"]))})
modules_classes <- lapply(pathway_info, function(x){paste0(unlist(x["class"]),  collapse = "_")})
modules_names <- lapply(pathway_info, function(x){unlist(x["name"])})
modules_labels <- mapply(function(x,y){paste(x, y, sep = "_")}, x = modules_names, y = modules_classes)
names(modules) <- modules_labels
modules_df <- bind_rows(modules, .id = "names_class") 
colnames(modules_df) = c("names_class", "desc") 
modules_df <- modules_df %>%
	separate(names_class, into = c("pathway", "class", "subclass"), sep = "_") %>%
	separate(desc, into = c("module_id", "desc"), sep = "  ") %>%
	unique() %>% 
	filter(!is.na(desc)) %>%
	select(class, subclass, pathway, module_id, desc)
rownames(modules_df) <- 1:nrow(modules_df)

#
write_csv(modules_df, "data/module_info.csv")
#

module_ids <- modules_df$mo_id
# Match KOs to Enzymes
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
colnames(modules_ko_df) <- c("ko_head", "class", "children")

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

modules_ec <- modules_to_ko %>%
	select(ko_id, ec_id) %>%
	separate(ko_id, into = c("ko_id", "ec_desc"), sep = "  ") %>%
	left_join(modules_ko) %>%
	left_join(modules_df) %>%
	rename(mo_desc = "desc") %>%
	unique() %>%
	filter(!is.na(module_id)) %>%
	rename(ec_number = ec_id) %>%
	mutate(ec_number3 = gsub("\\.[0-9]+\\Z", ".-", ec_number, perl = TRUE))
saveRDS(modules_ec, "modules_ec.RDS")

ec_function <- modules_ec %>%
	select(ko_id, ec_desc, ec_number, ec_number3) %>%
	unique() %>%
	mutate(ec_desc = gsub("\\[|\\]|EC:", "", ec_desc)) %>%
	mutate(ec_desc = trimws(gsub("[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+", "", ec_desc, perl = TRUE)))
write_csv(ec_function, "ec_function.csv")

module_enzymes <- modules_ec %>%
	select(ec_number, module_id) %>%
	rename(enzymes_a = "ec_number") %>%
	mutate(enzymes_b = enzymes_a) %>%
	group_by(module_id) %>%
	complete(enzymes_a, enzymes_b) %>%
	unique()


enzyme_interactions <- list(ec_enzymes, module_enzymes)
names(enzyme_interactions) <- c("pathway", "module") 
enzyme_interactions <- bind_rows(enzyme_interactions, .id = "group") 
write_csv(enzyme_interactions, "data/enzyme_interactions.csv")
# Identify EC Pathways specific to protists
# --------------------------- 
# read in Taxa file from KEGG
kegg_taxa <- fromJSON(file= "data/br08611.json")

# Separate out Bacter and Archae sections
bacteria <- kegg_taxa[[2]][[5]]

# function to extract information from list format to dataframes
custom_extract <- function(.x){
  res <- data.frame(matrix(unlist(.x), ncol=1, byrow=T))
  return(res)
  }

# list of dataframes
bact_df <- lapply(bacteria, custom_extract)[[2]]
colnames(bact_df) = "name"
# name lists

# convert named lists to dataframes
bact_df <- bact_df %>%
 	separate(name, into = c("abr", "name"), sep = "  ") %>%
 	filter(!is.na(name))

# get bacterial abbreivations
abbr <- unique(pull(bact_df, abr))

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

saveRDS(paths, "data/protist_pathways.RDS")
saveRDS(assemblies, "data/protist_geneinfo.RDS")

# Read in data
seguid_to_ec <- read_tsv("data/uniprot-pe1-exp-ec.extended_by_hfsp.mapping") 

ec_queries <- seguid_to_ec %>%
	pull(ec_number) %>%
	unique()


quit()




