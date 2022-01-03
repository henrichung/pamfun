library(taxonomizr)
library(data.table)
library(tidyverse)

# Read in Fusion Data
data = data.table::fread("fusion_data.tsv") 

# Reshape 
data2 <- select(data, c("organism_name", "organism_taxid", "species_taxid", "assembly_accession")) %>%
  mutate(is_strain = ifelse("organism_taxid" == "species_taxid", F, T)) %>%
  unique() %>%
  mutate(organism_taxid = as.character(organism_taxid))

# Prepare Taxa Accession database.
prepareDatabase('accessionTaxa.sql')

# Retrieve taxid's
taxaIDs = unique(data$organism_taxid)
# Look up Taxonomy
temp = getTaxonomy(taxaIDs,'accessionTaxa.sql')

# Reshape taxonomy data from list to data.frame
temp2 <- tibble::rownames_to_column(as.data.frame(temp), "organism_taxid") %>%
  mutate(organism_taxid = gsub("[ ]+", "", organism_taxid))

# Join taxonomy data to fusion data 
newdata <- left_join(data2, temp2, by = "organism_taxid") %>%
  rowwise() %>%
  mutate(strain = gsub(pattern = species, replacement = "", x = organism_name)) %>%
  rowwise() %>%
  mutate(strain = trimws(strain))

# Write to file
readr::write_csv(newdata, "taxaid_ref.csv")

#####
 
library(rentrez)

# Read in fusion data with Phylum to species taxonomy annotations
tax_data <- read_csv("data/taxaid_ref.csv")

# Custom function to pull out "infraspecieslist" information with rentrez
entrez_infraspecies <- function(x){
  a <- entrez_search(db = "assembly", term = x)
  b <- entrez_summary(db = "assembly", id = a$ids)
  c <- b$biosource$infraspecieslist
  return(c)
}
# Pull out assembly accession numbers
accs <- as.list(tax_data$assembly_accession)

# apply infraspecies function to accessions (2-3 hours~)
strain_list <- lapply(accs, entrez_infraspecies)
names(strain_list) <- accs

# Convert list to dataframe
strains <- bind_rows(strain_list, .id = "assembly_accession")

#Join strain information to fusion w/ taxonomy annotations.
tax_data_strains <- left_join(tax_data, strains, by = "assembly_accession")

# Write to file
write_csv(tax_data_strains, "outputs/taxa_reference.csv")

# Write entries without strain to file
write_csv(filter(tax_data_strains, is.na(sub_value)), "outputs/taxa_reference_NA.csv")


#####

library(tidyverse)

ec_numbers <- read_csv("data/ec_numbers.csv", col_names = FALSE) %>%
  separate(X1, into = c("ec_number", "desc"), sep = "- ") %>%
  mutate(ec_number = paste(gsub(" ", "", ec_number), "-", sep = "")) %>%
  mutate(desc = trimws(tolower(desc))) 
write_csv(ec_numbers, "data/ec_numbers2.csv")


class <- filter(ec_numbers, str_detect(ec_number, '\\.-\\.-\\.-')) %>%
  rename(class_desc = "desc", class = "ec_number"); class
subclass <- filter(ec_numbers, str_detect(ec_number, '[0-9]+\\.[0-9]+\\.-\\.-')) %>%
  rename(subclass_desc = "desc", subclass = "ec_number"); subclass
subsubclass  <- filter(ec_numbers, str_detect(ec_number, '[0-9]+\\.[0-9]+\\.[0-9]+\\.-')) %>%
  rename(subsubclass_desc = "desc", subsubclass = "ec_number"); subsubclass

ec_classes <- list(class, subclass, subsubclass)
saveRDS(ec_classes, "data/ec_classes.RDS")
