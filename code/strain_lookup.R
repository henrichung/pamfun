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