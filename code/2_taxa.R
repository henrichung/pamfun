# Set working directory
setwd("E:/Projects/pamfun")

# Load necessary Libraries
library(data.table)
library(tidyverse)
library(dtplyr)
library(Matrix)
library(Matrix.utils)
library(quanteda)
library(Rcpp)
rm(list = ls())

#Subset larger data into working set (once)
fusion_data <- data.table::fread("data/fusion_data.tsv")
#acc_to_pfam <- fread("data/fusion_pfam_scan_matches_domain_arrangement_preserved.tsv")
#seguid_to_ec <- fread("data/uniprot-pe-exp_seguid_to_ec_mapping.tsv")
taxa_ref <- readr::read_csv("data/taxa_reference.csv")

mydata = fusion_data
mydata2 <- dtplyr::lazy_dt(mydata)

keep <- mydata2 %>%
  select(ncbi_accession, fusion_lvl_1) %>%
  group_by(fusion_lvl_1) %>%
  summarize(n = n()) %>%
  filter(n > 10) %>%
  pull(fusion_lvl_1)

mydata3 <- mydata2 %>%
  filter(fusion_lvl_1 %in% keep)

# Calculate number of unique taxa levels.
t1 <- mydata3 %>%
  select(assembly_accession, fusion_lvl_1) %>%
  filter(!grepl("s|S", fusion_lvl_1)) %>%
  unique() %>%
  left_join(taxa_ref) 

t2 <- t1 %>%
  group_by(fusion_lvl_1) %>%
  summarize(phylum_n = length(unique(phylum)),
            class_n = length(unique(class)),
            order_n = length(unique(order)),
            family_n = length(unique(family)),
            genus_n = length(unique(genus)),
            species_n = length(unique(genus)),
            assembly_n = length(unique(assembly_accession))) %>%
  arrange(desc(assembly_n)) %>%
  as.data.table()

write_csv(t2, "clusters_by_taxa.csv")

# Pivot Phylum to wide
pivot_taxa <- function(x){
  count <- t1 %>%
    select(eval(x), fusion_lvl_1) %>%
    mutate(val = 1) %>%
    as.data.table() %>%
    data.table::dcast(fusion_lvl_1 ~ get(x), value.var = "val", drop = FALSE, fill = 0) %>%
    column_to_rownames("fusion_lvl_1") %>%
    rename("unknown" = "NA")

  prop <- count / rowSums(count)

  prop_long <- prop %>%
    rownames_to_column("fusion_lvl_1") %>%
    reshape2::melt(id.vars = "fusion_lvl_1") %>%
    filter(value != 0) %>%
    arrange(desc(value))

  ref <- stack(table(taxa_ref[[x]])) %>%
    rename("ref_count" = values)

  res <- list(prop_long, ref)
  return(res)
}
# Phylum
phylum <- pivot_taxa("phylum")
t3_phylum <- phylum[[1]]; ref_phylum <- phylum[[2]]
# Class
class <- pivot_taxa("class")
t3_class <- class[[1]]; ref_class <- class[[2]]
# Order
order <- pivot_taxa("order")
t3_order <- order[[1]]; ref_order <- order[[2]]
# Family
family <- pivot_taxa("family")
t3_family <- family[[1]]; ref_family <- family[[2]]
# Genus
#genus <- pivot_taxa("genus")
#t3_genus <- genus[[1]]; ref_genus <- genus[[2]]
# Species
#species <- pivot_taxa("species")
#t3_species <- species[[1]]; ref_species <- species[[2]]
# Strain
#strain <- pivot_taxa("strain")
#t3_strain <- strain[[1]]; ref_strain <- strain[[2]]

########
filter(t3_order, variable != "unknown") %>%
  filter(value > 0.50) %>%
  pull(fusion_lvl_1) %>%
  unique() %>%
  length()

filter(t3_order, variable != "unknown") %>%
  filter(value > 0.75) %>%
  pull(fusion_lvl_1) %>%
  unique() %>%
  length()

filter(t3_order, variable != "unknown") %>%
  filter(value > 0.90) %>%
  pull(fusion_lvl_1) %>%
  unique() %>%
  length()

filter(t3_order, variable != "unknown") %>%
  filter(value > 0.99) %>%
  pull(fusion_lvl_1) %>%
  unique() %>%
  length()

#family :  97479,  76300,  63248,  48482


#calculate clusters why high shannon diversity
t3_phylum_wide <- t3_phylum %>%
  filter(variable != "unknown") %>%
  pivot_wider(names_from = variable, values_from = value, values_fill = 0) %>%
  column_to_rownames("fusion_lvl_1")

# function to calculate shannon diversity
t3_phylum_shannon <- t3_phylum_wide %>%
  apply(1, function(x){
    y <- x[x!=0]
    z <- -sum(y * log(y))
  })
# remove 0s
t3_phylum_shannon = t3_phylum_shannon[t3_phylum_shannon != 0]

t3_phylum_shannon_long <- t3_phylum_shannon %>%
  sort() %>%
  stack() %>%
  as_tibble() %>%
  arrange(desc(values))
write_csv(t3_phylum_shannon_long, "outputs/clusters_shannon.csv")


test <- table(t3_class_long$variable) %>% 
  sort() %>% 
  rev() %>%
  stack() %>%
  rename("cluster_count" = values) %>%
  full_join(ref_class)

  ####
# Calculate phylum proportions 
t3b <- rowSums(t3a)
t3c <- t3a / t3b

#############
test <- filter(t3_class_long, variable != "unknown")

a <- filter(mydata2, fusion_lvl_1 == 102858) %>%
  left_join(taxa_ref) %>%
  select(class, fusion_lvl_1, protein_name)



# function to calculate shannon diversity
t3d <- t3c %>%
  apply(1, function(x){
    y <- x[x!=0]
    z <- -sum(y * log(y))
  })
# remove 0s
t3e = t3d[t3d != 0]


# Peek at clusters with high shannon diversity
t4 <- names(sort(t3e) %>% head(n = 100))
t4a <- t3a %>%
  rownames_to_column("name") %>%
  filter(name %in% t4)
t4b <- t4a %>%
  reshape2::melt(id.vars = "name") %>%
  ggplot(aes(x = name, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "fill")
t4b

# Peek at clusters with low shannon diversity
t5 <- names(sort(t3e) %>% rev() %>% head(n = 100))
t5a <- t3a %>%
  rownames_to_column("name") %>%
  filter(name %in% t5)
t5a

t5b <- t5a %>%
  reshape2::melt(id.vars = "name") %>%
  ggplot(aes(x = name, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "fill")
t5b

# Calculate most prevalent phylum over threshold

t6a <- t3c %>%
  rownames_to_column("fusion_lvl_1") %>%
  reshape2::melt(id.vars = "fusion_lvl_1")

t6b <- t6a %>%
  filter(value >= 0.5)

t6c <- stack(table(taxa_ref$phylum)) %>%
  rename("ref_count" = values)

t6d <- table(t6b$variable) %>% 
  sort() %>% 
  rev() %>%
  stack() %>%
  rename("cluster_count" = values) %>%
  full_join(t6c)

t6e <- t6d %>%
  ggplot(aes(x = cluster_count, y = ref_count)) + geom_point(size = 0.5) +
  xlab("Number of clusters") + ylab("Number of Assembly Accessions") +
  labs(title = "Correlation between # of over 50% Phylum clusters with # of Assemblies",
       subtitle = "Pearson correlation 0.99")

pdf("taxa_correlation_50.pdf")
t6e
dev.off()

## 
t7a <- t3c %>%
  rownames_to_column("fusion_lvl_1") %>%
  reshape2::melt(id.vars = "fusion_lvl_1")

t7b <- t7a %>%
  filter(value >= 0.7)
length(unique(t7b$fusion_lvl_1))
       
t7c <- stack(table(taxa_ref$phylum)) %>%
  rename("ref_count" = values)

t7d <- table(t7b$variable) %>% 
  sort() %>% 
  rev() %>%
  stack() %>%
  rename("cluster_count" = values) %>%
  full_join(t7c)

cor(t7d$cluster_count, t7d$ref_count, method = "pearson")

t7e <- t7d %>%
  ggplot(aes(x = cluster_count, y = ref_count)) + geom_point(size = 0.5) +
  xlab("Number of clusters") + ylab("Number of Assembly Accessions") +
  labs(title = "Correlation between # of over 50% Phylum clusters with # of Assemblies",
       subtitle = "Pearson correlation 0.99")

pdf("taxa_correlation_70.pdf")
t7e
dev.off()
## 
t8a <- t3c %>%
  rownames_to_column("fusion_lvl_1") %>%
  reshape2::melt(id.vars = "fusion_lvl_1")

t8b <- t8a %>%
  filter(value >= 0.9)
length(unique(t8b$fusion_lvl_1))
       
t8c <- stack(table(taxa_ref$phylum)) %>%
  rename("ref_count" = values)

t8d <- table(t8b$variable) %>% 
  sort() %>% 
  rev() %>%
  stack() %>%
  rename("cluster_count" = values) %>%
  full_join(t7c)

cor(t8d$cluster_count, t8d$ref_count, method = "pearson")

t8e <- t8d %>%
  ggplot(aes(x = cluster_count, y = ref_count)) + geom_point(size = 0.5) +
  xlab("Number of clusters") + ylab("Number of Assembly Accessions") +
  labs(title = "Correlation between # of over 90% Phylum clusters with # of Assemblies",
       subtitle = "Pearson correlation 0.99")

pdf("taxa_correlation_90.pdf")
t8e
dev.off()

###

t8a <- t3c %>%
  rownames_to_column("fusion_lvl_1") %>%
  reshape2::melt(id.vars = "fusion_lvl_1")

t8b <- t8a %>%
  filter(value >= 0.99)
length(unique(t8b$fusion_lvl_1))

t8c <- stack(table(taxa_ref$phylum)) %>%
  rename("ref_count" = values)

t8d <- table(t8b$variable) %>% 
  sort() %>% 
  rev() %>%
  stack() %>%
  rename("cluster_count" = values) %>%
  full_join(t7c)
write_csv(t8d, "clusters_taxa_99.csv")
cor(t8d$cluster_count, t8d$ref_count, method = "pearson")

t8e <- t8d %>%
  ggplot(aes(x = cluster_count, y = ref_count)) + geom_point(size = 0.5) +
  xlab("Number of clusters") + ylab("Number of Assembly Accessions") +
  labs(title = "Correlation between # of over 99% Phylum clusters with # of Assemblies",
       subtitle = "Pearson correlation 0.99")

pdf("taxa_correlation_99.pdf")
t8e
dev.off()


t6 <- t1 %>%
  filter(fusion_lvl_1 %in% names(t3e)) %>%
  group_by(fusion_lvl_1)%>%
  summarize(assembly_n = length(unique(assembly_accession))) %>%
  rename(ind = "fusion_lvl_1") %>%
  left_join(stack(t3e)) %>%
  rename(shannon = "values") %>%
  as.data.table() 

ggplot(t6, aes(x = assembly_n, y = shannon)) + geom_point(pch = ".") +
  scale_x_log10()


t3f <- stack(t3e) %>%
  rename(shannon = "values", fusion_lvl_1 = "ind")

t7 <- t2 %>%
  mutate(phylum_ratio = assembly_n / phylum_n) %>%
  arrange(phylum_ratio) %>%
  filter(phylum_ratio == 1) %>%
  pull(fusion_lvl_1)


t7a <- mydata2 %>%
  filter(fusion_lvl_1 %in% t7) %>%
  select(protein_name, fusion_lvl_1) %>%
  arrange(fusion_lvl_1) %>%
  as.data.table() 
View(t7a)
