# Set working directory
setwd("E:/Projects/pamfun")

# Load necessary Libraries
library(data.table)
library(tidyverse)
library(dtplyr)
library(Matrix)
rm(list = ls())

#Subset larger data into working set (once)
fusion_data <- data.table::fread("data/fusion_data.tsv")
# Sample 1mil rows
#ind <- sample(1:nrow(fusion_data), 1000000)
#working <-  fusion_data[ind,]
#readr::write_tsv(working, "fusion_data_working2.tsv")
#rm(fusion_data)
# Read in working data
#mydata <- readr::read_tsv("data/fusion_data_working2.tsv")

acc_to_pfam <- fread("data/fusion_pfam_scan_matches_domain_arrangement_preserved.tsv")
seguid_to_ec <- fread("data/uniprot-pe-exp_seguid_to_ec_mapping.tsv")
taxa_ref <- readr::read_csv("taxa_reference.csv")

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
# Combining Phylogeny information

u1 <- mydata2 %>%
  left_join(taxa_ref)
write_csv(as.data.table(u1), "fusion_taxa_data.csv")

# Summary Statistics with and without singletons
###################

# Number of unique Fusion clusters 1335972
length(unique(mydata$fusion_lvl_1))

# Number of unique organism_taxid 5280
length(unique(mydata$organism_taxid))

# Number of unique proteins 31566498
length(unique(mydata$ncbi_accession)) 

# Number of unique proteins by name 745118
length(unique(mydata$protein_name))

# Number of unique Assemblies 8906
length(unique(mydata$assembly_accession))


# Number of unique organisms in a cluster
#########################################

# Group by fusion cluster and count unique number of proteins.
t1 <- mydata2 %>%
  group_by(fusion_lvl_1) %>%
  summarize(n = length(unique(ncbi_accession))) %>%
  arrange(desc(n)) %>%
  rename(n_proteins = "n") %>%
  as.data.table()
t1

head(t1)
write_csv(t1, "outputs/fusion_cluster_size_organisms.csv")

# Calculate Frequency of unique number of proteins by cluster.
t1a <- as.data.frame(table(t1$n_proteins)) %>%
  rename("prots_in_cluster" = Var1, "n" = Freq) %>%
  mutate(prots_in_cluster = as.numeric(as.character(prots_in_cluster)))
head(t1a)

# What do the largest clusters comprise of function wise?
t1b <- mydata2 %>%
  filter(fusion_lvl_1 %in% head(t1)$fusion_lvl_1) %>%
  select(protein_name, fusion_lvl_1) %>%
  unique() %>%
  arrange(fusion_lvl_1, protein_name) %>%
  as.data.table()
head(t1b)

# Filter out singletons
not_single <- filter(t1, n_proteins > 1) %>% pull(fusion_lvl_1)
t1c <- mydata2 %>%
  filter(fusion_lvl_1 %in% not_single)

t1c

# Plot on Log Scale

#png("outputs/histogram_proteins_in_clusters.png", height = 720, width = 1080)
pdf("outputs/histogram_proteins_in_clusters.pdf", height = 8, width = 10)
ggplot(t1, aes( x = n_proteins)) + 
  geom_histogram() +
  scale_y_log10() +
  labs(title = "Histogram for # of proteins in a cluster") +
  xlab("Number of Proteins") + ylab("Frequency") 
dev.off()

# Print Summary Statistics
message("Summary Statistics for # of # of Proteins in a cluster")
print(summary(t1$n))

message("Summary Statistics for # of # of Proteins in a cluster [no singletons]")
print(summary(filter(t1, n_proteins != 1)%>% pull(n_proteins)))

# Number of unique organisms in a cluster?
#########################################

# Group by fusion cluster and count unique number of organisms.
t2 <- mydata2 %>%
  group_by(fusion_lvl_1) %>%
  summarize(n = length(unique(assembly_accession))) %>%
  arrange(desc(n)) %>%
  rename(n_organisms = "n") %>%
  as.data.table()
head(t2)
write_csv(t2, "outputs/fusion_cluster_by_organism.csv")

# Group by fusion cluster and count unique number of organisms (no singletons)
t2c <- t2 %>%
  filter(fusion_lvl_1 %in% not_single) %>%
  as.data.table()
head(t2c)

write_csv(t2, "fusion_cluster_size_proteins.csv")

# Calculate Frequency of unique number of organisms by cluster.
t2a <- as.data.frame(table(t2$n_organisms)) %>%
  rename("orgs_in_cluster" = Var1, "n" = Freq) %>%
  mutate(orgs_in_cluster = as.numeric(as.character(orgs_in_cluster))) %>%
  arrange(desc(n))

# Plot on Log Scale
pdf("outputs/histogram_assemblies_in_clusters.pdf", height = 8, width = 10)
ggplot(t2, aes(x = n_organisms)) + 
  geom_histogram() +
  scale_y_log10() +
  labs(title = "Histogram for # of Unique Assemblies in a cluster") +
  xlab("Number of Assemblies") + ylab("Frequency")
dev.off()

# Print Summary Statistics
message("Summary Statistics for # of Organisms in a cluster")
print(summary(t2$n_organisms))

message("Summary Statistics for # of Organisms in a cluster [no singletons]")
summary(t2c$n_organisms)

# Is there a relationship between # of organisms and # or proteins in a cluster?
eq <- function(x,y) {
  m <- lm(y ~ x)
  as.character(as.expression(
    substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                 list(a = format(coef(m)[1], digits = 4),
                      b = format(coef(m)[2], digits = 4),
                      r2 = format(summary(m)$r.squared, digits = 3)))))
  }


t3 <- t1 %>%
  left_join(t2, by = "fusion_lvl_1") %>% 
  as.data.table()

options(bitmapType='cairo')  
x11(type='cairo')

t3p <- t3 %>%  
  ggplot(aes(x = n_organisms, y = n_proteins)) +
  geom_point(shape = ".") + geom_smooth(method = "lm") + 
  labs(title = "Scatterplot of # of organisms by # of proteins in a cluster") +
  xlab("# of Organisms") + ylab("# of Proteins") +
  geom_text(x = 6000, y = 300, label = eq(t3$n_organisms, t3$n_proteins), parse = TRUE)

pdf("outputs/scatter_organism_protein.pdf", width = 8, height = 6)
scattermore::scattermoreplot(x = t3$n_organisms, y = t3$n_proteins,
                             xlab = "# of Assembly Accessions", ylab = "# of NCBI Accessions",
                             main = "Clusters plotted by their Assembly vs NCBI Accessions.")
dev.off()

cor(t3$n_organisms, y = t3$n_proteins, method = "pearson")


# What clusters are common across all organisms?
################################################

# clusters, # of organisms with that cluster, % coverage of that cluster over all organisms.
t4 <- mydata3 %>%
  select(fusion_lvl_1, assembly_accession) %>%
  unique() %>%
  pull(fusion_lvl_1) %>%
  table() %>%
  as.data.frame() %>%
  rename(fusion_lvl_1 = ".", freq = "Freq") %>%
  mutate(perc = freq / length(unique(mydata$assembly_accession))) %>%
  arrange(desc(perc))
head(t4)
write_csv(t4, "outputs/coverage_clusters.csv")
# Clusters with the highest % coverage contain the most functions.
# Clusters with the lowest % contain the least functions (singletons)

# Plot on Log Scale
#png("outputs/histogram_clusters_of_percentcoverage.png")
pdf("outputs/histogram_clusters_of_percentcoverage.pdf")
ggplot(t4, aes(x = perc)) + 
  geom_histogram() +
  scale_y_log10() +
  labs(title = "Histogram for # of Clusters of a Percent Coverage") + 
  xlab("Percent Coverage") + ylab("Frequency")
dev.off()

message("Summary Statistics for % coverage of Clusters")
print(summary(t4$perc))


highcoverage_clusters <- filter(t4, perc > 0.99) %>% pull(fusion_lvl_1)
t4b <- mydata2 %>%
  filter(fusion_lvl_1 %in% highcoverage_clusters) %>%
  select(fusion_lvl_1, protein_name) %>%
  unique() %>%
  filter(protein_name != "") %>%
  arrange(fusion_lvl_1, protein_name) %>%
  as.data.table()
t4b

write_csv(t4b, "outputs/99_coverage_clusters.csv")

# What about clusters that have a high number of organisms relative to their number of unique proteins?
t4 <- t3 %>%
  mutate(ratio = n_organisms / n_proteins) %>%
  arrange(desc(ratio)) %>%
  as.data.table()
t4  

write_csv(t4, "outputs/cluster_ratio.csv")

# Clusters with a high organism:protein ratio (many organisms, few unique proteins) 
# contain housekeeping genes including DNA repair and ribosomal subunits
a <- t4[1:30,]$fusion_lvl_1
t4a <- mydata2 %>%  
  filter(fusion_lvl_1 %in% a) %>%
  select(fusion_lvl_1, protein_name) %>%
  unique() %>%
  mutate(protein_name = tolower(protein_name), n = 1) %>%
  group_by(fusion_lvl_1, protein_name) %>%
  summarize(n = sum(n)) %>%
  arrange(fusion_lvl_1, protein_name, desc(n)) %>%
  as.data.table()
t4a
write_csv(t4a, "outputs/high_ratio_clusters.csv")

# Clusters with a low organism:protein ratio (few organisms, many unique proteins) 
# contain more unique functions.
t4b <- filter(mydata, fusion_lvl_1 %in% tail(t4,n = 30)$fusion_lvl_1) %>%
  select(fusion_lvl_1, protein_name) %>%
  unique() %>%
  mutate(protein_name = tolower(protein_name), n = 1) %>%
  group_by(fusion_lvl_1, protein_name) %>%
  summarize(n = sum(n)) %>%
  arrange(fusion_lvl_1, protein_name, desc(n)) %>%
  as.data.table()

write_csv(t4b, "outputs/low_ratio_clusters.csv")


pdf("outputs/histogram_ratio.pdf", height = 6, width = 8)
ggplot(t4, aes(x = ratio)) + geom_histogram() +
  scale_y_log10() +
  labs(title = "Histogram of organism:protein ratio for all clusters") +
  xlab("Organism:Protein Ratio") + ylab("Frequency")
dev.off()


# How many clusters is an organism spread over?
t5 <- mydata2 %>%
  select(assembly_accession, fusion_lvl_1) %>%
  unique() %>%
  mutate(s = ifelse(grepl("S", fusion_lvl_1), 1, 0)) %>%
  group_by(assembly_accession) %>%
  summarize(n = length(unique(fusion_lvl_1)), n_s = sum(s)) %>%
  arrange(desc(n)) %>%
  mutate(perc_single = n_s/n) %>%
  as.data.table() 
head(t5)


pdf("outputs/histogram_clusters_per_assembly.pdf", height = 6, width = 8)
ggplot(t5, aes(x = n)) + 
  geom_histogram() +
  scale_y_log10() +
  labs(title = "Histogram of the number of clusters assigned to an assembly") +
  xlab("Number of Clusters") + ylab("Frequency")
dev.off()

pdf("outputs/histogram_single_clusters_per_assembly.pdf", height = 6, width = 8)
ggplot(t5, aes(x = n_s)) + 
  geom_histogram() +
  scale_y_log10() +
  labs(title = "Histogram of the number of single-clusters assigned to an assembly") +
  xlab("Number of Clusters") + ylab("Frequency")
dev.off()

pdf("outputs/histogram_multi_clusters_per_assembly.pdf", height = 6, width = 8)
ggplot(t5, aes(x = n - n_s)) + 
  geom_histogram() +
  scale_y_log10() +
  labs(title = "Histogram of the number of multi-clusters assigned to an assembly") +
  xlab("Number of Clusters") + ylab("Frequency")
dev.off()

# Does number of assemblies correlate with number of clusters?
t6 <- mydata2 %>%
  select(organism_name, fusion_lvl_1, organism_taxid) %>%
  unique() %>%
  group_by(organism_name, organism_taxid) %>%
  summarize(n_clusters = length(unique(fusion_lvl_1))) %>%
  as.data.table()
head(t6)

t7 <- mydata2 %>%
  select(assembly_accession, organism_name, organism_taxid) %>%
  unique() %>%
  group_by(organism_name, organism_taxid) %>%
  summarize(n_assemblies = length(unique(assembly_accession))) %>%
  left_join(t6, by = c("organism_name", "organism_taxid")) %>%
  as.data.table()
head(t7)

pdf("outputs/scatterplot_clusters_assemblies.pdf", height = 6, width = 8)
ggplot(t7, aes(y = n_clusters, x = n_assemblies, group = n_assemblies)) +
  geom_boxplot() + 
  scale_x_log10() +
  labs(title = "Plot of # of assigned clusters by # of assemblies") +
  xlab("# of Assemblies") + ylab("# of Clusters")
dev.off()

lm(n_clusters ~ n_assemblies, t7)
cor(t7$n_assemblies,t7$n_clusters, method = "pearson")
### 

t9 <- mydata2 %>%
  select(organism_name, fusion_lvl_1) %>%
  unique() %>%
  mutate(s = ifelse(grepl("S", fusion_lvl_1), 1, 0)) %>%
  group_by(organism_name) %>%
  summarize(n = length(unique(fusion_lvl_1)), n_s = sum(s)) %>%
  arrange(desc(n)) %>%
  mutate(perc_single = n_s/n) %>%
  as.data.table() 
t9

write_csv(t9, "outputs/num_clusters_organism.csv")


pdf("outputs/histogram_clusters_per_assembly.pdf", height = 6, width = 8)
ggplot(t9, aes(x = n)) + 
  geom_histogram() +
  scale_y_log10() +
  labs(title = "Histogram of the number of clusters assigned to an assembly") +
  xlab("Number of Clusters") + ylab("Frequency")
dev.off()

pdf("outputs/histogram_single_clusters_per_assembly.pdf", height = 6, width = 8)
ggplot(t5, aes(x = n_s)) + 
  geom_histogram() +
  scale_y_log10() +
  labs(title = "Histogram of the number of single-clusters assigned to an assembly") +
  xlab("Number of Clusters") + ylab("Frequency")
dev.off()

pdf("outputs/histogram_multi_clusters_per_assembly.pdf", height = 6, width = 8)
ggplot(t5, aes(x = n - n_s)) + 
  geom_histogram() +
  scale_y_log10() +
  labs(title = "Histogram of the number of multi-clusters assigned to an assembly") +
  xlab("Number of Clusters") + ylab("Frequency")
dev.off()


# What is the EC composition of a cluster?

#ec_numbers <- readr::read_csv("data/ec_numbers2.csv", col_names = TRUE)

# Read in EC class data
ec_classes <- readRDS("data/ec_classes.RDS")

class= ec_classes[[1]]
subclass= ec_classes[[2]]
subsubclass = ec_classes[[3]]

t10 <- mydata2 %>%
  select(seguid, fusion_lvl_1) %>%
  unique() %>% 
  left_join(seguid_to_ec) %>%
  mutate(subsubclass = paste(ec_number, " ", sep = "")) %>%
  mutate(subsubclass = gsub(".[n|0-9]+ ", ".-", subsubclass)) %>%
  mutate(subclass = paste(ec_number, " ", sep = "")) %>%
  mutate(subclass = gsub(".[0-9]+.[n|0-9]+ ", ".-.-", subclass)) %>%
  mutate(class = paste(ec_number, " ", sep = "")) %>%
  mutate(class = gsub(".[0-9]+.[0-9]+.[n|0-9]+ ", ".-.-.-", class)) %>%
  left_join(class) %>%
  left_join(subclass) %>%
  left_join(subsubclass) %>%
  select(seguid, fusion_lvl_1, class, subclass, subsubclass, ec_number, class_desc, subclass_desc, subsubclass_desc)

# EC class of clusters
t10a <- t10 %>%
  group_by(fusion_lvl_1, class, class_desc) %>%
  summarize(n = n()) %>% 
  left_join(t1) %>%
  arrange(desc(n_proteins), fusion_lvl_1) %>%
  as.data.table()

head(t10a)

# EC subclass of clusters
t10b <- t10 %>%
  group_by(fusion_lvl_1, subclass, subclass_desc) %>%
  summarize(n = n()) %>%
  as.data.table()
head(t10b)

# EC subsubclass of clusters
t10c <- t10 %>%
  group_by(fusion_lvl_1, class, class_desc) %>%
  summarize(n = n()) %>%
  as.data.table()
head(t10c)

#
test <- mydata2 %>%
  select(organism_name, fusion_lvl_1) %>%
  unique() %>%
  mutate(n = 1) %>%
  Matrix.utils::dMcast(organism_name + n ~ fusion_lvl_1, drop = FALSE)
class(test)

test <- mydata2 %>%
  select(seguid, fusion_lvl_1) %>%
  unique() %>% 
  as.data.table()

length(unique(test$seguid))
length(unique(test$fusion_lvl_1))
