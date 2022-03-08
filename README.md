# PAMFUN

PAMFUN is a microbial protein function analysis project utilizing FusionDB. The primary goal for this project is to develop methods that can identify pathways in both individual microorganisms and across entire microbiomes. 

## Phylogentic Profile Analysis.

Phylogenetic profiles are a method of protein analysis where each protein can be represented as a N-dimensional binary vector where N represents the number of unique organism genomes. Each N-th element in a vector represents whether a given protein is present in the N-th genome. The reasoning for this approach is that functionally related proteins are likely to have evolved in the same groups of organisms over time, and thus share similiar vectors. In this analysis, we apply phylogenetic profile analysis to FusionDB clusters. Because fusion clusters are groups of sequence similiar proteins believed to have the same function, fusion cluster and fusion function may be used interchangeably.

The full FusionDB dataset (March 2022) consists of 8,907 distance microorganism assemblies and 31,566,498 proteins. The size of this dataset will typically necessitate use of a high performance computing (HPC) cluster and the provided R code is written for HPC usage. The order and usage of the script is explained below.

### Scripts

 - 0_subset.R - This script subsets the full fusion dataset into a desired workable subset. For this analysis, the largest 50k fusion clusters from assmeblies in balanced organism dataset were selected. Functions were further filtered to those which have a EC annotation mapped to them.

 - 0_ec_fusion.R - This scripts annotates fusion clusters to an EC annotation resolved at 4 levels under certain conditions*.

 - 0_fetch_kegg.R - To determine relatedness of fusion functions, information from the Kyoto Encyclopedia of genes and genomes was collected for each unique EC annotation mapped to a fusion function. Pathway information is limited to only pathways which exist in bacteria. This script scrapes this information from the KEGG API. 

 - 1_binDTM.R - Fusion data is presented as a long dataframe with variable information for each unique protein. This script reshapes this data into a N by M matrix where each row N represents a fusion function and each column M represents an organism assembly which contains a protein present in the fusion cluster. Profiles are only created for functions which are represented in >5 assemblies.

 - 2_distance.R - Preliminary analysis was heavily influeced by pairs of functions which were mapped to the same EC annotation but did were not found in the same organisms. These pairs are believed to be orthologous functions, and their phylogenetic profiles are merged. This script first calculates the Jaccard distance between any two functions and merges the profiles of any two functions which are mapped to the same EC annotation and have a distance > 0.80. Pairwise distance is then recalculated between any and all single or merged fusion clusters. This process is repeated until no more pairs of functions satisify merging conditions. As outputs, it returns the final merged phylogenetic profile matrix, distance matrix, a reference table mapping the newly merged clusters to their EC annotation, and a null distance matrix created by randomly permuting the elements of the profile matrix before calculating distance.

 - 3_analyze.R - This script joins the distance data between any two functions with the KEGG information collected previously. Figures are generated to then examine the distribution of distance values between functions involved in the same pathway (defined by KEGG) and further subsets these distances based on specific modules.

### Data

Necessary data/files to run this analysis is listed below.

- balanced_organism_set.treemmer-RTL-0.9.gtdb.incl_taxonomy.tsv - Assemblies from the balanced organism dataset.

- final_fusion_datafile.reduced_organism_set.treemmer.gtdb.top_50k.txt - list of largest 50k fusion clusters.

- uniprot-pe-exp_seguid_to_ec_mapping.tsv - matches each fusion protein to an EC annotation resolved at 4 levels.

- br08902.json - KEGG organisms in taxonomic groups (Jan 2022)

- fusion_data.tsv - full fusion dataset.
