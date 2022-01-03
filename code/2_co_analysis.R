library(Matrix)
library(quanteda)
library(tidyverse)

coocCounts <- readRDS("data/bin_cooc_sp.RDS")
binDTM <- readRDS("data/binDTM.RDS")
colnames(coocCounts) <- colnames(binDTM)
rownames(coocCounts) <- colnames(binDTM)
colsum <- colSums(binDTM)

#binDTM <- readr::read_csv("binDTM.csv") %>%
#	mltools::sparsify() %>%
#	t()


coverage <- read_csv("outputs/coverage_clusters.csv")
shannon <- read_csv("outputs/clusters_shannon.csv") %>%
	rename("fusion_lvl_1" = "ind", "shannon" = "values")

int <- left_join(coverage, shannon) 

filter(int, perc > 0.99) %>% arrange(desc(shannon))
cor(int$perc, int$shannon, use="complete.obs")

# 786
# 14720 pick a node with a random number of connections.
coocTerm <- "14720"
k <- nrow(binDTM)
ki <- sum(binDTM[, coocTerm])
kj <- colSums(binDTM)
names(kj) <- colnames(binDTM)
kij <- coocCounts[coocTerm, ]

#
# log (nrow(binDTM) * coocCounts[coocTerm] / ())
mutualInformationSig <- log(k * kij / (ki * kj))
mutualInformationSig <- mutualInformationSig[order(mutualInformationSig, decreasing = TRUE)]

resultOverView <- data.frame(
  names(sort(kij, decreasing=T)[1:10]), sort(kij, decreasing=T)[1:10],
  names(mutualInformationSig[1:10]), mutualInformationSig[1:10])
colnames(resultOverView) <- c("Freq-terms", "Freq", "MI-terms", "MI")


#
coocTerm <- "14720"
numberOfCoocs <- 50
k <- nrow(binDTM)
ki <- sum(binDTM[, coocTerm])
kj <- colSums(binDTM)
names(kj) <- colnames(binDTM)
kij <- coocCounts[coocTerm, ]

#
# log (nrow(binDTM) * coocCounts[coocTerm] / ())
mutualInformationSig <- log(k * kij / (ki * kj))
mutualInformationSig <- mutualInformationSig[order(mutualInformationSig, decreasing = TRUE)]
coocs <- head(mutualInformationSig, numberOfCoocs)
resultGraph <- data.frame(from = character(), to = character(), sig = numeric(0))
# The structure of the temporary graph object is equal to that of the resultGraph
tmpGraph <- data.frame(from = character(), to = character(), sig = numeric(0))

# Fill the data.frame to produce the correct number of lines
tmpGraph[1:numberOfCoocs, 3] <- coocs[1:numberOfCoocs]
# Entry of the search word into the first column in all lines
tmpGraph[, 1] <- coocTerm
# Entry of the co-occurrences into the second column of the respective line
tmpGraph[, 2] <- names(coocs)[1:numberOfCoocs]
# Set the significances
tmpGraph[, 3] <- coocs[1:numberOfCoocs]

# Attach the triples to resultGraph
resultGraph <- rbind(resultGraph, tmpGraph)

MIS <- function(x){
	k <- nrow(binDTM)
	ki <- sum(binDTM[, x])
	kj <- colSums(binDTM)
	names(kj) <- colnames(binDTM)
	kij <- coocCounts[x, ]

	mutualInformationSig <- log(k * kij / (ki * kj))
	mutualInformationSig <- mutualInformationSig[order(mutualInformationSig, decreasing = TRUE)]
	coocs <- head(mutualInformationSig, numberOfCoocs)
	return(coocs)
}

connect_no <- 2
for(j in 1:connect_no){
# Iteration over the most significant numberOfCoocs co-occurrences of the search term
	for (i in 1:numberOfCoocs){
  		print(i)
  		# Calling up the co-occurrence calculation for term i from the search words co-occurrences
  		newCoocTerm <- names(coocs)[i]
  		coocs2 <- MIS(newCoocTerm)
  		
  		# Structure of the temporary graph object
  		tmpGraph <- data.frame(from = character(), to = character(), sig = numeric(0))
  		tmpGraph[1:numberOfCoocs, 3] <- coocs2[1:numberOfCoocs]
  		tmpGraph[, 1] <- newCoocTerm
  		tmpGraph[, 2] <- names(coocs2)[1:numberOfCoocs]
  		tmpGraph[, 3] <- coocs2[1:numberOfCoocs]
  
  		#Append the result to the result graph
  		resultGraph <- rbind(resultGraph, tmpGraph[2:length(tmpGraph[, 1]), ])
	}
	temp <-  unique(c(resultGraph$from, resultGraph$to))[1:numberOfCoocs]
	names(coocs) <- unique(c(resultGraph$from, resultGraph$to))
	resultGraph <- unique(resultGraph)
}



#
library(igraph)
# set seed for graph plot
set.seed(1)

# Create the graph object as undirected graph
graphNetwork <- graph.data.frame(resultGraph, directed = F)

# Identification of all nodes with less than 2 edges
verticesToRemove <- V(graphNetwork)[degree(graphNetwork) < 2]

# These edges are removed from the graph
graphNetwork <- delete.vertices(graphNetwork, verticesToRemove) 

# Assign colors to nodes (search term blue, others orange)
V(graphNetwork)$color <- ifelse(V(graphNetwork)$name == coocTerm, 'cornflowerblue', 'orange') 

# Set edge colors
E(graphNetwork)$color <- adjustcolor("DarkGray", alpha.f = .5)

# scale significance between 1 and 10 for edge width
E(graphNetwork)$width <- scales::rescale(E(graphNetwork)$sig, to = c(1, 10))

# Set edges with radius
E(graphNetwork)$curved <- 0.15 
# Size the nodes by their degree of networking (scaled between 5 and 15)
V(graphNetwork)$size <- scales::rescale(log(degree(graphNetwork)), to = c(5, 15))


pdf("14720.pdf")
# Define the frame and spacing for the plot
par(mai=c(0,0,1,0)) 

# Final Plot
plot(
  graphNetwork,             
  layout = layout.fruchterman.reingold, # Force Directed Layout 
  main = paste(coocTerm, ' Graph'),
  vertex.label.family = "sans",
  vertex.label.cex = 0.3,
  vertex.shape = "circle",
  vertex.label.dist = 0.1,          # Labels of the nodes moved slightly
  vertex.frame.color = adjustcolor("darkgray", alpha.f = .5),
  vertex.label.color = 'black',     # Color of node names
  vertex.label.font = 0.5,            # Font of node names
  vertex.label = V(graphNetwork)$name,      # node names
  vertex.label.cex = 0.5 # font size of node names
)
dev.off()

# 

#fusion_data <- data.table::fread("data/fusion_data.tsv")

clusters <- unique(c(resultGraph$from, resultGraph$to))

mydata2 <- dtplyr::lazy_dt(fusion_data)

fusion_int <- mydata2 %>% 
	filter(fusion_lvl_1 %in% clusters) %>% 
	select(fusion_lvl_1, protein_name) %>%
	unique() %>% 
	arrange(fusion_lvl_1) %>%
	data.table::as.data.table()
write_csv(fusion_int, "26554.csv")

#1743 - highest prevalence 0.998 | 1.46
#4321  - highest prevalence 0.997 | 1.46
#6076 - highest shannon 0.00797 | 2.52
# 26554 - highest shannon 0.00449 | 2.52

# melt dataframe to long?


test <- mefa4::Melt(coocCounts)
library(data.table)
test_dt <- data.table::as.data.table(test)
test2 <- unique(test_dt)