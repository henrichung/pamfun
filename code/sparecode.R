# Before

quit()
#
a <- filter(p1_data, group == "shared_post")  %>% pull(value) %>% summary() %>% as.matrix()
b <- filter(p1_data, group == "unshared_post") %>% pull(value) %>% summary() %>% as.matrix()
c <- filter(p1_data, group == "random") %>% pull(value) %>% summary() %>% as.matrix()
d <- cbind(a, b, c)
colnames(d) <- c("shared_post", "unshared_post", "random")
write_csv(as.data.frame(d), "Table of Vector Distance between functions.csv")


# look at weird distances

pathways_dist2 %>%
	arrange(desc(value))

# ec00650: 5938, 8005 | look at functions
a = "8005"
b = "5938"
grouped_fusion <- readRDS("grouped_fusion_all.RDS") 
filter(grouped_fusion, fusion_lvl_1 == a | fusion_lvl_1 == b)

# look at taxa representation
taxa_ref <- read_csv("data/taxa_reference.csv")
filter(taxa_ref)

# identify what assemblies they participate in
reduced_binDTM <- readRDS("outputs/4_reduced_binDTM1.RDS")
int <- reduced_binDTM[c("8005", "5938"),]

a_assemblies <- names(int)[int[1,] == 1] 
b_assemblies <- names(int)[int[2,] == 1] 

a_taxa <- filter(taxa_ref, assembly_accession %in% a_assemblies)
b_taxa <- filter(taxa_ref, assembly_accession %in% b_assemblies)

temp <- list(a_taxa, b_taxa)
names(temp) = c(a, b)
write_csv(bind_rows(temp, .id = "fusion_lvl_1"), "test.csv")
#


#
binDTM <- readRDS("subsample_binDTM.RDS")
int <- binDTM[c("2771", "54452"),]

t2771 <- int[1,][int[1,] == 1] %>% names()
t54452 <- int[2,][int[2,] == 1] %>% names()
taxa_ref <- read_csv("data/taxa_reference.csv")

a2771 <- filter(taxa_ref, assembly_accession %in% t2771)
a54452 <-  filter(taxa_ref, assembly_accession %in% t54452)

b <- list(a2771, a54452)
names(b) <- c("2771", "54452")
c <- bind_rows(b, .id = "fusion_lvl_1")
write_csv(c, "test3.csv")
#
# Z test
# Set difference to be tested
d0<-0
# Set standard deviation of sample with status 0
sigma0<-sd(other_dist$value)
# Set standard deviation of sample with status 1
sigma1<-sd(pathways_dist$value)

# Calculate the two means 
mean_status0<-mean(other_dist$value)
mean_status1<-mean(pathways_dist$value)

# Calculate both the sample sizes 
n_status0<-length(other_dist$value)
n_status1<-length(pathways_dist$value)

# Calculate test statistic and two-sided p-value
z<-((mean_status0-mean_status1)-d0)/sqrt(sigma0^2/n_status0+sigma1^2/n_status1)
p_value=2*pnorm(-abs(z))

p_value

custom_clt <- function(.x, .y, n = 100, rep = 1000){
	x_means <- list()
	y_means <- list()
	for(i in 1:rep){
		x_means[[i]] <- mean(sample(.x, size = n, replace = TRUE))
		y_means[[i]] <- mean(sample(.y, size = n, replace = TRUE))
	}
	return(list(x_means, y_means))
}

custom_z <- function(.x, .y, x_sd, y_sd, d0 = 0){
	# Calculate the two sds.
	sigma0<-x_sd
	sigma1<-y_sd

	# Calculate the two means. 
	mean_status0<-mean(.x)
	mean_status1<-mean(.y)

	# Calculate both the sample sizes. 
	n_status0<-length(.x)
	n_status1<-length(.y)

	# Calculate test statistic and two-sided p-value
	z<-((mean_status0-mean_status1)-d0)/sqrt(sigma0^2/n_status0+sigma1^2/n_status1)
	p_value=2*pnorm(-abs(z))
	return(p_value)
}

sample_means <- custom_clt(other_dist$value, pathways_dist$value)
t <- custom_z(unlist(sample_means[[1]]), unlist(sample_means[[2]]), x_sd = sd(other_dist$value), y_sd = sd(pathways_dist$value) )
# simulation analysis?
N = 10000
k = length(pathways_dist$value)
pop = ec_molten$value
pvals <- list()
for(i in 1:N){
	ind = sample(seq_len(length(pop)), size = k, replace = FALSE)
	group1 = pop[ind]
	group2  = pop[-ind]
	temp_means <- custom_clt(group1, group2)
	pvals[[i]] <- custom_z(unlist(temp_means[[1]]), unlist(temp_means[[2]]), x_sd = sd(group1), y_sd = sd(group2))
}

percentile = ecdf(unlist(pvals))
percentile(t)

# Compare against random
random_means <- custom_clt(other_dist$value, rand_molten$value)
t <- custom_z(unlist(random_means[[1]]), unlist(random_means[[2]]), x_sd = sd(other_dist$value), y_sd = sd(rand_molten$value))
#
library(fastcluster)
library(parallelDist)
m = 1618
n = 1618
rand_binDTM <- sample.int(2, m*n, TRUE)-1L; dim(rand_binDTM) <- c(m,n)
message(Sys.time(), " | RAND DISTANCE CALCULATION START")
rand_dist <- as.matrix(parDist(rand_binDTM, method = "binary"))
message(Sys.time(), " | RAND DISTANCE CALCULATION END")

rand_molten <- reshape2::melt(rand_dist) %>%
	mutate(group = "random", ec_id = NA, combo = "NA", Var1 = 1:nrow(.), Var2 = 1:nrow(.), val = value)

#

binDTM <- readRDS("subsample_binDTM.RDS")
ec_binDTM <- binDTM[rownames(binDTM) %in% seguid_to_fusion$fusion_lvl_1, ]
row_sums <- rowSums(ec_binDTM)
Q = K = ncol(ec_binDTM) 

test = combn(row_sums, 2)
test_list <- as.list(data.frame(test))
test_res <- lapply(test_list, function(.x){(.x[1]/Q) * (.x[2]/Q)})

#
res <- list()
for(i in 1:10000){
rand_vector1 <- sample.int(2, 48718, TRUE)-1L;
rand_vector2 <- sample.int(2, 48718, TRUE)-1L;
data <- rbind(rand_vector1, rand_vector2)
res[[i]] <- dist(data, method = "binary")
}




#
other_pop <- other_dist$val
pathways_sample <- pathways_dist$val
sample_n <- length(pathways_sample)
N = 10000
temp_list <- list()
for(i in 1:N){
temp_list[[i]] <- t.test(sample(other_pop, size = sample_n, replace= TRUE), pathways_sample)$p.value
}
ecdf(unlist(temp_list))(mean(pathways_sample)
###
pathways_samples <- list()
other_samples <- list()
for(i in 1:N){
pathways_samples[[i]] <- sample(pathways_sample, size = 1000, replace = TRUE)
other_samples[[i]] <- sample(other_pop, size = 1000, replace = TRUE)
}
t.test(unlist(pathways_sample), unlist(other_samples))
#######
binDTM <- readRDS("subsample_binDTM.RDS")
pca <- prcomp(t(binDTM), center = TRUE,scale. = TRUE)
########
z.test = function(x,mu,popvar){ one.tail.p <- NULL
	z.score <- round((mean(x)-mu)/(popvar/sqrt(length(x))),3) one.tail.p <- round(pnorm(abs(z.score),lower.tail = FALSE),3) cat(" z =",z.score,"\n",
"one-tailed probability =", one.tail.p,"\n",
"two-tailed probability =", 2*one.tail.p )}



 N = K = 500
test <-  matrix(rbinom(N * K, 1, 0.5), ncol = K, nrow = N)
saveRDS(test, "outputs/test_binDTM.RDS")
#

# filter to pathways specific from protists

# Find pathways found in protists from KEGG
protist_pathways <- readRDS("outputs/protist_pathways.RDS")

protist_pathways_df <- protist_pathways %>%
	bind_rows(.id = "abbr")

common <- protist_pathways_df %>%
	unique() %>%
	group_by(num) %>%
	summarize(n = n()) %>%
	arrange(desc(n)) %>%
	filter(num != "")