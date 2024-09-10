
# linkage disequilibrium tutorial -----------------------------------------------------------------

# tutorial https://grunwaldlab.github.io/Population_Genetics_in_R/Linkage_disequilibrium.html----------------------------------------------------------------

#This test is useful to determine if populations are clonal (where significant disequilibrium is expected due to linkage 
#among loci) or sexual (where linkage among loci is not expected). The null hypothesis tested is that alleles observed at 
#different loci are not linked if populations are sexual while alleles recombine freely into new genotypes during the process
#of sexual reproduction. In molecular ecology we typically use the index of association or related indices to test this phenomenon.


# library("poppr")
# library("magrittr")
# data("Pinf") # Load the data
# MX <- popsub(Pinf, "North America")
# ia(MX, sample = 999)
# MX %>% clonecorrect(strata= ~Continent/Country) %>% ia(sample = 999)
# SA <- popsub(Pinf, "South America")
# ia(SA, sample = 999)
# SA %>% clonecorrect(strata= ~Continent/Country) %>% ia(sample=999)
# #Both clone-corrected (N=29) and uncorrected data (N=38) reject the hypothesis of no linkage among markers. We thus have support 
# #for populations in Mexico being sexual while those in South America are clonal.
# mxpair <- MX %>% clonecorrect(strata = ~Continent/Country) %>% pair.ia
# sapair <- SA %>% clonecorrect(strata = ~Continent/Country) %>% pair.ia
# head(mxpair, 10) # Mexico
# head(sapair, 10) # South America
# plotrange <- range(c(mxpair, sapair), na.rm = TRUE)
# plot(mxpair, limits = plotrange)
# plot(sapair, limits = plotrange)



# linkage disequilibrium my data (working)-----------------------------------------------------------------

library("poppr")
library("magrittr")
library("adegenet")


# Calculate index of association
ia(data_genind_adult, sample = 999)   
#likely clones because of reps

#remove adults rep duplicates
# genotypes <- data_genind_adult$other$ind.metrics$genotype
# ind_names <- indNames(data_genind_adult)
# geno_df <- data.frame(individual = ind_names, genotype = genotypes, stringsAsFactors = FALSE)
# unique_geno_df <- geno_df %>% distinct(genotype, .keep_all = TRUE)
# unique_ind_names <- unique_geno_df$individual
# data_genind_unique <- data_genind_adult[unique_ind_names, ]


data_genind_adult_unique %>% clonecorrect() %>% ia(sample = 999)  #check without clones
#the continued sig. effect might indicated that population is clonal structured. 


# Pairwise index of association
pairwise_ia <- data_genind_adult %>% clonecorrect() %>% pair.ia

# Display first 10 pairs
head(pairwise_ia, 10)

# Plot results  (this just standarised the colours amoung populations)
plotrange <- range(pairwise_ia, na.rm = TRUE)
plot(pairwise_ia, limits = plotrange)


# subset by cluster
#by 2
# clusters <- data_genind_unique$other$ind.metrics$cluster
# data_genind_adult_cluster2 <- data_genind_unique[clusters == 2, ]
# data_genind_adult_cluster2$other$ind.metrics$cluster
# data_genind_adult_cluster2$other$ind.metrics$cluster <- droplevels(data_genind_adult_cluster2$other$ind.metrics$cluster)
data_genind_adult_subset1 %>% clonecorrect() %>% ia(sample = 999)  #check without clones
#by 1
# clusters <- data_genind_unique$other$ind.metrics$cluster
# data_genind_adult_cluster1 <- data_genind_unique[clusters == 1, ]
# data_genind_adult_cluster1$other$ind.metrics$cluster
data_genind_adult_subset2 %>% clonecorrect() %>% ia(sample = 999)  #check without clones



# LD using pegas ----------------------------------------------------------

library(adegenet) # Load adegenet package
library(pegas)    # Load pegas package
genotype_data <- tab(data_genind_adult, NA.method = "mean") # Extracting genotype data
genotype_data_df <- as.data.frame(genotype_data) # Convert to data frame
genotype_data_df[] <- lapply(genotype_data_df, factor) # Ensure all data is factor
biallelic_loci <- sapply(genotype_data_df, function(x) length(unique(x[!is.na(x)])) == 2) # Identify biallelic loci
genotype_data_clean <- genotype_data_df[, biallelic_loci] # Filter data
genotype_data_pegas <- as.loci(genotype_data_clean) # Convert to loci format
ld_results <- pegas::LD(genotype_data_pegas) # Calculate LD
correlation_matrix <- ld_results$`Correlations among alleles` # Extract correlations
correlation_matrix <- as.matrix(correlation_matrix) # Convert to matrix
str(correlation_matrix) # Check structure
heatmap(correlation_matrix, main = "Linkage Disequilibrium Heatmap", xlab = "Loci", ylab = "Loci") # Plot heatmap



# number of alleles per locus ---------------------------------------------

# Extract genotype data from genind object
genotype_data <- tab(data_genind_adult, NA.method = "mean")

# Convert to data frame
genotype_data_df <- as.data.frame(genotype_data)

# Check the number of unique alleles per locus
num_alleles_per_locus <- apply(genotype_data_df, 2, function(x) length(unique(na.omit(x))))

# Filter loci that have exactly two alleles
biallelic_loci <- num_alleles_per_locus == 2

# Summary of number of alleles per locus
summary(num_alleles_per_locus)

# Filter the genotype data to retain only biallelic loci
genotype_data_biallelic <- genotype_data_df[, biallelic_loci]

# Check the structure of the filtered data
str(genotype_data_biallelic)



# clones ------------------------------------------------------------------

# Perform MLG analysis
mlg_analysis <- mlg(data_genind, quiet = FALSE)
#63 distinct individuals - note that this includes indivd reps

#Compare relatedness on single indivdual vs all
data_genind_adult@pop <- factor(rep("population1", nrow(data_genind_adult@tab))) #combine all
genetic_dist_matrix <- gd.smouse(data_genind_adult, verbose = TRUE)  #Smouse and Peakall (1999) is a method used to quantify the
genetic_dist_df <- as.data.frame(as.matrix(genetic_dist_matrix))
genetic_dist_df <- tibble::rownames_to_column(genetic_dist_df, "Individual1")
adult_colonies <- pivot_longer(genetic_dist_df, cols = -Individual1, names_to = "Individual2", values_to = "Distance") %>% data.frame()
adult_colonies_sort <- adult_colonies %>% arrange(Individual1, Distance)
hist(adult_colonies_sort$Distance, main = "Genetic Distance Distribution", xlab = "Genetic Distance", ylab = "Frequency")

first_group_data <- adult_colonies_sort %>%
  filter(Individual1 == first_group) %>%
  mutate(Individual2 = factor(Individual2, levels = Individual2[order(Distance)]))
p1 <- ggplot(first_group_data, aes(x = Individual2, y = Distance)) +
  geom_point() +
  facet_wrap(~ Individual1, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Genetic Distances for Adult Colonies", x = "Individual2", y = "Genetic Distance")
p1  #some error between reps

# Calculate Nei's distance for additional insights (doesnt work)
nei_dist <- nei.dist(data_genind,warning = TRUE)
nei_dist_df <- as.data.frame(as.matrix(nei_dist))
print(nei_dist_df)

# Calculate Roger's distance for additional insights (deosnt work)
rogers_dist <- rogers.dist(data_genind)
rogers_dist_df <- as.data.frame(as.matrix(rogers_dist))
print(rogers_dist_df)



# MLL ---------------------------------------------------------------------
# Identify and analyze multilocus lineages (MLLs)
(mll_data <- mll(data_genind_adult))
length(mll_data)

# Calculate allelic richness and observed/expected heterozygosity  (not working)
# Compare genetic diversity before and after clone correction
devtools::install_github("kkeenan02/diveRsity", upgrade = 'never')  #
genetic_diversity <- diveRsity::divBasic(data_genind_adult)
print(genetic_diversity)






