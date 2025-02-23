

library(poppr)
library(magrittr)
library(adegenet) # Load adegenet package
library(pegas)    # Load pegas package

















ia(data_genind_adult, sample = 999)   


data_genind_adult_unique %>% clonecorrect() %>% ia(sample = 999)  #check without clones

pairwise_ia <- data_genind_adult %>% clonecorrect() %>% pair.ia

plotrange <- range(pairwise_ia, na.rm = TRUE)
plot(pairwise_ia, limits = plotrange)

data_genind_adult_subset1 %>% clonecorrect() %>% ia(sample = 999)  #check without clones
data_genind_adult_subset2 %>% clonecorrect() %>% ia(sample = 999)  #check without clones


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



genotype_data <- tab(data_genind_adult, NA.method = "mean")
genotype_data_df <- as.data.frame(genotype_data)
num_alleles_per_locus <- apply(genotype_data_df, 2, function(x) length(unique(na.omit(x))))
biallelic_loci <- num_alleles_per_locus == 2
summary(num_alleles_per_locus)
genotype_data_biallelic <- genotype_data_df[, biallelic_loci]
str(genotype_data_biallelic)


mlg_analysis <- mlg(data_genind, quiet = FALSE)

data_genind_adult@pop <- factor(rep("population1", nrow(data_genind_adult@tab))) #combine all
genetic_dist_matrix <- gd.smouse(data_genind_adult, verbose = TRUE)  #Smouse and Peakall (1999) is a method used to quantify the
genetic_dist_df <- as.data.frame(as.matrix(genetic_dist_matrix))
genetic_dist_df <- tibble::rownames_to_column(genetic_dist_df, "Individual1")
adult_colonies <- pivot_longer(genetic_dist_df, cols = -Individual1, names_to = "Individual2", values_to = "Distance") %>% data.frame()
adult_colonies_sort <- adult_colonies %>% arrange(Individual1, Distance)
hist(adult_colonies_sort$Distance, main = "Genetic Distance Distribution", xlab = "Genetic Distance", ylab = "Frequency")

first_group = 'pd1.a.1'
first_group_data <- adult_colonies_sort %>%
  filter(Individual1 == first_group) %>%
  mutate(Individual2 = factor(Individual2, levels = Individual2[order(Distance)]))
p1 <- ggplot(first_group_data, aes(x = Individual2, y = Distance)) +
  geom_point() +
  facet_wrap(~ Individual1, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Genetic Distances for Adult Colonies", x = "Individual2", y = "Genetic Distance")
p1  #some error between reps

nei_dist <- nei.dist(data_genind,warning = TRUE)
nei_dist_df <- as.data.frame(as.matrix(nei_dist))
print(nei_dist_df)

rogers_dist <- rogers.dist(data_genind)
rogers_dist_df <- as.data.frame(as.matrix(rogers_dist))
print(rogers_dist_df)


(mll_data <- mll(data_genind_adult))
length(mll_data)

genetic_diversity <- diveRsity::divBasic(data_genind_adult)
print(genetic_diversity)






