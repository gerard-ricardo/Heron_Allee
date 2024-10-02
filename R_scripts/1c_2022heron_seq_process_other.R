# This script is used to create various data outputs to be used in non-R software

# load libraries ----------------------------------------------------------
library(readr)  
library(dplyr)  
library(tidyr)

# id row, meta and loci cols ----------------------------------------------

# Load your raw data
raw_snp_data <- read_csv('./data/tranpose_platy_snp.csv')

colnames(raw_snp_data)[1] <- 'id'  # Rename the first column to 'id'

raw_metadata <- read_csv('./data/meta_platy_ordered.csv')

# Merge data
merged_data <- inner_join(raw_metadata, raw_snp_data, by = 'id')

# Replace missing values and filter the data
snp_columns <- names(merged_data)[9:length(names(merged_data))]
merged_data[snp_columns] <- replace(merged_data[snp_columns], is.na(merged_data[snp_columns]), -99)

# Filter individuals based on missing data
missing_data_percentage <- apply(merged_data[snp_columns], 1, function(row) mean(row == -99))
filtered_data <- merged_data[missing_data_percentage <= 0.5, ]

# Filter loci based on missing data
missing_data_percentage_loci <- colMeans(filtered_data[snp_columns] == -99)
valid_loci <- names(missing_data_percentage_loci)[missing_data_percentage_loci <= 0.3]
filtered_data <- filtered_data[c('id', 'sex', 'stage', 'genotype', 'rep', 'lat', 'lon', valid_loci)]
str(filtered_data)

# Export filtered data to CSV
#write_csv(filtered_data, './data/filtered_data_for_admixture.csv')



# id row, loci col, 2 id row format (this isn't great because not includ filtering ------------------------------------------

raw_snp_data <- read_csv('./data/Report_DPlatyg23-7805_SNP_2 - Copy corrected.csv', skip = 6)
loci = data.frame(loci = raw_snp_data[,2])
head(loci)
id = raw_snp_data[25:ncol(raw_snp_data)]
names = colnames(id)
#id = data.frame(id = raw_snp_data[25:ncol(raw_snp_data)])
colnames(id) <- names
head(id)
snp_data = cbind(loci, id)
head(snp_data)
nrow(snp_data)
snp_data = snp_data %>% arrange(CloneID)

#test wrangle
(test = snp_data[1:8, 1:8])


# Handling duplicate CloneID entries by keeping only the first two entries when there are more than two
snp_data1 <- snp_data %>%
  group_by(CloneID) %>%
  filter(n() == 2) %>%
  ungroup()
  

library(stringr)
snp_data_long <- snp_data1 %>%
  pivot_longer(-CloneID, names_to = "ID", values_to = "allele") %>%
  group_by(CloneID, ID) %>%
  arrange(ID, CloneID) %>% 
  mutate(allele_row = row_number()) %>%
  ungroup() %>%
  pivot_wider(names_from = CloneID, values_from = allele) %>% 
  mutate(ID = str_replace_all(ID, "\\.", "_")) %>% 
  dplyr::select(.,-c(allele_row)) %>% 
  mutate(across(everything(), ~str_replace_all(as.character(.), "-", "-99")))

head(snp_data_long)

#filter for adults
snp_data_adults<- snp_data_long %>%
  filter(str_detect(ID, "_a_")) %>% data.frame()# Uses 'stringr' to detect '_a_' in the ID column
head(snp_data_adults)
nrow(snp_data_adults)/2  #38
ncol(snp_data_adults) - 1
head(snp_data_adults)
  
write.table(snp_data_adults, file = file.path("./data", "STRUCT_plty_adults.txt"), row.names = FALSE, col.names = F, quote = FALSE)









# COLONY formatting -------------------------------------------------------

genotype_matrix <- as.matrix(tab(data_genind_adult))
ncol(genotype_matrix)
column_names <- colnames(genotype_matrix)
locus_identifiers <- sapply(strsplit(column_names, "/"), `[`, 1)
locus_counts <- table(locus_identifiers)
valid_loci <- names(locus_counts[locus_counts == 2])
valid_columns <- locus_identifiers %in% valid_loci
genotype_matrix <- genotype_matrix[, valid_columns]
ncol(genotype_matrix)
num_columns_to_keep = ncol(genotype_matrix)

##  make marker file ##
(num_loci <- ncol(genotype_matrix)  /2)
marker_names <- paste0("locus", 1:num_loci)

marker_types <- rep(0, num_loci)

# Use default error rates for allelic dropout and other errors (can be modified)
allelic_dropout_rates <- rep(0.01, num_loci)
genotyping_error_rates <- rep(0.01, num_loci)

# Create a matrix for easier handling
marker_info <- rbind(
  marker_names,
  marker_types,
  allelic_dropout_rates,
  genotyping_error_rates
)

# Set the output filename
output_filename <- "./data/marker_info.txt"

write.table(marker_info, file = output_filename, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

marker_info  #this creates a col for each loci. Should be half colony_data_lines -1
ncol(marker_info)

## make offspring file ##

# Extract individual names
individual_names <- indNames(data_genind_adult)

# Prepare for COLONY format: allele pairs per locus
colony_data_lines <- character(nrow(genotype_matrix))
for (i in 1:nrow(genotype_matrix)) {
  individual_data <- c(individual_names[i])  # Start with individual ID
  
  # Loop through loci, taking two columns at a time
  for (j in seq(1, ncol(genotype_matrix), by = 2)) {
    # Extract the two columns representing the dosage for alleles at each locus
    allele1_dosage <- genotype_matrix[i, j]
    allele2_dosage <- genotype_matrix[i, j + 1]
    
    # Determine the COLONY-compatible alleles based on dosage information
    if (is.na(allele1_dosage) || is.na(allele2_dosage)) {
      # Missing data
      alleles <- c(0, 0)
    } else if (allele1_dosage == 2 && allele2_dosage == 0) {
      # Homozygous for first allele, i.e., both alleles are type "C"
      alleles <- c(1, 1)
    } else if (allele1_dosage == 1 && allele2_dosage == 1) {
      # Heterozygous, one "C" and one "G"
      alleles <- c(1, 2)
    } else if (allele1_dosage == 0 && allele2_dosage == 2) {
      # Homozygous for second allele, i.e., both alleles are type "G"
      alleles <- c(2, 2)
    } else {
      # Handle any unexpected cases by assigning missing (failsafe)
      alleles <- c(0, 0)
    }
    
    # Add alleles to the individual's data
    individual_data <- c(individual_data, alleles)
  }
  
  # Combine all elements into a space-separated string
  colony_data_lines[i] <- paste(individual_data, collapse = " ")
}




# Write to a text file for COLONY input
writeLines(colony_data_lines, con = "./data/colony_input.txt")

str(colony_data_lines)
length(unlist(strsplit(colony_data_lines[1], " ")))   #this creates 200 +1 col. 2 col per loci

# split_lines <- strsplit(colony_data_lines, " ")
# colony_df <- do.call(rbind, lapply(split_lines, function(x) x))
# colony_df <- as.data.frame(colony_df, stringsAsFactors = FALSE)



