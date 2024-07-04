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



# id row, loci col, 2 row format ------------------------------------------

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
  filter(str_detect(ID, "_a_")) # Uses 'stringr' to detect '_a_' in the ID column
head(snp_data_adults)
nrow(snp_data_adults)/2  #38
ncol(snp_data_adults) - 1
  
write.table(snp_data_adults, file = file.path("./data", "STRUCT_plty_adults.txt"), row.names = FALSE, col.names = F, quote = FALSE)


