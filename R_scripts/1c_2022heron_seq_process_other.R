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




# input cervus outputs ----------------------------------------------------

#pd_afa_out2.txt
# Define the path to the input file
input_file <- "C:\\Users\\gerar\\OneDrive\\1_Work\\4_Writing\\1_Allee_effects\\3_Heron_Platy_ms\\Cervus\\pd_afa_out2.txt"

# Read the entire file into R
lines <- readLines(input_file)

# Extract lines 11 to 797
start_line <- 11
end_line <- 797
extracted_lines <- lines[start_line:end_line]

# Split each line by spaces
split_lines <- strsplit(extracted_lines, " +")

# Extract the column names from the first element
column_names <- split_lines[[1]]

# Convert the remaining split lines into a dataframe
data1 <- do.call(rbind, lapply(split_lines[-1], function(x) {
  length(x) <- length(column_names) # Pad with NA to match the number of columns
  return(x)
}))

# Convert to dataframe and assign column names
data1 <- as.data.frame(data1, stringsAsFactors = FALSE)
colnames(data1) <- column_names
head(data1)
str(data1)
data1$HObs   <- as.numeric(as.character(data1$HObs))
data1$HExp   <- as.numeric(as.character(data1$HExp))
data1$`F(Null)` <- as.numeric(data1$`F(Null)`)
nrow(data1)
length(which(data1$`F(Null)` > 0.05))

## filter out high null alleles - CAUTION this can reduce homo excess
# Define a threshold for filtering out null alleles
threshold <- 0.3   #
data1 <- subset(data1, `F(Null)` <= threshold)
nrow(data1)





