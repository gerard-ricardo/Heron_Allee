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



# Cervus Platy mapping extraction (working)------------------------------------------------------

#theere were some issues running all larvae at once so had split them up to run

####

#trying to extract from filtered object rather than raw (below) so already filtered.
snp_data_list <- data_gl_filtered@gen
snp_data_matrix <- do.call(cbind, lapply(snp_data_list, as.integer))
data1 <- as.data.frame(snp_data_matrix)
#rownames(data1) <- data_gl_filtered@loc.names  # Locus names as rownames
colnames(data1) <- data_gl_filtered@ind.names  # Individual IDs as colnames
data1$AlleleID = data_gl_filtered@loc.names
data1$CloneID = data_gl_filtered@other$loc.metrics$CloneID
data1$SNP = data_gl_filtered@other$loc.metrics$SNP
head(data1)

####

# SNP 1 Row Mapping Format: "0" = Reference allele homozygote, "1" = SNP allele homozygote, "2"= heterozygote and "-" = double null/null allele homozygote (absence of fragment with SNP in genomic representation)
# data2 <- read.csv(skip = 6, header = T, "./data/Report_DPlatyg23-7805_SNP_mapping_2 - Copy.csv")
# head(data2)
# mean(data2$CallRate) # mean locus callrate is very low. create threshold for above 80%
# data2 <- data2[data2$CallRate > 0.75, ] # subset numerically

####

#wrangling
names(data1) <- gsub("\\.", "_", names(data1))  #use underscores instead

# remove all eggs
# columns_with_e <- grep("_e_", names(data1), value = FALSE)
# data1 <- data1[, -columns_with_e]
# head(data1)


# clean and simplyfy dataframe
data2 <- data1 %>%
  mutate(LocusID = CloneID) %>%
  dplyr::select(., c(LocusID, SNP), starts_with(c("pd")))

# Replacing '-' with NA and converting to character to prevent conversion to factor
data2[data2 == "-"] <- NA
data2 <- mutate_all(data2, as.character)
data2$ref <- substr(data2$SNP, start = nchar(data2$SNP) - 2, stop = nchar(data2$SNP) - 2) # extract ref
data2$var <- substr(data2$SNP, start = nchar(data2$SNP), stop = nchar(data2$SNP)) # extract var
nrow(data2) # no of alleles
length(unique(data2$LocusID)) # distint alleles
head(data2)
duplicated_rows <- data2[duplicated(data2$LocusID), ]
# Remove duplicate rows based on LocusID.Not sure why duplicates but will remove
data2 <- data2 %>% distinct(LocusID, .keep_all = TRUE)
nrow(data2)

# long format prep for two columns
data3a <- data2 %>%
  pivot_longer(
    cols = starts_with(c("pd")),
    names_to = "sample",
    values_to = "genotype"
  ) %>%
  mutate(., rowid = "a") %>%
  data.frame()

data3b <- data2 %>%
  pivot_longer(
    cols = starts_with(c("pd")),
    names_to = "sample",
    values_to = "genotype"
  ) %>%
  mutate(., rowid = "b") %>%
  data.frame()
nrow(data3b)

# Code binary to two letter. NOTE THAT THE SECOND DF USES A DIFFERENT VALUE FOR '2' TO ALLOW FOR heterozygotes
data3a$base <- ifelse(data3a$genotype == 0, data3a$ref,
                      ifelse(data3a$genotype == 1, data3a$var,
                             ifelse(data3a$genotype == 2, data3a$ref,
                                    data3a$genotype
                             )
                      )
)

data3b$base <- ifelse(data3b$genotype == 0, data3b$ref,
                      ifelse(data3b$genotype == 1, data3b$var,
                             ifelse(data3b$genotype == 2, data3b$var,
                                    data3b$genotype
                             )
                      )
)

data4 <- rbind(data3a, data3b) # join back
data4$base <- ifelse(is.na(data4$base), 0, data4$base) # need missing values to be 0 for nalysis.
head(data4)
data5 <- data4 %>%
  select(c(LocusID, sample, base, rowid)) %>%
  dplyr::arrange(., sample, LocusID)
head(data5)
unique(data5$sample)

data5$LocusID <- paste0(data5$LocusID, data5$rowid)
data6 <- data5 %>% dplyr::select(., c(LocusID, sample, base))
data6$base
data_wide <- data6 %>%
  tidyr::pivot_wider(names_from = LocusID, values_from = base) %>%
  data.frame()

rownames(data_wide) <- NULL
head(data_wide)
data_wide$sample
(no_loc <- (ncol(data_wide) - 1) / 2) # no of distict locii
# data_wide$sample <- gsub("\\.", "_", data_wide$sample)

# min typed loci
typed_loci_per_individual <- apply(data_wide, 1, function(x) sum(x != 0))
(min_typed_loci <- min(typed_loci_per_individual))
data_wide <- data_wide[typed_loci_per_individual >= 500, ] # remove all <500 loci
data_wide$sample


# split larvae in ~half (there might be some issues running altogether)
#data_wide1 <- data_wide[grep("_a_|^pd13|^pd14", data_wide$sample), ]
#data_wide2 <- data_wide[grep("_a_|^pd5|^pd9|^pd15", data_wide$sample), ]
#str(data_wide2)
#nrow(data_wide2)
#(ncol(data_wide2)-1)/2  #no of distict locii
#data_wide2$sample
# trimmed test
# quarter_loc = ceiling(no_loc / 4)
# data_wide = data_wide[, 1:(quarter_loc+1)]


#old
# write.csv(data_wide, row.names = FALSE, file = file.path("C:/Users/gerar/OneDrive/1 Work/3 Results/11 Allee effects/3 field experiments/2022_12 Heron/genetics/Cervus", "platy_map_letters_code.csv"))
# write.csv(data_wide1, row.names = FALSE, file = file.path("C:/Users/gerar/OneDrive/1 Work/3 Results/11 Allee effects/3 field experiments/2022_12 Heron/genetics/Cervus", "platy_map_letters_code1.csv"))
# write.csv(data_wide2, row.names = FALSE, file = file.path("C:/Users/gerar/OneDrive/1 Work/3 Results/11 Allee effects/3 field experiments/2022_12 Heron/genetics/Cervus", "platy_map_letters_code2.csv"))
#write.table(data_wide, row.names = FALSE, file = file.path("C:/Users/gerar/OneDrive/1 Work/3 Results/11 Allee effects/3 field experiments/2022_12 Heron/genetics/Cervus", "platy_map_letters_code.txt"))

#new
write.csv(data_wide, row.names = FALSE,
          file = file.path("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects/3_Heron_Platy_ms/Cervus",
                           "platy_map_letters_code2_2.csv"))







### Offspring file
indices_with_l <- grep("l", data_wide$sample)
labels_with_l <- data_wide$sample[indices_with_l]
df1 <- data.frame(offspring = labels_with_l)
dam <- sub("_.*", "", df1$offspring)
#dam_counts <- table(df1$dam)
#known_dam <- rep(names(dam_counts), dam_counts)
# known_dam <- gsub("_(\\d)$", "_0\\1", known_dam)
#df1$known_dam <- known_dam
# df1$Offspring <- gsub("pd(\\d)_", "pd0\\1_", df1$Offspring)
# df1$known_dam <- gsub("pd(\\d)_", "pd0\\1_", df1$known_dam)

indices_with_a <- grep("_a_", data_wide$sample)
labels_with_a <- data_wide$sample[indices_with_a]
first_matching_label <- sapply(dam, function(dam) {
  matching_labels <- grep(dam, labels_with_a, value = TRUE)
  if (length(matching_labels) > 0) {
    return(matching_labels[1])
  } else {
    return(NA)
  }
})
df1$known_dam  = first_matching_label
find_error_dam = which(df1$known_dam == 'pd5_a_1')
df1[find_error_dam, 2] = 'pd13_a_1'
len = length(unique(labels_with_a))
not_known_dam <- labels_with_a[!labels_with_a %in% first_matching_label]
length(not_known_dam)
cands <- rep(labels_with_a, length(df1$known_dam)) %>% sort(.)
cands_df <- matrix(cands, nrow = length(df1$known_dam), ncol = len, byrow = FALSE) %>% data.frame()
colnames(cands_df) <- rep("candidate", len)
offspring_df <- cbind(df1, cands_df)
# for (col in names(cands_df)) {
#   # Add leading zero to single-digit numbers at the end of the strings
#   offspring_df[[col]] <- gsub("_(\\d)$", "_0\\1", offspring_df[[col]])
#
#   # Add leading zero to single-digit numbers following 'pd'
#   offspring_df[[col]] <- gsub("pd(\\d)_", "pd0\\1_", offspring_df[[col]])
# }
offspring_df

#some reason i removed this
#offspring_df <- subset(offspring_df, offspring != "pd15_l_14_10") # remove factor treatment level. Use '%in%' to keep.


# split
#offspring_df1 <- offspring_df[grep("^pd13|^pd14", offspring_df$Offspring), ]
#offspring_df2 <- offspring_df[grep("^pd5|^pd9|^pd15", offspring_df$Offspring), ]

#str(offspring_df)
# write.csv(offspring_df, row.names = FALSE, file = file.path("C:/Users/gerar/OneDrive/1 Work/3 Results/11 Allee effects/3 field experiments/2022_12 Heron/genetics/Cervus", "offspring_platy.csv"))
# write.csv(offspring_df1, row.names = FALSE, file = file.path("C:/Users/gerar/OneDrive/1 Work/3 Results/11 Allee effects/3 field experiments/2022_12 Heron/genetics/Cervus", "offspring_platy1.csv"))
# write.csv(offspring_df2, row.names = FALSE, file = file.path("C:/Users/gerar/OneDrive/1 Work/3 Results/11 Allee effects/3 field experiments/2022_12 Heron/genetics/Cervus", "offspring_platy2.csv"))

write.csv(offspring_df, row.names = FALSE, file = file.path("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects/3_Heron_Platy_ms/Cervus", "offspring_platy.csv"))


# Cervus Acro mapping extraction (working)------------------------------------------------------
# SNP 1 Row Mapping Format: "0" = Reference allele homozygote, "1" = SNP allele homozygote, "2"= heterozygote and "-" = double null/null allele homozygote (absence of fragment with SNP in genomic representation)
data1 <- read.csv(
  skip = 6, header = T,
  file = file.path("C:/Users/gerar/OneDrive/1_Work/3_Results/11 Allee effects/3 field experiments/2022_12 Heron/genetics/Report-DAc23-7804", "Report_DAc23-7804_SNP_mapping_2.csv")
)
head(data1)

data2 <- data1 %>%
  mutate(LocusID = CloneID) %>%
  dplyr::select(., c(LocusID, SNP), starts_with(c("sp", "at")))

data2 <- gl.filter.callrate(data2, method = "ind", threshold = 0.5, v = 3) # filter by ind callrate
missing <- apply(data2[, -c(1, 2)], 2, function(x) sum(x == "-"))
df1 <- data.frame(id = colnames(data2)[3:length(colnames(data2))], callrate_ind = 1 - (missing / nrow(data2)))
df2 <- df1[df1$callrate_ind > 0.45, ] # subset numerically

# 29337
data2 <- gl.filter.callrate(data2, method = "loc", threshold = 0.7, v = 3) # filter by loci callrate
zero_counts <- apply(data2[, -c(1, 2)], 1, function(x) sum(x == 0))


data2 <- gl.filter.callrate(data2, method = "loc", threshold = 0.7, v = 3) # filter by loci callrate
# 552
data2 <- gl.filter.callrate(data2, method = "ind", threshold = 0.7, v = 3) # filter by ind callrate
# 552
data2 <- gl.filter.reproducibility(data2, t = 0.7, v = 3) # filter out loci with limited reproducibility
# 552
data2 <- gl.filter.monomorphs(data2, v = 3) # remove monomorphic loci (loci with 1 fixed allele across the entire dataset (no differences) )
# 249
data2 <- gl.filter.hwe(data2, alpha_val = 0.05, subset = "each", multi_comp_method = "bonferroni", v = 3) # filter out loci that depart from H-W proportions
# 161
data2 <- gl.filter.secondaries(data2, method = "random", verbose = 3) # remove loci fragment that shared SNPs. Only keep 1
# 151
list.match <- data2$loc.names[which(data2$other$loc.metrics$OneRatioSnp > 0.05 & data2$other$loc.metrics$OneRatioSnp < 0.95 & data2$other$loc.metrics$OneRatioRef < 0.95 & data2$other$loc.metrics$OneRatioRef > 0.05 & data2$other$loc.metrics$coverage > 10)] # remove loci based on minor allele frequency and low data coverage
data2 <- data2[, match(list.match, data2$loc.names)] # keep only loci in the list above
# 141 loci left



# write.csv(data2, file = file.path("C:/Users/gerar/OneDrive/1 Work/3 Results/11 Allee effects/3 field experiments/2022_12 Heron/genetics/Cervus", "mapping_cleaned.csv"))


# Replacing '-' with NA and converting to character to prevent conversion to factor
data2[data2 == "-"] <- NA
data2 <- mutate_all(data2, as.character)
data2$ref <- substr(data2$SNP, start = nchar(data2$SNP) - 2, stop = nchar(data2$SNP) - 2)
data2$var <- substr(data2$SNP, start = nchar(data2$SNP), stop = nchar(data2$SNP))
nrow(data2)
length(unique(data2$LocusID))
head(data2)
duplicated_rows <- data2[duplicated(data2$LocusID), ]
# Remove duplicate rows based on LocusID.Not sure why duplicates but will remove
data2 <- data2 %>%
  distinct(LocusID, .keep_all = TRUE)
nrow(data2)


data3a <- data2 %>%
  pivot_longer(
    cols = starts_with(c("sp", "at")),
    names_to = "sample",
    values_to = "genotype"
  ) %>%
  mutate(., rowid = "a") %>%
  data.frame()
data3b <- data2 %>%
  pivot_longer(
    cols = starts_with(c("sp", "at")),
    names_to = "sample",
    values_to = "genotype"
  ) %>%
  mutate(., rowid = "b") %>%
  data.frame()
nrow(data3b)
# data3a <- data3a %>%
#   mutate(LocusID = paste0(LocusID, "a"))
# data3b <- data3a %>%
#   mutate(LocusID = paste0(LocusID, "b"))
# head(data3a)
# NOTE THAT THE SECOND DF USES A DIFFERENT VALUE FOR '2' TO ALLOW FOR ALLELES
data3a$base <- ifelse(data3a$genotype == 0, data3a$ref,
                      ifelse(data3a$genotype == 1, data3a$var,
                             ifelse(data3a$genotype == 2, data3a$ref,
                                    data3a$genotype
                             )
                      )
)
data3b$base <- ifelse(data3b$genotype == 0, data3b$ref,
                      ifelse(data3b$genotype == 1, data3b$var,
                             ifelse(data3b$genotype == 2, data3b$var,
                                    data3b$genotype
                             )
                      )
)

data4 <- rbind(data3a, data3b)
data4$base <- ifelse(is.na(data4$base), 0, data4$base)
head(data4)
data5 <- data5 %>%
  select(c(LocusID, sample, base, rowid)) %>%
  dplyr::arrange(., sample, LocusID)
head(data5)


# # Pivot to wide format
# data5 <- data4 %>%
#   group_by(LocusID, sample)  %>% data.frame()
# head(data5)

data5$LocusID <- paste0(data5$LocusID, data5$rowid)
data6 <- data5 %>% dplyr::select(., c(LocusID, sample, base))
data_wide <- data6 %>%
  tidyr::pivot_wider(names_from = LocusID, values_from = base) %>%
  data.frame()
rownames(data_wide) <- NULL
(ncol(data_wide) - 1) / 2

# make to remove excess column in output
# write.csv(data_wide,  row.names = FALSE, file = file.path("C:/Users/gerar/OneDrive/1 Work/3 Results/11 Allee effects/3 field experiments/2022_12 Heron/genetics/Cervus", "mapping_letters_coding.csv"))



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





