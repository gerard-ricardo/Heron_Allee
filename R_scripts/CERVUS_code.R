#CERVUS and parent-offspring mismatches


#there is an error between the original files and the conversion from data_gl_filtered. Tryng to work on genind.
# - Convert genind to bases
#seems to work but has empty values
#Using CERVUS to check mismatches

#Ive tunred of CERVUS code until I resolve this

# might want to check that there are two col for each locus and rerun

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



# -------------------------------------------------------------------------




#inputs
check = 100156704  #loci to check from cervus output
offspring_id <- 'pd13.l.14.11'
mother_id <- 'pd13.a.1'


# compare SNP_2 file with SNP_mapping  --------------------------------

# first check raw CSV shows mismatch
raw_snp_data <- read_csv('./data/Report_DPlatyg23-7805_SNP_2 - Copy corrected.csv', skip = 6)
raw_snp_data2 <- read_csv('./data/Report_DPlatyg23-7805_SNP_mapping_2 - Copy.csv', skip = 6)

colnames(raw_snp_data)
loci = select(raw_snp_data, CloneID, AlleleSequence, SNP)
colnames(raw_snp_data2)
loci2 = select(raw_snp_data2, CloneID, AlleleSequenceRef, AlleleSequenceSnp)

head(loci)
id = raw_snp_data[25:ncol(raw_snp_data)]
id2 = raw_snp_data2[27:ncol(raw_snp_data2)]
names = colnames(id)
snp_data = cbind(loci, id)
snp_data2 = cbind(loci2, id2)
snp_data = snp_data %>% arrange(CloneID)
snp_data %>%  dplyr::select(., CloneID, AlleleSequence , SNP,pd13.l.14.11, pd13.a.1) %>% filter(., CloneID == check )
#1/1 indicates hetero
snp_data$AlleleSequence[1] == snp_data$AlleleSequence[2]
snp_data2 %>%  dplyr::select(., CloneID, AlleleSequenceRef , AlleleSequenceSnp,pd13.l.14.11, pd13.a.1) %>% filter(., CloneID == check )
#2 indicates hetero

#conclusion, the SNP_2 file matches the SNP_mapping_2 file. Both show adult  hetero, and the offspring homo ref


# filter for parent-offspring mismatch (not working yet) ------------------------------------

ind_metrics <- data_genind$other$ind.metrics

#data1 = data_genind$other$ind.metrics
#data1$id2 <- gsub("\\.", "_", data1$id)  #use underscores instead
#data_genind$other$ind.metrics$genotype <- sub("\\..*", "", data1$id)

# Extract the meta and loc metrics
ind_metrics <- data_genind$other$ind.metrics
loc_metrics = data_genind@tab %>%  data.frame()
data_gl_filtered@
  
  # Ensure the genotype column is present
  ind_metrics$genotype <- sub("\\..*", "", ind_metrics$id)

# Identify mother-offspring pairs
mothers <- ind_metrics %>% filter(stage == "adults")
offspring <- ind_metrics %>% filter(stage == "larvae")

# Find the first matching mother id for each offspring genotype
offspring$mother_id <- sapply(offspring$genotype, function(genotype) {
  match_indices <- grep(genotype, mothers$id)
  if(length(match_indices) > 0) {
    return(mothers$id[match_indices[1]])  # Return the first match
  } else {
    return(NA)  # Return NA if no match is found
  }
})

#create data frame of larvae and their mums
pairs <- data.frame(offspring = offspring$id, mother = offspring$mother_id)

#initialise emtpy df
results <- data.frame(IndividualID = character(), SNP = character(), MendelianInconsistency = logical())

##compare single parent-offspring loci
offspring_id <- pairs$offspring[1]
mother_id <- pairs$mother[1]

#compare a locus - note that 0 (homo for reference) and 2 (homo for alternative) pairs indicates mismatch.
#Also, issues comparing loci with missing data (secondary issue)
loc_metrics
comp2 <- loc_metrics %>%
  dplyr::filter(row.names(.) %in% c(offspring_id, mother_id))
comp2[1:2, 1:5]
# Extract the parent and offspring rows
parent <- comp2[1, ] # Assuming first row is parent
offspring <- comp2[2, ] # Assuming second row is offspring

# Remove columns with NA values in either parent or offspring rows
complete_cases <- complete.cases(parent, offspring)
parent_filtered <- parent[complete_cases]
offspring_filtered <- offspring[complete_cases]

# Define a function to check for mismatches where one value is 0 and the other is 2
is_mismatch <- function(parent_allele, offspring_allele) {
  return((parent_allele == 0 & offspring_allele == 2) | (parent_allele == 2 & offspring_allele == 0))
}

# Apply the mismatch function to each locus
mismatches <- mapply(is_mismatch, parent_filtered, offspring_filtered)

# Display the mismatches
mismatched_loci <- names(mismatches)[mismatches]

# Display the parent and offspring genotypes at mismatched loci
parent_mismatches <- parent_filtered[mismatches]
offspring_mismatches <- offspring_filtered[mismatches]

# Create a data frame to show the mismatches
mismatch_df <- data.frame(
  Locus = mismatched_loci,
  Parent_Genotype = parent_mismatches,
  Offspring_Genotype = offspring_mismatches
)

mismatch_df

# Remove columns with NA values in either parent or offspring rows
complete_cases <- complete.cases(parent, offspring)
parent_filtered <- parent[complete_cases]
offspring_filtered <- offspring[complete_cases]

is_mismatch <- function(parent_allele, offspring_allele) {
  return((parent_allele == 0 & offspring_allele == 2) | (parent_allele == 2 & offspring_allele == 0))
}
mismatches <- mapply(is_mismatch, parent, offspring)
mismatched_loci <- names(mismatches)[mismatches]


for (i in 1:nrow(pairs)) {
  
  # offspring_id <- pairs$offspring[i]
  # mother_id <- pairs$mother[i]
  if (is.na(mother_id)) next  # Skip if no matching mother

  for (snp in locNames(data_genind)) {
    mother_genotype <- paste(data_genind@tab[mother_id, "85078316-7-G/A"], collapse = "")
    
    ff = data_genind@tab
    
    mother_genotype <- paste(data_genind@tab[mother_id, snp], collapse = "")
    child_genotype <- paste(data_genind@tab[offspring_id, snp], collapse = "")
    
    if (any(is.na(c(mother_genotype, child_genotype)))) {
      inconsistency <- NA  # Missing data
    } else {
      mother_alleles <- unlist(strsplit(mother_genotype, ""))
      child_alleles <- unlist(strsplit(child_genotype, ""))
      possible_alleles <- unique(mother_alleles)
      inconsistency <- !(all(child_alleles %in% possible_alleles))
    }
    
    results <- rbind(results, data.frame(IndividualID = offspring_id, SNP = snp, MendelianInconsistency = inconsistency))
  }
}

results
data_genind


# Cervus Platy mapping extraction  ------------------------------------------------------

#using inputs above to check mismatches

### extracting from the genlight dart file is causing errors
## SNP 1 Row Mapping Format: "0" = Reference allele homozygote, "1" = SNP allele homozygote, "2"= heterozygote 
#and "-" = double null/null allele homozygote (absence of fragment with SNP in genomic representation)
# snp_data_list <- data_gl_filtered@gen
# as.integer(snp_data_list[[1]])
# snp_data_matrix <- do.call(cbind, lapply(snp_data_list, as.integer))
# data1 <- as.data.frame(snp_data_matrix)
# #rownames(data1) <- data_gl_filtered@loc.names  # Locus names as rownames
# colnames(data1) <- data_gl_filtered@ind.names  # Individual IDs as colnames
# data1$AlleleID = data_gl_filtered@loc.names
# data1$CloneID = data_gl_filtered@other$loc.metrics$CloneID
# data1$SNP = data_gl_filtered@other$loc.metrics$SNP
# head(data1)
# names(data1) <- gsub("\\.", "_", names(data1))  #use underscores instead
# 
# #parent-offspring mismatch test
# data1 %>%  dplyr::select(., CloneID, SNP , pd13_l_14_11, pd13_a_1) %>% filter(., CloneID == check )  #
#CloneID    SNP       pd13_l_14_11     pd13_a_1
#100156704  25:T>C     0               1
#indicates adult is homo for var (CC), and offspring is homo for ref (TT). So there is a mismatch between here and original files


#try genind
# Extracting SNP data list from genind object
data1 <- data_genind@tab  # This assumes SNP data is stored in the 'tab' slot of the genind object
head(data1)

data1 = t(data1) %>% data.frame()
#data1$SNP = data_genind@other$loc.metrics$SNP
colnames(data1) <- gsub("\\.", "_", colnames(data1))
data1$rownames <- rownames(data1)
data1 <- data1 %>% mutate(CloneID = sub("-.*", "", rownames), allele = sub(".*\\.(.)$", "\\1", rownames)) %>% 
  select(CloneID, allele, everything()) %>% select(-rownames) %>% arrange(., CloneID)
data1 %>%  dplyr::select(., CloneID, pd13_l_14_11, pd13_a_1) %>% filter(., CloneID == check )  #

#this is correct, so there must be an issue with extracting from the genlight

#data2 = data1[8:16, 1:5]
data1_long = data1 %>% tidyr::pivot_longer(-c(CloneID, allele) ,  names_to = "id" ,values_to = "counts") %>% arrange(., CloneID, id) %>% data.frame()  #keep vec.x, add all other columns to factors , add all their values to meas)


nrow(data1_long)  #91413
row_counts = data1_long %>%
  group_by(CloneID, id) %>%
  summarise(row_count = n()) %>%
  ungroup()
min(row_counts$row_count)

# Filter out groups that do not have exactly 2 rows
data1_long = data1_long %>%
  group_by(CloneID, id) %>%
  filter(n() == 2) %>%
  ungroup() %>%
  data.frame()
nrow(data1_long)  #81144
any(data1_long == "", na.rm = TRUE) # This will return TRUE if there are any empty strings in the dataframe

#data1_long = data1_long[129:140,]

# Function to repeat allele with a separator
repeat_with_separator <- function(allele, counts, sep = "_") {
  if (is.na(allele) | is.na(counts)) {
    return("NA")
  } else {
    return(paste(rep(allele, counts), collapse = sep))
  }
}

# Apply the function to the dataframe
data1_long <- data1_long %>%
  mutate(allele_multiplied = mapply(repeat_with_separator, allele, counts, sep = ",")) # you can change sep to "_" if you want underscores
any(data1_long == "", na.rm = TRUE) # This will return TRUE if there are any empty strings in the dataframe

# Display the intermediate result
print(data1_long)

# Group and summarise the data
summary_data <- data1_long %>%
  group_by(CloneID, id) %>%
  summarise(
    id = first(id),
    base = paste(allele_multiplied, collapse = ",")
  ) %>%
  arrange(id) %>%
  data.frame()

# Remove leading or trailing commas in summary_data$base
summary_data$base <- gsub("^,+|,+$", "", summary_data$base)


# Separate the 'base' column into 'base1' and 'base2'
summary_data <- summary_data %>%
  separate(base, into = c("a", "b"), sep = ",", extra = "drop", fill = "right")

data1_long = summary_data %>% tidyr::pivot_longer(-c(CloneID, id) ,  names_to = "base" ,values_to = "code") %>% arrange(., CloneID, id) %>% data.frame()  #keep vec.x, add all other columns to factors , add all their values to meas)
data1_long$LocusID = paste0(data1_long$CloneID, data1_long$base)
data1_long = data1_long %>% select(-c(CloneID , base))
data1_long <- data1_long %>%
  mutate(code = ifelse(code == "NA", "*", code))
data_wide <- data1_long %>% tidyr::pivot_wider(names_from = LocusID, values_from = code, names_prefix = "X") %>% 
  data.frame()#year goes to columns, their areas go as the values, area is the prefix

  
#finished


###############
# Parent-offspring mismatch test
data1 %>%
  dplyr::select(CloneID, SNP, pd13_l_14_11, pd13_a_1) %>%  # Select specified columns
  filter(CloneID == check)  # Filter rows where CloneID equals the value of 'check'





is_mismatch <- function(parent_allele, offspring_allele) {
  return((parent_allele == 0 & offspring_allele == 2) | (parent_allele == 2 & offspring_allele == 0))
}
mismatches <- mapply(is_mismatch, test_df[2], test_df[3])
test_df[mismatches, ]

####

# SNP 1 Row Mapping Format: "0" = Reference allele homozygote, "1" = SNP allele homozygote, "2"= heterozygote and "-" = double null/null allele homozygote (absence of fragment with SNP in genomic representation)
# data2 <- read.csv(skip = 6, header = T, "./data/Report_DPlatyg23-7805_SNP_mapping_2 - Copy.csv")
# head(data2)
# mean(data2$CallRate) # mean locus callrate is very low. create threshold for above 80%
# data2 <- data2[data2$CallRate > 0.75, ] # subset numerically

####

#wrangling

# remove all eggs
# columns_with_e <- grep("_e_", names(data1), value = FALSE)
# data1 <- data1[, -columns_with_e]
# head(data1)


# clean and simplyfy dataframe
# data2 <- data1 %>%
#   mutate(LocusID = CloneID) %>%
#   dplyr::select(., c(LocusID, SNP), starts_with(c("pd")))
# 
# # Replacing '-' with NA and converting to character to prevent conversion to factor
# data2[data2 == "-"] <- NA
# data2 <- mutate_all(data2, as.character)
# data2$ref <- substr(data2$SNP, start = nchar(data2$SNP) - 2, stop = nchar(data2$SNP) - 2) # extract ref
# data2$var <- substr(data2$SNP, start = nchar(data2$SNP), stop = nchar(data2$SNP)) # extract var
# nrow(data2) # no of alleles
# length(unique(data2$LocusID)) # distint alleles
# head(data2)
# duplicated_rows <- data2[duplicated(data2$LocusID), ]
# # Remove duplicate rows based on LocusID.Not sure why duplicates but will remove
# data2 <- data2 %>% distinct(LocusID, .keep_all = TRUE)
# nrow(data2)
# head(data2)
# data2 %>%  dplyr::select(., LocusID, SNP , pd13_l_14_11, pd13_a_1) %>% filter(., LocusID == check )  #this is same as above
# 
# 
# # long format prep for two columns
# data3a <- data2 %>%
#   pivot_longer(
#     cols = starts_with(c("pd")),
#     names_to = "sample",
#     values_to = "genotype"
#   ) %>%
#   mutate(., rowid = "a") %>%
#   data.frame()
# 
# data3a %>% filter(., sample %in% c('pd13_l_14_11', 'pd13_a_1') & LocusID == check)  #same as above
# 
# data3b <- data2 %>%
#   pivot_longer(
#     cols = starts_with(c("pd")),
#     names_to = "sample",
#     values_to = "genotype"
#   ) %>%
#   mutate(., rowid = "b") %>%
#   data.frame()
# nrow(data3b)
# 
# #description: SNP 1 Row Mapping Format: "0" = Reference allele homozygote (ref/ref), "1" = SNP allele homozygote (ref/var), 
# #"2"= heterozygote (var/var) and "-" = double null/null allele homozygote (absence of fragment with SNP in genomic representation)
# data3a$base <- ifelse(data3a$genotype == 0, data3a$ref,
#                       ifelse(data3a$genotype == 1, data3a$var,
#                              ifelse(data3a$genotype == 2, data3a$ref,
#                                     data3a$genotype
#                              )
#                       )
# )
# data3a
# 
# data3b$base <- ifelse(data3b$genotype == 0, data3b$ref,
#                       ifelse(data3b$genotype == 1, data3b$var,
#                              ifelse(data3b$genotype == 2, data3b$var,
#                                     data3b$genotype
#                              )
#                       )
# )
# 
# data4 <- rbind(data3a, data3b) # join back
# data4$base <- ifelse(is.na(data4$base), '*', data4$base) # need missing values to be 0 for nalysis.
# head(data4)
# data4 %>% filter(., sample %in% c('pd13_l_14_11', 'pd13_a_1') & LocusID == check)  #issue here
# 
# 
# data5 <- data4 %>%
#   select(c(LocusID, sample, base, rowid)) %>%
#   dplyr::arrange(., sample, LocusID)
# head(data5)
# unique(data5$sample)
# 
# data5$LocusID <- paste0(data5$LocusID, data5$rowid)
# data6 <- data5 %>% dplyr::select(., c(LocusID, sample, base))
# data6$base
# data_wide <- data6 %>%
#   tidyr::pivot_wider(names_from = LocusID, values_from = base) %>%
#   data.frame()
# 
# rownames(data_wide) <- NULL
# head(data_wide)
# data_wide$sample
# (no_loc <- (ncol(data_wide) - 1) / 2) # no of distict locii
# # data_wide$sample <- gsub("\\.", "_", data_wide$sample)
# 
# # min typed loci
# typed_loci_per_individual <- apply(data_wide, 1, function(x) sum(x != 0))
# (min_typed_loci <- min(typed_loci_per_individual))
# data_wide <- data_wide[typed_loci_per_individual >= 500, ] # remove all <500 loci
# data_wide$sample


####checking mismatch from cervus output
# which(data_wide$sample == 'pd13_l_14_11')
# which(data_wide$sample == 'pd13_a_1')
# data_wide$'X100156704a'[9]
# data_wide$'X100156704b'[9]
# data_wide$'X100156704a'[7]
# data_wide$'X100156704b'[7]
# 
# 
# head(data4)
# data4$sample
# data4$LocusID <- paste0(data4$LocusID, data4$rowid)
# 
# data4_pd13_l_14_11 = data4[grep("pd13_l_14_11", data4$sample),]
# data4_pd13_a_1 = data4[grep("pd13_a_1", data4$sample),]
# 
# sing_join = left_join(data4_pd13_a_1, data4_pd13_l_14_11,  by = 'LocusID')
# nrow(sing_join)
# sing_join[which(sing_join$SNP.x == sing_join$SNP.y),]
# 
# #checking mismatch from cervus output
# which(sing_join$sample == 'pd13_l_14_11')
# which(sing_join$sample == 'pd13_a_1')
# 
# sing_join %>% filter(LocusID == "85069500a" & sample.x == "pd13_a_1")
# sing_join %>% filter(LocusID == "85069500b" & sample.x == "pd13_a_1")
# 
# #ok seems to be mismatch here from. Parent is AA so child cant be GG

#######

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

