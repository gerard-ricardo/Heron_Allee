





input_file <- "C:\\Users\\gerar\\OneDrive\\1_Work\\4_Writing\\1_Allee_effects\\3_Heron_Platy_ms\\Cervus\\pd_afa_out2.txt"

lines <- readLines(input_file)

start_line <- 11
end_line <- 797
extracted_lines <- lines[start_line:end_line]

split_lines <- strsplit(extracted_lines, " +")

column_names <- split_lines[[1]]

data1 <- do.call(rbind, lapply(split_lines[-1], function(x) {
  length(x) <- length(column_names) # Pad with NA to match the number of columns
  return(x)
}))

data1 <- as.data.frame(data1, stringsAsFactors = FALSE)
colnames(data1) <- column_names
head(data1)
str(data1)
data1$HObs   <- as.numeric(as.character(data1$HObs))
data1$HExp   <- as.numeric(as.character(data1$HExp))
data1$`F(Null)` <- as.numeric(data1$`F(Null)`)
nrow(data1)
length(which(data1$`F(Null)` > 0.05))

threshold <- 0.3   #
data1 <- subset(data1, `F(Null)` <= threshold)
nrow(data1)







check = 100156704  #loci to check from cervus output
offspring_id <- 'pd13.l.14.11'
mother_id <- 'pd13.a.1'



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
snp_data$AlleleSequence[1] == snp_data$AlleleSequence[2]
snp_data2 %>%  dplyr::select(., CloneID, AlleleSequenceRef , AlleleSequenceSnp,pd13.l.14.11, pd13.a.1) %>% filter(., CloneID == check )




ind_metrics <- data_genind$other$ind.metrics


ind_metrics <- data_genind$other$ind.metrics
loc_metrics = data_genind@tab %>%  data.frame()
data_gl_filtered@
  
  ind_metrics$genotype <- sub("\\..*", "", ind_metrics$id)

mothers <- ind_metrics %>% filter(stage == "adults")
offspring <- ind_metrics %>% filter(stage == "larvae")

offspring$mother_id <- sapply(offspring$genotype, function(genotype) {
  match_indices <- grep(genotype, mothers$id)
  if(length(match_indices) > 0) {
    return(mothers$id[match_indices[1]])  # Return the first match
  } else {
    return(NA)  # Return NA if no match is found
  }
})

pairs <- data.frame(offspring = offspring$id, mother = offspring$mother_id)

results <- data.frame(IndividualID = character(), SNP = character(), MendelianInconsistency = logical())

offspring_id <- pairs$offspring[1]
mother_id <- pairs$mother[1]

loc_metrics
comp2 <- loc_metrics %>%
  dplyr::filter(row.names(.) %in% c(offspring_id, mother_id))
comp2[1:2, 1:5]
parent <- comp2[1, ] # Assuming first row is parent
offspring <- comp2[2, ] # Assuming second row is offspring

complete_cases <- complete.cases(parent, offspring)
parent_filtered <- parent[complete_cases]
offspring_filtered <- offspring[complete_cases]

is_mismatch <- function(parent_allele, offspring_allele) {
  return((parent_allele == 0 & offspring_allele == 2) | (parent_allele == 2 & offspring_allele == 0))
}

mismatches <- mapply(is_mismatch, parent_filtered, offspring_filtered)

mismatched_loci <- names(mismatches)[mismatches]

parent_mismatches <- parent_filtered[mismatches]
offspring_mismatches <- offspring_filtered[mismatches]

mismatch_df <- data.frame(
  Locus = mismatched_loci,
  Parent_Genotype = parent_mismatches,
  Offspring_Genotype = offspring_mismatches
)

mismatch_df

complete_cases <- complete.cases(parent, offspring)
parent_filtered <- parent[complete_cases]
offspring_filtered <- offspring[complete_cases]

is_mismatch <- function(parent_allele, offspring_allele) {
  return((parent_allele == 0 & offspring_allele == 2) | (parent_allele == 2 & offspring_allele == 0))
}
mismatches <- mapply(is_mismatch, parent, offspring)
mismatched_loci <- names(mismatches)[mismatches]


for (i in 1:nrow(pairs)) {
  
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






data1 <- data_genind@tab  # This assumes SNP data is stored in the 'tab' slot of the genind object
head(data1)

data1 = t(data1) %>% data.frame()
colnames(data1) <- gsub("\\.", "_", colnames(data1))
data1$rownames <- rownames(data1)
data1 <- data1 %>% mutate(CloneID = sub("-.*", "", rownames), allele = sub(".*\\.(.)$", "\\1", rownames)) %>% 
  select(CloneID, allele, everything()) %>% select(-rownames) %>% arrange(., CloneID)
data1 %>%  dplyr::select(., CloneID, pd13_l_14_11, pd13_a_1) %>% filter(., CloneID == check )  #


data1_long = data1 %>% tidyr::pivot_longer(-c(CloneID, allele) ,  names_to = "id" ,values_to = "counts") %>% arrange(., CloneID, id) %>% data.frame()  #keep vec.x, add all other columns to factors , add all their values to meas)


nrow(data1_long)  #91413
row_counts = data1_long %>%
  group_by(CloneID, id) %>%
  summarise(row_count = n()) %>%
  ungroup()
min(row_counts$row_count)

data1_long = data1_long %>%
  group_by(CloneID, id) %>%
  filter(n() == 2) %>%
  ungroup() %>%
  data.frame()
nrow(data1_long)  #81144
any(data1_long == "", na.rm = TRUE) # This will return TRUE if there are any empty strings in the dataframe


repeat_with_separator <- function(allele, counts, sep = "_") {
  if (is.na(allele) | is.na(counts)) {
    return("NA")
  } else {
    return(paste(rep(allele, counts), collapse = sep))
  }
}

data1_long <- data1_long %>%
  mutate(allele_multiplied = mapply(repeat_with_separator, allele, counts, sep = ",")) # you can change sep to "_" if you want underscores
any(data1_long == "", na.rm = TRUE) # This will return TRUE if there are any empty strings in the dataframe

print(data1_long)

summary_data <- data1_long %>%
  group_by(CloneID, id) %>%
  summarise(
    id = first(id),
    base = paste(allele_multiplied, collapse = ",")
  ) %>%
  arrange(id) %>%
  data.frame()

summary_data$base <- gsub("^,+|,+$", "", summary_data$base)


summary_data <- summary_data %>%
  separate(base, into = c("a", "b"), sep = ",", extra = "drop", fill = "right")

data1_long = summary_data %>% tidyr::pivot_longer(-c(CloneID, id) ,  names_to = "base" ,values_to = "code") %>% arrange(., CloneID, id) %>% data.frame()  #keep vec.x, add all other columns to factors , add all their values to meas)
data1_long$LocusID = paste0(data1_long$CloneID, data1_long$base)
data1_long = data1_long %>% select(-c(CloneID , base))
data1_long <- data1_long %>%
  mutate(code = ifelse(code == "NA", "*", code))
data_wide <- data1_long %>% tidyr::pivot_wider(names_from = LocusID, values_from = code, names_prefix = "X") %>% 
  data.frame()#year goes to columns, their areas go as the values, area is the prefix

  


data1 %>%
  dplyr::select(CloneID, SNP, pd13_l_14_11, pd13_a_1) %>%  # Select specified columns
  filter(CloneID == check)  # Filter rows where CloneID equals the value of 'check'





is_mismatch <- function(parent_allele, offspring_allele) {
  return((parent_allele == 0 & offspring_allele == 2) | (parent_allele == 2 & offspring_allele == 0))
}
mismatches <- mapply(is_mismatch, test_df[2], test_df[3])
test_df[mismatches, ]














write.csv(data_wide, row.names = FALSE,
          file = file.path("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects/3_Heron_Platy_ms/Cervus",
                           "platy_map_letters_code2_2.csv"))


indices_with_l <- grep("l", data_wide$sample)
labels_with_l <- data_wide$sample[indices_with_l]
df1 <- data.frame(offspring = labels_with_l)
dam <- sub("_.*", "", df1$offspring)

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
offspring_df





write.csv(offspring_df, row.names = FALSE, file = file.path("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects/3_Heron_Platy_ms/Cervus", "offspring_platy.csv"))


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

data2 <- gl.filter.callrate(data2, method = "loc", threshold = 0.7, v = 3) # filter by loci callrate
zero_counts <- apply(data2[, -c(1, 2)], 1, function(x) sum(x == 0))


data2 <- gl.filter.callrate(data2, method = "loc", threshold = 0.7, v = 3) # filter by loci callrate
data2 <- gl.filter.callrate(data2, method = "ind", threshold = 0.7, v = 3) # filter by ind callrate
data2 <- gl.filter.reproducibility(data2, t = 0.7, v = 3) # filter out loci with limited reproducibility
data2 <- gl.filter.monomorphs(data2, v = 3) # remove monomorphic loci (loci with 1 fixed allele across the entire dataset (no differences) )
data2 <- gl.filter.hwe(data2, alpha_val = 0.05, subset = "each", multi_comp_method = "bonferroni", v = 3) # filter out loci that depart from H-W proportions
data2 <- gl.filter.secondaries(data2, method = "random", verbose = 3) # remove loci fragment that shared SNPs. Only keep 1
list.match <- data2$loc.names[which(data2$other$loc.metrics$OneRatioSnp > 0.05 & data2$other$loc.metrics$OneRatioSnp < 0.95 & data2$other$loc.metrics$OneRatioRef < 0.95 & data2$other$loc.metrics$OneRatioRef > 0.05 & data2$other$loc.metrics$coverage > 10)] # remove loci based on minor allele frequency and low data coverage
data2 <- data2[, match(list.match, data2$loc.names)] # keep only loci in the list above





data2[data2 == "-"] <- NA
data2 <- mutate_all(data2, as.character)
data2$ref <- substr(data2$SNP, start = nchar(data2$SNP) - 2, stop = nchar(data2$SNP) - 2)
data2$var <- substr(data2$SNP, start = nchar(data2$SNP), stop = nchar(data2$SNP))
nrow(data2)
length(unique(data2$LocusID))
head(data2)
duplicated_rows <- data2[duplicated(data2$LocusID), ]
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



data5$LocusID <- paste0(data5$LocusID, data5$rowid)
data6 <- data5 %>% dplyr::select(., c(LocusID, sample, base))
data_wide <- data6 %>%
  tidyr::pivot_wider(names_from = LocusID, values_from = base) %>%
  data.frame()
rownames(data_wide) <- NULL
(ncol(data_wide) - 1) / 2


