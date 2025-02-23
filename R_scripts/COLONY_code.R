data_genind_adult


genotype_data_means <- data_genind_adult@other$loc.metrics %>%
  dplyr::select(CallRate:coverage) %>%  # Select the desired range of columns
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))  # Use anonymous function to pass na.rm

genotype_data_means_sub1 <- data_genind_adult_subset1@other$loc.metrics %>%
  dplyr::select(CallRate:coverage) %>%  # Select the desired range of columns
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))  # Use anonymous function to pass na.rm


genotyping_error_rate <- 1 - genotype_data$CallRate
med_genotyping_error_rate <- median(genotyping_error_rate, na.rm = TRUE)





colony_data <- function(genind_obj) {
  
  input_name <- deparse(substitute(genind_obj))
  
  
  
  hist(genind_obj@other$loc.metrics$AvgPIC )
  
  genotype_matrix <- as.matrix(tab(genind_obj))
  ncol(genotype_matrix)
  column_names <- colnames(genotype_matrix)
  locus_identifiers <- sapply(strsplit(column_names, "/"), `[`, 1)
  locus_counts <- table(locus_identifiers)
  valid_loci <- names(locus_counts[locus_counts == 2])
  valid_columns <- locus_identifiers %in% valid_loci
  genotype_matrix <- genotype_matrix[, valid_columns]
  ncol(genotype_matrix)
  num_columns_to_keep = ncol(genotype_matrix)
  
  
  
  (num_loci <- ncol(genotype_matrix)  /2)
  marker_names <- paste0("locus", 1:num_loci)
  
  marker_types <- rep(0, num_loci)
  
  allelic_dropout_rates <- rep(0.01, num_loci)   #based on REPave  being >99
  genotyping_error_rates <- rep(0.01, num_loci)  #checked
  
  marker_info <- rbind(
    marker_names,
    marker_types,
    allelic_dropout_rates,
    genotyping_error_rates
  )
  
  write.table(marker_info, file = paste0("./data/marker_info_", input_name, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
  
  marker_info  #this creates a col for each loci. Should be half colony_data_lines -1
  ncol(marker_info)
  
  
  individual_names <- indNames(genind_obj)
  
  colony_data_lines <- character(nrow(genotype_matrix))
  for (i in 1:nrow(genotype_matrix)) {
    individual_data <- c(individual_names[i])  # Start with individual ID
    
    for (j in seq(1, ncol(genotype_matrix), by = 2)) {
      allele1_dosage <- genotype_matrix[i, j]
      allele2_dosage <- genotype_matrix[i, j + 1]
      
      if (is.na(allele1_dosage) || is.na(allele2_dosage)) {
        alleles <- c(0, 0)
      } else if (allele1_dosage == 2 && allele2_dosage == 0) {
        alleles <- c(1, 1)
      } else if (allele1_dosage == 1 && allele2_dosage == 1) {
        alleles <- c(1, 2)
      } else if (allele1_dosage == 0 && allele2_dosage == 2) {
        alleles <- c(2, 2)
      } else {
        alleles <- c(0, 0)
      }
      
      individual_data <- c(individual_data, alleles)
    }
    
    colony_data_lines[i] <- paste(individual_data, collapse = " ")
  }
  
  
  writeLines(colony_data_lines, con = paste0("./data/col_input_", input_name, ".txt"))
  
  length(unlist(strsplit(colony_data_lines[1], " ")))   #this creates 200 +1 col. 2 col per loci
  ((length(unlist(strsplit(colony_data_lines[1], " ")))) - 1)/2
}

colony_data(data_genind_adult)
colony_data(data_genind_adult_subset1)
colony_data(data_genind_adult_subset2)





data1 <- read.table(file = "./data/marker_info_data_genind_adult_subset1.txt", header = TRUE, dec = ",", na.strings = c("", ".", "na")) ## replace XXX to document name.
data2 <- read.table(file = "C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects/3_Heron_Platy_ms/Colony/heron_sub1_inb_clones/heron_sub1_inb_clones.ErrorRate", 
                    header = TRUE, sep = ",", na.strings = c("", ".", "na")) ## replace XXX to document name.

error_rates <- data2 %>% select(MarkerID, DropRateEst)

data1_t <- as.data.frame(t(data1))
data1_t <- tibble::rownames_to_column(data1_t, "MarkerID")

updated_data <- data1_t %>% left_join(error_rates, by = "MarkerID") %>% mutate(V2 = ifelse(!is.na(DropRateEst), DropRateEst, V2)) %>% 
  select(-DropRateEst)


data1_updated <- as.data.frame(t(updated_data))
colnames(data1_updated) <- NULL
rownames(data1_updated) <- NULL
write.table(data1_updated, file = './data/marker_info_data_genind_adult_subset1_corr.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")


