# Fst (between pops) and FIS (Inbreeding Coefficient) ---------------------------------------------------------------------

##Note that many inbreeding calcs here might be biased because of genotyping errors increase homozygotes. 


#basic relatedness stats
bs.nc <- basic.stats(data_genind_adult)
bs.nc
#Ho(0-1)  Hs(0-1) Ht(0-1) Dst(-)  Htp(0-1) Dstp    Fst(0-1) Fstp(0-1)   Fis(-1-1)  Dest(0-1) 
#0.0842   0.0533  0.1553  0.1020  0.1619   0.1086  0.6569   0.6708     -0.5808     0.1147 
#bit odd FIS is so low here

# Weir and Cockerham estimates
wc(data_genind_adult[, -2])  #Computes Weir and Cockerham estimates of Fstatistics
#agrees with above
# The high FST suggests a high level of differentiation among populations, while the negative FIS suggests a
# possible excess of heterozygotes (outcrossing) within the populations.

## Fst
# Calculate population-specific Fst values for the filtered genetic data
betas(data_genind_adult)
# Extract the population-specific Fst values
betas_values <- betas(data_genind_adult)$betaiovl
sorted_betas <- sort(betas_values)
sorted_betas
barplot(sorted_betas,
        main = "Population-specific Fst Values",
        ylab = "Fst", xlab = "Population", col = "blue",
        las = 2
)

# Overall Inbreeding Coefficient.0 means random mating.  Positive values indicate a deficiency of heterozygotes,
# suggesting inbreeding. Negative values indicate an excess of heterozygotes, suggesting outcrossing.
# Pos values might be related to clones and self fert (but probably not)




# inbreeding coefs --------------------------------------------------------

## Fis  (-1 to 1)
# Calculate inbreeding coefficients (FIS)
fis_values <- inbreeding(data_genind_adult)
#str(fis_values)
(median_fis_values <- lapply(fis_values, median))


ids <- rep(names(fis_values), times = sapply(fis_values, length))  # Repeat each name according to the length of each list element
fis_values_flat <- unlist(fis_values, use.names = FALSE)  # Unlist without preserving names to avoid auto-generated names
fis_df <- data.frame(id = ids, fis = fis_values_flat)
str(fis_df)
# Define a function to calculate mode
calculate_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Compute mode for each id and create a data frame of modes
mode_df <- fis_df %>%
  group_by(id) %>%
  summarise(mode = calculate_mode(fis), .groups = 'drop') %>% data.frame()

# Join the mode back to the original dataframe
fis_df <- fis_df %>%
  left_join(mode_df, by = "id") %>%  arrange(mode)  # Joining the mode values back to the original dataframe based on id

range(fis_df$fis)
#remotes::install_github("R-CoderDotCom/ridgeline@main")
library(ridgeline)
ridgeline(fis_df$fis, fis_df$id, mode = T) 

# Sort the Fst values from low to high
sorted_fis <- arrange(mode_df, mode)
sorted_fis

# Extract the mode values from the sorted data frame
height <- sorted_fis$mode

# Create the bar plot with sorted Fis values
barplot(height,
        names.arg = sorted_fis$id,
        main = "Population-specific Fis Values",
        ylab = "Fis",  col = "blue",
        las = 2
)


#########(Identify disequilibrium (ranges 0 - 1)
library(inbreedR)
data('mouse_snps')
mat_0_1_coded
str(mat_0_1_coded)
check_data(mat_0_1_coded, num_ind = 35, 1454)
g2 = g2_snps(mat_0_1_coded, nperm = 100, nboot =100, CI = 0.95)
plot(g2)
r2_hf(mat_0_1_coded, nboot = 100, type = 'msats')
r2_Wf(mat_0_1_coded, nboot = 100, type = 'msats')


# Hardy-Weinberg equilibrium and heterozygote excess----------------------------------------------


##NOTE: Ho and He might be affect by high null alleles. Might be best to use values created from Cervus as this adjust for nul alleles. 

# Perform HWE test for each locus
hwe_results <- hw.test(data_genind, B = 0)  # B is the number of permutations
#Significant deviations can indicate factors such as inbreeding, genetic drift, selection, or self-fertilisation.

##Identify loci with heterozygote excess:
# Extract p-values and heterozygote excess information
hwe_pvalues <- hwe_results[, 3]
# Adjust p-values for multiple testing using Bonferroni correction
p_adjusted <- p.adjust(hwe_pvalues, method = "fdr")
# Identify loci with significant heterozygote excess after adjustment
significant_loci <- which(p_adjusted < 0.05)
# Extract observed and expected heterozygosity for these loci
obs_het <- summary(data_genind)$Hobs
exp_het <- summary(data_genind)$Hexp
# Check for heterozygote excess
heterozygote_excess <- obs_het > exp_het
# Loci with significant heterozygote excess
loci_het_excess <- which(heterozygote_excess & (p_adjusted < 0.05))
# Print loci with heterozygote excess
print(loci_het_excess)
length(loci_het_excess)

# Check for heterozygote deficit
heterozygote_deficit <- obs_het < exp_het
# Loci with significant heterozygote deficit
loci_het_deficit <- which(heterozygote_deficit & (p_adjusted < 0.05))
# Print loci with heterozygote deficit
print(loci_het_deficit)
# Number of loci with heterozygote deficit
num_loci_het_deficit <- length(loci_het_deficit)
print(num_loci_het_deficit)


# Plot observed vs expected heterozygosity
# plot(obs_het ~exp_het, xlab = "Expected Heterozygosity", ylab = " Observed Heterozygosity")
# abline(0, 1, col = "red")
# points(obs_het[loci_het_excess]~ exp_het[loci_het_excess], col = "blue", pch = 19)
# Points below the red line indicate loci with a deficiency of heterozygotes (observed < expected), suggesting 
#inbreeding or other factors.Points above the red line indicate loci with an excess of heterozygotes (observed > 
#expected), which can suggest outcrossing, heterozygote advantage, or self-fertilisation.

# Create a data frame from the observed and expected heterozygosity
data_plot <- data.frame(Expected = exp_het, Observed = obs_het, Color = "black")
# Mark the significant loci with heterozygote excess
data_plot$Color[loci_het_excess] <- "steelblue4"
# Mark the significant loci with heterozygote deficit
data_plot$Color[loci_het_deficit] <- "orchid4"

# Create the plot
ggplot(data_plot, aes(x = Expected, y = Observed, color = Color)) +
  geom_point(alpha = 0.6) +
  scale_color_identity() +  # Use actual colors stored in 'Color' column
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = "Observed vs Expected Heterozygosity",
    x = "Expected Heterozygosity",
    y = "Observed Heterozygosity"
  ) +
  theme_sleek2() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )



# A primer of conservation genetics equations pg 186 (working but issues)----------------------

#the F and S are implausible. Possible genotyping errors 

# Assuming 'data_genind' is your genind object containing all individuals
ind_names <- indNames(data_genind)
# Find the index of the maternal plant 'pd15.a.1'
mum_index <- which(ind_names == 'pd15.a.1')
# Extract the maternal plant
maternal_plant <- data_genind[mum_index, ]

# Extract progeny of 'pd15.a.1'
# Assuming progeny are named in a pattern like 'pd15.a.1.<suffix>'
progeny_indices <- grep('^pd15.l\\.', ind_names, value = FALSE) # Use correct pattern to match progeny names
# Exclude the maternal plant from progeny
progeny_indices <- progeny_indices[progeny_indices != mum_index]
# Extract the progeny
progeny_genind <- data_genind[progeny_indices, ]

# Calculate observed heterozygosity (Ho) for progeny
summary_progeny <- summary(progeny_genind)
Ho <- summary_progeny$Hobs

# Calculate expected heterozygosity (He) for progeny
He <- summary_progeny$Hexp

# Calculate average observed and expected heterozygosity
Ho_avg <- mean(Ho, na.rm = TRUE) # Average observed heterozygosity
He_avg <- mean(He, na.rm = TRUE) # Average expected heterozygosity

# Check if Ho_avg and He_avg are correctly calculated
cat("Average observed heterozygosity (Ho):", Ho_avg, "\n")
cat("Average expected heterozygosity (He):", He_avg, "\n")

# Compute the inbreeding coefficient (F)
F <- 1 - (Ho_avg / He_avg)

# Determine the selfing rate (S)
S <- (2 * F) / (1 + F)

# Print results
cat("Inbreeding coefficient (F):", F, "\n")
cat("Selfing rate (S):", S, "\n")



# cervus analysis ---------------------------------------------------------

#Run output  cervus script first and ensure the most up to data pd_afa_out2.txt file in 1b_2022heron_seq_process_other.R
#These values are senstive to null allelel filtering

(Ho_avg = mean(data1$HObs, na.rm = T)) 
hist(data1$HObs)
(Ho_med = median(data1$HObs, na.rm = T)) 

(He_avg = mean(data1$HExp, na.rm = T))
hist(data1$HExp)
(He_med = median(data1$HExp, na.rm = T))

## this indicates hetero excess, but might also be from filtering null allees


# Calculate F_IS (inbreeding coefficient) for each locus
data1$F_IS <- (data1$HExp - data1$HObs) / data1$HExp
hist(data1$F_IS)
# Calculate the mean F_IS
med_F_IS <- median(data1$F_IS, na.rm = TRUE)
# Display the mean F_IS
med_F_IS




# Determine the selfing rate (S) - only work with positive F (i.e evidence of  inbreeding)
(S1 <- (2 * F1) / (1 + F1))

