# Hardy-Weinberg equilibrium and heterozygote excess----------------------------------------------

# Perform HWE test for each locus
hwe_results <- hw.test(data_genind_filtered, B = 0)  # B is the number of permutations
#Significant deviations can indicate factors such as inbreeding, genetic drift, selection, or self-fertilisation.

##Identify loci with heterozygote excess:
# Extract p-values and heterozygote excess information
hwe_pvalues <- hwe_results[, 3]
# Adjust p-values for multiple testing using Bonferroni correction
p_adjusted <- p.adjust(hwe_pvalues, method = "fdr")
# Identify loci with significant heterozygote excess after adjustment
significant_loci <- which(p_adjusted < 0.05)
# Extract observed and expected heterozygosity for these loci
obs_het <- summary(data_genind_filtered)$Hobs
exp_het <- summary(data_genind_filtered)$Hexp
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
