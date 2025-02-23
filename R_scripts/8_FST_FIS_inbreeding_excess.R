









pop_info <- data_genind_adult_unique@other$ind.metrics$cluster # population information for individuals
genotypes <- tab(data_genind_adult_unique)  # genotype matrix (one locus per column)
bs_df <- data.frame(Population = pop_info, genotypes)
(bs.nc <- basic.stats(bs_df))


bs.nc <- basic.stats(data_genind_adult_subset1)
bs.nc


gl.report.heterozygosity(data_gl_filtered_adult)


wc(data_genind_adult[, -2])  #Computes Weir and Cockerham estimates of Fstatistics

betas(data_genind_adult)
betas_values <- betas(data_genind_adult)$betaiovl
sorted_betas <- sort(betas_values)
sorted_betas
barplot(sorted_betas,
        main = "Population-specific Fst Values",
        ylab = "Fst", xlab = "Population", col = "steelblue",
        las = 2
)


(bs.nc <- basic.stats(data_genind_adult_subset1))
gl.report.heterozygosity(data_genind_adult_subset1)

(bs.nc <- basic.stats(data_genind_adult_cluster2))
(bs.nc <- basic.stats(data_genind_adult_subset3))



data_genind_adult_unique@pop <- data_genind_adult_unique@other$ind.metrics$cluster
values <- inbreeding(data_genind_adult_unique)
(median__values <- lapply(values, median))
(df <- do.call(rbind, median__values)) 

values_all <- inbreeding(data_genind)
(median__values_all <- lapply(values_all, median))

values_1 <- inbreeding(data_genind_adult_subset1)
(median_values_1 <- lapply(values_1, median))
(df1 <- do.call(rbind, median_values_1))  #

values_2 <- inbreeding(data_genind_adult_subset2)
(median_values_2 <- lapply(values_2, median))
(df2 <- do.call(rbind, median_values_2))  #

ids <- rep(names(values), times = sapply(values, length))  # Repeat each name according to the length of each list element
values_flat <- unlist(values, use.names = FALSE)  # Unlist without preserving names to avoid auto-generated names
df <- data.frame(id = ids, bbb = values_flat)
df$id <- as.factor(as.character(df$id))
calculate_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

mode_df <- df %>% group_by(id) %>% summarise(mode = calculate_mode()) %>% data.frame()
str(df)
levels(df$id)

(med_df <- df %>% dplyr::group_by(id) %>% dplyr::summarise(med = median(bbb, na.rm = TRUE)) %>% data.frame())
median(med_df$med)
plot(density(med_df$med), xlab  = 'Inbreeding coef')

str(df)
str(mode_df)
df <- left_join(df, med_df, by = "id") %>%  arrange(med)  # Joining the mode values back to the original dataframe based on id
range(df$bbb)
library(ridgeline)
ridgeline(df$bbb, df$id, mode = T) 

sorted_ <- arrange(med_df, med)
sorted_

height <- sorted_$med

barplot(height,names.arg = sorted_$id, main = "Individual-specific  Values", ylab = "",  col = "blue", las = 2)



(stats_adult <- basic.stats(data_genind_parents))
stats_adult$Fis
(indiv_mean_adult <- colMeans(stats_adult$Fis, na.rm = TRUE))
mean(indiv_mean_adult, na.rm = TRUE)
(stats_progeny <- basic.stats(data_genind_progeny))
stats_progeny$Fis
(indiv_mean_prog <- colMeans(stats_progeny$Fis, na.rm = TRUE))
mean(indiv_mean_prog, na.rm = TRUE)

shapiro.test(indiv_mean_adult)  # For adults
shapiro.test(indiv_mean_prog)   # For progeny.
combined_data <- c(indiv_mean_adult, indiv_mean_prog)
group_factor <- factor(c(rep("Adults", length(indiv_mean_adult)), 
                         rep("Progeny", length(indiv_mean_prog))))
car::leveneTest(combined_data, group_factor)
t.test(indiv_mean_adult, indiv_mean_prog, var.equal = TRUE)


library(inbreedR)
data('mouse_snps')
mat_0_1_coded
str(mat_0_1_coded)
check_data(mat_0_1_coded, num_ind = 35, 1454)
g2 = g2_snps(mat_0_1_coded, nperm = 100, nboot =100, CI = 0.95)
plot(g2)
r2_hf(mat_0_1_coded, nboot = 100, type = 'msats')
r2_Wf(mat_0_1_coded, nboot = 100, type = 'msats')



hwe_results <- hw.test(data_genind_adult, B = 0)  # B is the number of permutations

hwe_pvalues <- hwe_results[, 3]
p_adjusted <- p.adjust(hwe_pvalues, method = "fdr")
significant_loci <- which(p_adjusted < 0.05)
obs_het <- summary(data_genind_adult)$Hobs
exp_het <- summary(data_genind_adult)$Hexp
heterozygote_excess <- obs_het > exp_het
loci_het_excess <- which(heterozygote_excess & (p_adjusted < 0.05))
loci_het_excess
length(loci_het_excess)

heterozygote_deficit <- obs_het < exp_het
loci_het_deficit <- which(heterozygote_deficit & (p_adjusted < 0.05))
print(loci_het_deficit)
num_loci_het_deficit <- length(loci_het_deficit)
print(num_loci_het_deficit)


data_plot <- data.frame(Expected = exp_het, Observed = obs_het, Color = "black")
data_plot$Color[loci_het_excess] <- "steelblue2"
data_plot$Color[loci_het_deficit] <- "orchid4"

ggplot(data_plot, aes(x = Expected, y = Observed, color = Color)) +
  geom_point(alpha = 0.6, size = 2) +
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
ggsave( filename = 'heron_pHo_vs_He.tiff',  path = "./plots", device = "tiff",  width = 5, height = 5)  #this often works better than pdf


Go <- apply(tab(data_genind_adult_unique), 1, function(x) length(unique(x))) # Number of unique genotypes per locus
Ge <- apply(tab(data_genind_adult_unique), 1, function(locus) {
  p <- sum(locus * (1:length(locus)-1)) / sum(locus) # Calculate allele frequency p for bi-allelic loci
  1 - (p^2 + (1-p)^2) # HW expected heterozygosity for a bi-allelic locus
})




ind_names <- indNames(data_genind)
mum_index <- which(ind_names == 'pd15.a.1')
maternal_plant <- data_genind[mum_index, ]

progeny_indices <- grep('^pd15.l\\.', ind_names, value = FALSE) # Use correct pattern to match progeny names
progeny_indices <- progeny_indices[progeny_indices != mum_index]
progeny_genind <- data_genind[progeny_indices, ]

summary_progeny <- summary(progeny_genind)
Ho <- summary_progeny$Hobs

He <- summary_progeny$Hexp

Ho_avg <- mean(Ho, na.rm = TRUE) # Average observed heterozygosity
He_avg <- mean(He, na.rm = TRUE) # Average expected heterozygosity

cat("Average observed heterozygosity (Ho):", Ho_avg, "\n")
cat("Average expected heterozygosity (He):", He_avg, "\n")

F <- 1 - (Ho_avg / He_avg)

S <- (2 * F) / (1 + F)

cat("Inbreeding coefficient (F):", F, "\n")
cat("Selfing rate (S):", S, "\n")





(Ho_avg = mean(data1$HObs, na.rm = T)) 
hist(data1$HObs)
(Ho_med = median(data1$HObs, na.rm = T)) 

(He_avg = mean(data1$HExp, na.rm = T))
hist(data1$HExp)
(He_med = median(data1$HExp, na.rm = T))



data1$F_IS <- (data1$HExp - data1$HObs) / data1$HExp
hist(data1$F_IS)
med_F_IS <- median(data1$F_IS, na.rm = TRUE)
med_F_IS




(S1 <- (2 * F1) / (1 + F1))

