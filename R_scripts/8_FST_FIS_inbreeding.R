# Fst (between pops) and FIS (Inbreeding Coefficient) ---------------------------------------------------------------------

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

##