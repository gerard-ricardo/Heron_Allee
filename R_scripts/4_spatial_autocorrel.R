# spatial autocorrelation -------------------------------------------------

##notes
# change to genid

## genetic distances

#remove 12 test  (removing 12 didn;t change the spatial relationship)
# id <- data_genind_adult$other$ind.metrics$id
# data_genind_adult <- data_genind_adult[id != c('pd12.a.1', 'pd12.a.2'), ]


# filter genID for single reps
# data_genind_adult$other$ind.metrics$stage
# # Check the population assignments
# table(data_genind_adult@pop)
# # Remove subsamples
# unique_indices <- !duplicated(data_genind_adult@pop)
# data_genind_adult_unique <- data_genind_adult[unique_indices, ]
# table(data_genind_adult_unique@pop)
# # Collapse population
#data_genind_adult_unique@pop <- factor(rep("population1", nrow(data_genind_adult_unique@tab))) #combine all

# Calculate the genetic distance matrix using Smouse method
data_genind_adult_unique@pop <- factor(rep("population1", nrow(data_genind_adult_unique@tab))) #combine all
genetic_dist_matrix <- gd.smouse(data_genind_adult_unique, verbose = TRUE)  #Smouse and Peakall (1999) is a method used to quantify the
#genetic distance (dissimilarity) between pairs of individuals. This approach is often used in population genetics to understand 
#genetic structure, diversity, and relatedness among individuals within and between populations.
genetic_dist_matrix1 = as.matrix(genetic_dist_matrix)


## spatial distances
# # filter data_gl_filtered_adult for single reps
# # Identify unique indices based on population
# unique_indices_gl <- !duplicated(data_gl_filtered_adult@pop)
# # Create a new filtered object with these unique indices
# data_gl_filtered_adult_unique <- data_gl_filtered_adult[unique_indices_gl, ]
# # Collapse populations into one
# #data_gl_filtered_adult@pop <- factor(rep("population1", nrow(data_gl_filtered_adult@tab)))
# # Check the population assignments after collapsing
# table(data_gl_filtered_adult_unique@pop)

# Extract spatial coordinates for adults only
coordinates <- data_genind_adult_unique@other$ind.metrics %>%
  dplyr::select(lat, lon) %>%
  dplyr::mutate(id = rownames(data_genind_adult_unique@other$ind.metrics)) %>%
  # dplyr::filter(id %in% data_gl_filtered_adult_unique@ind.names) %>%
  # dplyr::arrange(match(id, data_gl_filtered_adult_unique@ind.names)) %>%
  dplyr::select(lat, lon)
# Ensure the coordinates are in the same order as the individuals in the genind object
#coordinates <- coordinates[rownames(data_gl_filtered_adult_unique@tab), ]
colnames(coordinates) <- c("lat", "lon")
# Convert to sf object
library(sf)
coordinates_sf <- st_as_sf(coordinates, coords = c("lon", "lat"), crs = 4326)
# Transform to UTM (zone 56S for Heron Island)
utm_crs <- st_crs("+proj=utm +zone=56 +south +datum=WGS84 +units=m +no_defs")
coordinates_utm <- st_transform(coordinates_sf, utm_crs)
# Extract UTM coordinates as a matrix
coordinates_matrix <- st_coordinates(coordinates_utm)
# Calculate the Euclidean distance matrix based on coordinates
euclidean_dist_matrix <- dist(coordinates_matrix)
euclidean_dist_matrix <- as.matrix(euclidean_dist_matrix)

# Perform spatial autocorrelation analysis (on raw)
bin = 8
(spatial_autocor_results <- spautocor(gen.m = genetic_dist_matrix1, eucl.m = euclidean_dist_matrix, bins = bin, 
                                      shuffle = F))
#Moran's I coefficient, values between -1 and 1

# Plot the results
df1 = data.frame(y = spatial_autocor_results$r, x = spatial_autocor_results$bin)
plot(df1$y~ df1$x, type = "b", 
     xlab = "Distance Classes (m)", ylab = "Autocorrelation Coefficient (r)",
     main = "Spatial Autocorrelation Analysis")

ggplot(df1, aes(x = x, y = y)) +
  geom_line(color = "steelblue") +
  geom_point(color = "steelblue", size = 3) +
  labs(x = "Distance classes (m)", y = "Autocorrelation coefficient (r)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 1) +
  theme_sleek2()

# run on 1000 perms
# Function to calculate spatial autocorrelation with shuffling
calculate_spatial_autocorrelation <- function(genetic_matrix, spatial_matrix, bins) {
  spautocor(gen.m = genetic_matrix, eucl.m = spatial_matrix, bins = bin, shuffle = T)$r
}

# Run permutations and store the results
set.seed(123) # For reproducibility
num_permutations <- 5000
permutation_results <- replicate(num_permutations, {
  calculate_spatial_autocorrelation(genetic_dist_matrix1, euclidean_dist_matrix, bins = bin)
})

# Convert permutation results to a long format data frame for ggplot
permutation_df <- data.frame((permutation_results))
permutation_df$dist <- spatial_autocor_results$bin
head(permutation_df)
data_long <- permutation_df %>% tidyr::pivot_longer(-dist, names_to = "perm", values_to = "resp") %>% 
  arrange(perm, dist) %>% data.frame()

# Function to calculate the 95% confidence intervals
calculate_ci <- function(permutation_results, ci_level = 0.95) {
  apply(permutation_results, 1, function(x) {
    quantile(x, probs = c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2))
  })
}

# Calculate the 95% CI for the permutations
ci_level <- 0.95
ci <- calculate_ci(permutation_results, ci_level)
ci_df <- data.frame(bin = spatial_autocor_results$bin, 
                    lower = ci[1, ], 
                    upper = ci[2, ])

# plot
p4 = ggplot() +
  geom_ribbon(data = ci_df, mapping = aes(x = bin, ymin = lower, ymax = upper), fill = "lightgrey", alpha = 0.5) +
  geom_line(data_long, mapping = aes(x = dist, y = resp, group = perm),alpha = 0.05, colour = "steelblue") +
  geom_hline(aes(yintercept = 0), col = 'grey30', linetype = "dashed") +
  geom_line(data = df1, mapping = aes(x = x, y = y), color = "salmon", size = 1) +
  labs(
    x = "Distance classes (m)",
    y = "Autocorrelation coefficient (r)",
  ) +
  theme_sleek2()
p4



# ## attempt multiple individual runs
# # Function to calculate the spatial autocorrelation for a given permutation
# calculate_spatial_autocorrelation <- function(genetic_matrix, spatial_matrix, bins) {
#   shuffled_genetic_matrix <- genetic_matrix[sample(nrow(genetic_matrix)), sample(ncol(genetic_matrix))]
#   spautocor(gen.m = shuffled_genetic_matrix, eucl.m = spatial_matrix, bins = bins)$r
# }
# 
# # Run permutations and store the results
# set.seed(123) # For reproducibility
# num_permutations <- 1000
# permutation_results <- replicate(num_permutations, {
#   calculate_spatial_autocorrelation(genetic_dist_matrix1, euclidean_dist_matrix, bins = 8)
# })
# dim(permutation_results)
# 
# # Compare the observed results to the permutation results
# observed_r <- spatial_autocor_results$r
# p_values <- apply(permutation_results, 1, function(x) mean(x >= observed_r))
# 
# # Plot the observed vs. permuted distributions
# hist(permutation_results[1, ], breaks = 30, main = "Distribution of Permuted Autocorrelation Coefficients",
#      xlab = "Autocorrelation Coefficient (r)", col = "lightgrey", border = "white")
# abline(v = observed_r[1], col = "red", lwd = 2)
# legend("topright", legend = c("Observed r"), col = c("red"), lwd = 2)
# 
# # Convert permutation results to a long format data frame for ggplot
# permutation_df <- data.frame((permutation_results))
# permutation_df$bin <- 1:nrow(permutation_df)
# data1_long = permutation_df %>% tidyr::pivot_longer(-bin,  names_to = "perm" ,values_to = "resp") %>% 
#   arrange(perm, bin) %>% data.frame()  #keep vec.x, add all other columns to factors , add all their values to meas)
# 
# # Plot the observed vs. permuted distributions
# str(data1_long)
# data1_long$perm <- as.factor(as.character(data1_long$perm))
# library(ggplot2)
# ggplot() +
#   geom_line(data1_long, mapping = aes(x = bin, y = resp, group = perm), alpha = 0.9, colour = "grey") +
#   geom_line(data = data.frame(bin = 1:length(observed_r), value = observed_r), 
#             aes(x = bin, y = value), color = "red", size = 1.5) +
#   labs(
#     x = "Distance classes (m)",
#     y = "Autocorrelation coefficient (r)",
#     title = "Observed and Permuted Autocorrelation Coefficients"
#   ) +
#   theme_minimal()