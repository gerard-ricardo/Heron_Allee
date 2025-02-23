





data_genind_adult_unique@pop <- factor(rep("population1", nrow(data_genind_adult_unique@tab))) #combine all
genetic_dist_matrix <- gd.smouse(data_genind_adult_unique, verbose = TRUE)  #Smouse and Peakall (1999) is a method used to quantify the
genetic_dist_matrix1 = as.matrix(genetic_dist_matrix)



coordinates <- data_genind_adult_unique@other$ind.metrics %>%
  dplyr::select(lat, lon) %>%
  dplyr::mutate(id = rownames(data_genind_adult_unique@other$ind.metrics)) %>%
  dplyr::select(lat, lon)
colnames(coordinates) <- c("lat", "lon")
library(sf)
coordinates_sf <- st_as_sf(coordinates, coords = c("lon", "lat"), crs = 4326)
utm_crs <- st_crs("+proj=utm +zone=56 +south +datum=WGS84 +units=m +no_defs")
coordinates_utm <- st_transform(coordinates_sf, utm_crs)
coordinates_matrix <- st_coordinates(coordinates_utm)
euclidean_dist_matrix <- dist(coordinates_matrix)
euclidean_dist_matrix <- as.matrix(euclidean_dist_matrix)

bin = 8
(spatial_autocor_results <- spautocor(gen.m = genetic_dist_matrix1, eucl.m = euclidean_dist_matrix, bins = bin, 
                                      shuffle = F))

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

calculate_spatial_autocorrelation <- function(genetic_matrix, spatial_matrix, bins) {
  spautocor(gen.m = genetic_matrix, eucl.m = spatial_matrix, bins = bin, shuffle = T)$r
}

set.seed(123) # For reproducibility
num_permutations <- 5000
permutation_results <- replicate(num_permutations, {
  calculate_spatial_autocorrelation(genetic_dist_matrix1, euclidean_dist_matrix, bins = bin)
})

permutation_df <- data.frame((permutation_results))
permutation_df$dist <- spatial_autocor_results$bin
head(permutation_df)
data_long <- permutation_df %>% tidyr::pivot_longer(-dist, names_to = "perm", values_to = "resp") %>% 
  arrange(perm, dist) %>% data.frame()

calculate_ci <- function(permutation_results, ci_level = 0.95) {
  apply(permutation_results, 1, function(x) {
    quantile(x, probs = c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2))
  })
}

ci_level <- 0.95
ci <- calculate_ci(permutation_results, ci_level)
ci_df <- data.frame(bin = spatial_autocor_results$bin, 
                    lower = ci[1, ], 
                    upper = ci[2, ])

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



