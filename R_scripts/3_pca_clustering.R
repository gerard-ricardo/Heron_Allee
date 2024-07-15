

# objects -----------------------------------------------------------------

data_gl_filtered_adult  #only adults
data_gl_filtered  #all stages



# adult only ---------------------------------------------------------------------

#quick plot
pca = gl.pcoa(data_gl_filtered_adult)
gl.pcoa.plot(glPca = pca, data_gl_filtered_adult)

# PCA Analysis
pca_data <- tab(data_gl_filtered_adult, freq = TRUE, NA.method = "mean") %>% na.omit() # Convert to tabular format and omit NAs
pca <- dudi.pca(pca_data, center = TRUE, scale = FALSE, nf = 2, scannf = FALSE) # Perform PCA
pca_complete <- data.frame(pca$li, pop = data_gl_filtered_adult$pop) # Combine PCA results with population data
#use for adults
#pca_complete <- data.frame(pca$li) # Combine PCA results with population data

# Explained variance
(explained_variance <- pca$eig / sum(pca$eig) * 100)
scree_plot <- data.frame(PC = 1:length(explained_variance), Variance = explained_variance)

ggplot(scree_plot, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_line(aes(y = cumsum(Variance)), group = 1, color = "red") +
  geom_point(aes(y = cumsum(Variance)), color = "red") +
  labs(title = "Scree Plot", x = "Principal Component", y = "Percentage of Variance Explained") +
  theme_sleek2()

# Hopkins statistic
set.seed(123) # for reproducibility
(hopkins_stat <- hopkins(pca_data, n = nrow(pca_data) - 1))
# Calculated values 0-0.3 indicate regularly-spaced data. Values around 0.5 indicate random data. Values 0.7-1 indicate clustered data.
#PD = 0.22


# K-means clustering
set.seed(123) # for reproducibility
kmeans_result <- kmeans(pca_data, centers = 3, nstart = 25)
individuals_in_cluster3 <- which(kmeans_result$cluster == 3) #find indiv in each cluster
silhouette_score <- silhouette(kmeans_result$cluster, dist(pca_data))
summary(silhouette_score)
plot(silhouette_score)
pca_complete$Cluster <- as.factor(kmeans_result$cluster)
#PD: cluster 3 is quite strong, others poor to mod. 

# DBSCAN clustering
# Find the appropriate eps value using kNNdistplot
kNNdistplot(pca_data, k = 5)
elbow = 8.9 # Place this at the elbow of the line
abline(h = elbow, col = "red", lty = 2)  
library(dbscan)
# Function to perform DBSCAN clustering and plot results
perform_dbscan <- function(pca_data, pca_complete, eps_value, min_pts = 5) {
  dbscan_result <- dbscan(pca_data, eps = eps_value, minPts = min_pts)
  cluster_col_name <- paste0("Cluster_dbscan_", eps_value)
  pca_complete[[cluster_col_name]] <- as.factor(dbscan_result$cluster)
  plot <- ggplot(pca_complete, aes_string(x = "Axis1", y = "Axis2", color = cluster_col_name)) +
    geom_point(alpha = 0.6) +
    labs(title = paste("PCA Plot with DBSCAN Clusters (eps =", eps_value, ")"),
         x = "Principal Component 1",
         y = "Principal Component 2") +
    theme_minimal()
  silhouette_score <- silhouette(dbscan_result$cluster, dist(pca_data))
  print(dbscan_result)
  print(summary(silhouette_score))
  return(plot)
}

eps_values <- elbow 
for (eps in eps_values) {
  plot <- perform_dbscan(pca_data, pca_complete, eps)
  print(plot)
}
#The clustering contains 2 cluster(s) and 2 noise points (id = 12).


# plotting
#PD
pca_complete <- pca_complete %>%
  mutate(
    Stage = ifelse(str_detect(row.names(pca_complete), "\\.a\\."), "Adult", "Larva"),  #add stage
    RepID = str_extract(row.names(pca_complete), "(?<=\\.)\\d+$"),
    NewID = paste0(Stage,  pop, "_", RepID)
  )

data1 <- dplyr::arrange(pca_complete, Axis1) # 
pca_complete <- pca_complete %>% mutate(across(c(Stage, MumID, RepID, NewID), as.factor))
pca_complete1 = pca_complete %>% select(Cluster, pop) %>% rename(id = pop)  #for bathymetry file 
str(pca_complete)
my_palette <- c(
  "dodgerblue", "firebrick", "mediumseagreen", "orchid", "darkorange", "gold",
  "skyblue", "sandybrown", "palevioletred", "mediumturquoise", "khaki",
  "darkslategray", "plum", "lightslategray", "limegreen", "cornflowerblue",
  "tomato",'red'
)

#color individuals
t2 <- ggplot(pca_complete, aes(x=Axis1, y=Axis2, group=NewID)) +
  geom_point(aes(fill=pop, color=ifelse(grepl("Larva", Stage), "red", 'black'), alpha = 0.8), shape=21, size=4, stroke=1, alpha = 0.8) +
  scale_fill_manual(values=my_palette) +
  scale_color_manual(values=c("red"= "red", "black"= "black")) +
  theme_minimal() +
  labs(x = "PCA1", y = "PCA2", color = "Stage", fill = "Population")  # Add labels to the axes and legend
t2

#per cluster
t2 <- ggplot(pca_complete, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(color = factor(Cluster)), shape = 22, 
             size = 3, stroke = 1, alpha = 0.7, position = position_jitter(width = 0.1, height = 0.1)
  ) +
  geom_text_repel(aes(label = NewID), size = 3, max.overlaps = 38, point.padding = 0.5, box.padding = 0.5) +
  scale_color_manual(values = c("1" = "dodgerblue", "2" = "salmon", "3" = "mediumseagreen")) +
  stat_ellipse(aes(x = Axis1, y = Axis2, group = Cluster, color = Cluster), level = 0.95, linetype = 2, size = 1) + # Add ellipses around clusters
  theme_sleek2() +
  labs(
    x = paste0("PCA1 (", round(explained_variance[1], 2), "%)"),
    y = paste0("PCA2 (", round(explained_variance[2], 2), "%)"),
    color = "Cluster", fill = "Population", shape = "Stage"
  ) 
t2

# Convert the ggplot to an interactive plotly plot
t2_interactive <- ggplotly(t2)
t2_interactive



# 3d pca ------------------------------------------------------------------

# Conduct PCA with three components
pca1 <- dudi.pca(pca_data, center = TRUE, scale = FALSE, nf = 3, scannf = FALSE) # Calculate the first three PCs

# Combine PCA results with population data and extract three principal components
pca_complete1 <- data.frame(pca1$li, pop = data_gl_filtered$pop)
colnames(pca_complete1)[1:3] <- c("Axis1", "Axis2", "Axis3") # Rename columns to Axis1, Axis2, Axis3 for clarity

# Augment pca_complete1 with additional details
pca_complete1 <- pca_complete1 %>%
  mutate(
    Stage = ifelse(str_detect(row.names(pca_complete1), "\\.a\\."), "Adult", "Larva"),
    MumID = str_extract(row.names(pca_complete1), "(?<=pd)\\d+"),
    RepID = str_extract(row.names(pca_complete1), "(?<=\\.)\\d+$"),
    NewID = paste0(Stage, MumID, "_", RepID)
  ) %>%
  mutate(across(c(Stage, MumID, RepID, NewID), as.factor))  # Convert relevant columns to factors

# Map the Cluster to a specific color in my_palette
cluster_colors <- my_palette[pca_complete1$Cluster]

plot_3d <- plot_ly(data = pca_complete1, x = ~Axis1, y = ~Axis2, z = ~Axis3, type = 'scatter3d', mode = 'markers',
                   marker = list(size = 5, color = cluster_colors, colorscale = list(c(0, 1), c(min(cluster_colors), max(cluster_colors))), opacity = 0.7),
                   text = ~paste("ID:", NewID, "<br>Pop:", pop)) %>%
  add_trace(data = pca_complete1, x = ~Axis1, y = ~Axis2, z = ~Axis3, type = 'scatter3d', mode = 'text', text = ~NewID,
            textposition = "top center", textfont = list(color = '#000000')) %>%
  layout(title = "3D PCA Plot",
         scene = list(xaxis = list(title = paste0("PC1 (", round(explained_variance[1], 2), "%)")),
                      yaxis = list(title = paste0("PC2 (", round(explained_variance[2], 2), "%)")),
                      zaxis = list(title = paste0("PC3 (", round(explained_variance[3], 2), "%)")),
                      camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))))

plot_3d


# adult and larvae --------------------------------------------------------
#quick plot
pca = gl.pcoa(data_gl_filtered)
gl.pcoa.plot(glPca = pca, data_gl_filtered)

# PCA Analysis
pca_data <- tab(data_gl_filtered, freq = TRUE, NA.method = "mean") %>% na.omit() # Convert to tabular format and omit NAs
pca <- dudi.pca(pca_data, center = TRUE, scale = FALSE, nf = 2, scannf = FALSE) # Perform PCA
pca_complete2 <- data.frame(pca$li, pop = data_gl_filtered$pop) # Combine PCA results with population data
#use for adults
#pca_complete2 <- data.frame(pca$li) # Combine PCA results with population data

# Explained variance
(explained_variance <- pca$eig / sum(pca$eig) * 100)
scree_plot <- data.frame(PC = 1:length(explained_variance), Variance = explained_variance)

ggplot(scree_plot, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_line(aes(y = cumsum(Variance)), group = 1, color = "red") +
  geom_point(aes(y = cumsum(Variance)), color = "red") +
  labs(title = "Scree Plot", x = "Principal Component", y = "Percentage of Variance Explained") +
  theme_sleek2()

# Hopkins statistic
set.seed(123) # for reproducibility
(hopkins_stat <- hopkins(pca_data, n = nrow(pca_data) - 1))
# Calculated values 0-0.3 indicate regularly-spaced data. Values around 0.5 indicate random data. Values 0.7-1 indicate clustered data.
#all PD = 0.21

# K-means clustering
set.seed(123) # for reproducibility
kmeans_result <- kmeans(pca_data, centers = 3, nstart = 25)
individuals_in_cluster3 <- which(kmeans_result$cluster == 3) #find indiv in each cluster
silhouette_score <- silhouette(kmeans_result$cluster, dist(pca_data))
summary(silhouette_score)
plot(silhouette_score)
pca_complete2$Cluster <- as.factor(kmeans_result$cluster)
#PD: cluster 3 is quite strong, others poor to mod. 

pca_complete2 <- pca_complete2 %>%
  mutate(
    Stage = ifelse(str_detect(row.names(pca_complete2), "\\.a\\."), "Adu", "Lar"),
    MumID = str_extract(row.names(pca_complete2), "(?<=pd)\\d+"),
    RepID = str_extract(row.names(pca_complete2), "(?<=\\.)\\d+$"),
    NewID = paste0(Stage,  MumID, "_", RepID)
  )

my_palette <- c(
  "dodgerblue", "firebrick", "mediumseagreen", "orchid", "darkorange", "gold",
  "skyblue", "sandybrown", "palevioletred", "mediumturquoise", "khaki",
  "darkslategray", "plum", "lightslategray", "limegreen", "cornflowerblue",
  "tomato", 'red'
)

# Plot with ggrepel for label lines
t2 <- ggplot(pca_complete2, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(fill = pop, shape = Stage, color = ifelse(grepl("Lar", Stage), "red", "black")),
             size = 3, stroke = 1, alpha = 0.7, position = position_jitter(width = 0.1, height = 0.1)) +
  geom_text_repel(aes(label = NewID), size = 3, max.overlaps = 38, point.padding = 0.5, box.padding = 0.5) +
  stat_ellipse(aes(x = Axis1, y = Axis2, group = Cluster, color = Cluster), level = 0.95, linetype = 2, size = 1) + # Add ellipses around clusters
  scale_fill_manual(values = my_palette) +
  scale_color_manual(values = c("1" = "dodgerblue", "2" = "salmon", "3" = "mediumseagreen", "red" = "red", "black" = "black")) +
  scale_shape_manual(values = c("Adu" = 22, "Lar" = 21)) + # Set shapes: squares for adults and circles for larvae
  theme_sleek2() +
  labs(
    x = paste0("PCA1 (", round(explained_variance[1], 2), "%)"),
    y = paste0("PCA2 (", round(explained_variance[2], 2), "%)"),
    color = "Cluster", fill = "Population", shape = "Stage") # Add labels to the axes and legend
t2
#Note Larvae5_1 is likely from adult 13. Check notes. 




