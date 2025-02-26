










pca_data <- tab(data_genind_adult_unique, freq = TRUE, NA.method = "mean") %>% na.omit()

pca <- dudi.pca(pca_data, center = TRUE, scale = FALSE, nf = 2, scannf = FALSE) # Perform PCA
pca_complete <- data.frame(pca$li, pop = data_genind_adult_unique$pop)


(explained_variance <- pca$eig / sum(pca$eig) * 100)
scree_plot <- data.frame(PC = 1:length(explained_variance), Variance = explained_variance)

ggplot(scree_plot, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_line(aes(y = cumsum(Variance)), group = 1, color = "red") +
  geom_point(aes(y = cumsum(Variance)), color = "red") +
  labs(x = "Principal Component", y = "Percentage of Variance Explained") +
  theme_sleek2()

set.seed(123) # for reproducibility
(hopkins_stat <- hopkins(pca_data, n = nrow(pca_data) - 1))

kNNdistplot(pca_data, k = 5)  #k is no of nearest neighbors used, not clusters
elbow = 12.5 # Place this at the elbow of the line . Or 8 depending on run
abline(h = elbow, col = "red", lty = 2)  
dbscan_result <- dbscan(pca_data, eps = elbow, minPts = 5)
cluster_col_name <- paste0("cluster_dbscan_", elbow)
pca_complete[[cluster_col_name]] <- as.factor(dbscan_result$cluster)  #add cluster to pca_complete
ggplot(pca_complete, aes_string(x = "Axis1", y = "Axis2", color = cluster_col_name)) +
  geom_point(alpha = 0.6) +
  labs(title = paste("PCA Plot with DBSCAN clusters (eps =", elbow, ")"),
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal()
silhouette_score <- silhouette(dbscan_result$cluster, dist(pca_data))
plot(silhouette_score)  #see pca_complete to see which groups relate to what


set.seed(123) # for reproducibility
kmeans_result <- kmeans(pca_data, centers = 3, nstart = 25)
individuals_in_cluster3 <- which(kmeans_result$cluster == 3) #find indiv in each cluster
silhouette_score <- silhouette(kmeans_result$cluster, dist(pca_data))
summary(silhouette_score)
plot(silhouette_score)
pca_complete$cluster <- as.factor(kmeans_result$cluster)

pca_complete <- pca_complete %>%
  mutate(
    stage = ifelse(str_detect(row.names(pca_complete), "\\.a\\."), "Adult", "Larva"),  #add stage
    rep_id = str_extract(row.names(pca_complete), "(?<=\\.)\\d+$"),
    new_id = paste0(pop, "_", rep_id),
    id = rownames(pca_complete)
  )

ind_metrics <- data_genind_adult_unique@other$ind.metrics
ind_metrics_updated <- left_join(ind_metrics, pca_complete, by = 'id') %>%
  dplyr::select(-c(service, plate_location, stage.y))
data_genind_adult_unique@other$ind.metrics <- ind_metrics_updated
data_genind_adult_unique@other$ind.metrics <- data_genind_adult_unique@other$ind.metrics %>%
  mutate(cluster_colour = case_when(
    cluster == 1 ~ "mediumseagreen",
    cluster == 2 ~ "dodgerblue",
    cluster == 3 ~ "salmon"
  ))



clusters <- data_genind_adult@other$ind.metrics$cluster
data_genind_adult_subset1 <- data_genind_adult[clusters == "1", ]
data_genind_adult_subset2 <- data_genind_adult[clusters == "2", ]
data_genind_adult_subset3 <- data_genind_adult[clusters == "3", ]
genotype <- data_genind_adult@other$ind.metrics$genotype
data_genind_parents <- data_genind_adult[data_genind_adult@pop %in% c('9', '13', '14', '15'), ]
clusters <- data_genind_adult_unique@other$ind.metrics$cluster
data_genind_adult_subset1 <- data_genind_adult_unique[clusters == "1", ]
data_genind_adult_subset2 <- data_genind_adult_unique[clusters == "2", ]
data_genind_adult_subset3 <- data_genind_adult_unique[clusters == "3", ]

data1 <- dplyr::arrange(pca_complete, Axis1) # 
pca_complete <- pca_complete %>% mutate(across(c(stage, pop, rep_id, new_id), as.factor))
str(pca_complete)
my_palette <- c(
  "dodgerblue", "firebrick", "mediumseagreen", "orchid", "darkorange", "gold",
  "skyblue", "sandybrown", "palevioletred", "mediumturquoise", "khaki",
  "darkslategray", "plum", "lightslategray", "limegreen", "cornflowerblue",
  "tomato",'red'
)

t2 <- ggplot(pca_complete, aes(x=Axis1, y=Axis2, group=new_id)) +
  geom_point(aes(fill=pop), shape=21, size=4, stroke=1, alpha=0.8) +  # Points
  geom_text(aes(label=pop), vjust=1.5, hjust=0.5, color="black", size=3) +  # Add text labels
  scale_fill_manual(values=my_palette) +
  labs(x = "PCA1", y = "PCA2", color = "Population", fill = "Population") +
  theme_minimal() 
t2

t2 <- ggplot(pca_complete, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(color = factor(cluster)), shape = 22, 
             size = 3, stroke = 1, alpha = 0.7, position = position_jitter(width = 0.1, height = 0.1)
  ) +
  geom_text_repel(aes(label = new_id), size = 3, max.overlaps = 38, point.padding = 0.5, box.padding = 0.5, color = 'grey30') +
  scale_color_manual(values = c("3" = "salmon", "2" = "dodgerblue", "1" = "mediumseagreen")) +
  stat_ellipse(aes(x = Axis1, y = Axis2, group = cluster, color = cluster), level = 0.95, linetype = 2, size = 1) + # Add ellipses around clusters
  theme_sleek2() +
  labs(
    x = paste0("PCA1 (", round(explained_variance[1], 2), "%)"),
    y = paste0("PCA2 (", round(explained_variance[2], 2), "%)"),
    color = "cluster", fill = "Population", shape = "stage"
  ) 
t2






pca1 <- dudi.pca(pca_data, center = TRUE, scale = FALSE, nf = 3, scannf = FALSE) # Calculate the first three PCs

pca_complete1 <- data.frame(pca1$li, pop = data_gl_filtered$pop)
colnames(pca_complete1)[1:3] <- c("Axis1", "Axis2", "Axis3") # Rename columns to Axis1, Axis2, Axis3 for clarity

pca_complete1 <- pca_complete1 %>%
  mutate(
    stage = ifelse(str_detect(row.names(pca_complete1), "\\.a\\."), "Adult", "Larva"),
    mum_id = str_extract(row.names(pca_complete1), "(?<=pd)\\d+"),
    rep_id = str_extract(row.names(pca_complete1), "(?<=\\.)\\d+$"),
    new_id = paste0(stage, mum_id, "_", rep_id)
  ) %>%
  mutate(across(c(stage, mum_id, rep_id, new_id), as.factor))  # Convert relevant columns to factors

cluster_colors <- my_palette[pca_complete1$cluster]

plot_3d <- plot_ly(data = pca_complete1, x = ~Axis1, y = ~Axis2, z = ~Axis3, type = 'scatter3d', mode = 'markers',
                   marker = list(size = 5, color = cluster_colors, colorscale = list(c(0, 1), c(min(cluster_colors), max(cluster_colors))), opacity = 0.7),
                   text = ~paste("ID:", new_id, "<br>Pop:", pop)) %>%
  add_trace(data = pca_complete1, x = ~Axis1, y = ~Axis2, z = ~Axis3, type = 'scatter3d', mode = 'text', text = ~new_id,
            textposition = "top center", textfont = list(color = '#000000')) %>%
  layout(title = "3D PCA Plot",
         scene = list(xaxis = list(title = paste0("PC1 (", round(explained_variance[1], 2), "%)")),
                      yaxis = list(title = paste0("PC2 (", round(explained_variance[2], 2), "%)")),
                      zaxis = list(title = paste0("PC3 (", round(explained_variance[3], 2), "%)")),
                      camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))))

plot_3d



pca_data <- tab(data_gl_filtered, freq = TRUE, NA.method = "mean") %>% na.omit() # Convert to tabular format and omit NAs
pca <- dudi.pca(pca_data, center = TRUE, scale = FALSE, nf = 2, scannf = FALSE) # Perform PCA
pca_complete2 <- data.frame(pca$li, pop = data_gl_filtered$pop) # Combine PCA results with population data

(explained_variance <- pca$eig / sum(pca$eig) * 100)
scree_plot <- data.frame(PC = 1:length(explained_variance), Variance = explained_variance)

p0 = ggplot(scree_plot, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_line(aes(y = cumsum(Variance)), group = 1, color = "red") +
  geom_point(aes(y = cumsum(Variance)), color = "red") +
  labs(x = "Principal component axes", y = "Variance explained (%)") +
  theme_sleek1()+
  scale_y_continuous(limits = c(0, 100)) + 
  scale_x_continuous(limits = c(0, 30))
p0

set.seed(123) # for reproducibility
(hopkins_stat <- hopkins(pca_data, n = nrow(pca_data) - 1))

set.seed(123) # for reproducibility
kmeans_result <- kmeans(pca_data, centers = 3, nstart = 25)
individuals_in_cluster3 <- which(kmeans_result$cluster == 3) #find indiv in each cluster
silhouette_score <- silhouette(kmeans_result$cluster, dist(pca_data))
summary(silhouette_score)
plot(silhouette_score)
pca_complete2$cluster <- as.factor(kmeans_result$cluster)

pca_complete2 <- pca_complete2 %>%
  mutate(
    stage = ifelse(str_detect(row.names(pca_complete2), "\\.a\\."), "Adu", "Lar"),
    mum_id = pca_complete2$pop,
    rep_id = str_extract(row.names(pca_complete2), "(?<=\\.)\\d+$"),
    new_id = paste0(stage,  mum_id, "_", rep_id)
  )

my_palette <- c(
  "tomato", "mediumturquoise", "mediumseagreen", "orchid", "cornflowerblue", "gold",
  "skyblue", "darkorange", "palevioletred", "firebrick", "khaki",
  "darkslategray", "plum", "lightslategray", "red", "sandybrown",
  "dodgerblue", 'limegreen'
)
pca_complete2$pop_color <- my_palette[pca_complete2$pop]  # 

pca_complete2$color_by_stage <- ifelse(grepl("Lar", pca_complete2$stage), "indianred3", "grey60") # Assign colours based on the condition

t2 <- ggplot(pca_complete2, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(fill = pop_color, shape = stage, color = color_by_stage), 
             size = 3, stroke = 0.5, alpha = 0.7, position = position_jitter(width = 0.1, height = 0.1)) +
  geom_text_repel(aes(label = new_id, color = color_by_stage), 
                  size = 3, max.overlaps = 38, point.padding = 0.5, box.padding = 0.5, alpha = 1) +
  scale_fill_identity() +
  scale_shape_manual(values = c("Adu" = 22, "Lar" = 21)) +
  scale_color_identity() +
  theme_sleek2() +
  labs(
    x = paste0("PCA1 (", round(explained_variance[1], 2), "%)"),
    y = paste0("PCA2 (", round(explained_variance[2], 2), "%)"),
    color = "cluster", fill = "Population", shape = "stage")
t2





axis1_mid <- 0.5
axis2_mid <- 1

pca_complete2$quadrant <- with(pca_complete2, ifelse(Axis1 <= axis1_mid & Axis2 >= axis2_mid, "Top Left",
                                                     ifelse(Axis1 > axis1_mid & Axis2 >= axis2_mid, "Top Right",
                                                            ifelse(Axis1 <= axis1_mid & Axis2 < axis2_mid, "Bottom Left",
                                                                   "Bottom Right"))))

plots <- lapply(unique(pca_complete2$quadrant), function(q) {
  ggplot(pca_complete2[pca_complete2$quadrant == q,], aes(x = Axis1, y = Axis2)) +
    geom_point(aes(fill = pop_color, shape = stage, color = color_by_stage), 
               size = 4, stroke = 0.5, alpha = 0.7, position = position_jitter(width = 0.01, height = 0.01)) +
    geom_text_repel(aes(label = new_id, color = color_by_stage), 
                    size = 3, max.overlaps = 80, point.padding = 0.3, box.padding = 0.3, alpha = 1, segment.size = 0.5) +
    scale_fill_identity() +
    scale_color_identity() +
    scale_shape_manual(values = c("Adu" = 22, "Lar" = 21)) +
    theme_sleek2() +
    labs(x = "PCA1", y = "PCA2")  # Removed title
})
(t3 = plots[[3]])  #top right
(t4 = plots[[2]])  #bottom right



