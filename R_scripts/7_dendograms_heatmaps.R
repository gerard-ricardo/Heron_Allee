# dendogram and heatmaps ---------------------------------------------------------------

library(pheatmap)

# Compute the distance matrix (assuming this step is already done)
dist_matrix <- bitwise.dist(data_gl_filtered_adult)
dist_matrix_mat <- as.matrix(dist_matrix)

# Check and replace NA, NaN, and Inf values
dist_matrix_mat[is.na(dist_matrix_mat)] <- max(dist_matrix_mat, na.rm = TRUE)
dist_matrix_mat[is.nan(dist_matrix_mat)] <- max(dist_matrix_mat, na.rm = TRUE)
dist_matrix_mat[is.infinite(dist_matrix_mat)] <- max(dist_matrix_mat, na.rm = TRUE)

# Create a heatmap of the distance matrix
pheatmap(dist_matrix_mat, 
         cluster_rows = TRUE, # Cluster rows
         cluster_cols = TRUE, # Cluster columns
         display_numbers = FALSE, # Optionally display numbers
         main = "Distance Matrix Heatmap", # Title of the heatmap
         color = colorRampPalette(c("white", "blue", "red"))(50)) # Color scheme
#low values mean similarlity, so white boxes proaably indicate clones

# Perform Hierarchical Clustering
hc <- hclust(as.dist(dist_matrix_mat), method = "ward.D2")
# Plot the dendrogram
library(dendextend)

# Convert 'hc' to a dendrogram object
dend <- as.dendrogram(hc)
dend <- rotate(dend, 1:nleaves(dend))
# Set the labels for the dendrogram
dend <- set(dend, "labels", rownames(dist_matrix_mat))
plot(dend, horiz = TRUE, cex = 0.8)
mtext("Distance", side = 1, line = 2)

#Represents a distance matrix as a heatmap
D <- dist(as.matrix(data_gl_filtered),upper=TRUE,diag=TRUE)
gl.plot.heatmap(D)