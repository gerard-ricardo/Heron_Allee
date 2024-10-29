# dendogram and heatmaps ---------------------------------------------------------------



# Compute the distance matrix (better for large data sets)
dist_matrix <- bitwise.dist(data_gl_filtered_adult)
dist_matrix_mat <- as.matrix(dist_matrix)

# Check and replace NA, NaN, and Inf values
dist_matrix_mat[is.na(dist_matrix_mat)] <- max(dist_matrix_mat, na.rm = TRUE)
dist_matrix_mat[is.nan(dist_matrix_mat)] <- max(dist_matrix_mat, na.rm = TRUE)
dist_matrix_mat[is.infinite(dist_matrix_mat)] <- max(dist_matrix_mat, na.rm = TRUE)

#try smouse (better for smaller micrsat data)
data_genind_adult@pop <- factor(rep("population1", nrow(data_genind_adult@tab))) #combine all
dist_matrix_mat <- gd.smouse(data_genind_adult, verbose = TRUE)  #Smouse and Peakall (1999) is a method used to quantify the
dist_matrix_mat[is.na(dist_matrix_mat)] <- max(dist_matrix_mat, na.rm = TRUE)
dist_matrix_mat[is.nan(dist_matrix_mat)] <- max(dist_matrix_mat, na.rm = TRUE)
dist_matrix_mat[is.infinite(dist_matrix_mat)] <- max(dist_matrix_mat, na.rm = TRUE)

# Create a heatmap of the distance matrix
pheatmap(dist_matrix_mat, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         display_numbers = FALSE, 
         main = "Distance Matrix Heatmap", 
         color = colorRampPalette(c("white", "blue", "red"))(50)) 
#low values mean similarlity, so white boxes proaably indicate clones

# Perform Hierarchical Clustering
hc <- hclust(as.dist(dist_matrix_mat), method = "ward.D2")
# Plot the dendrogram

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



# phylo tree --------------------------------------------------------------
#evolutionary relationships among a set of organisms or genes

# Alternatively, you can construct a phylogenetic tree
genetic_distance_dist <- dist_matrix

phylo_tree <- nj(genetic_distance_dist) # Neighbor-joining method
plot(phylo_tree, main = "Phylogenetic Tree of SNP Data")
print(phylo_tree$edge.length)  #branch lengths


# You may also want to bootstrap the phylogenetic tree to assess the robustness of the branches
library(boot)
boot_phylo <- boot.phylo(phy = phylo_tree, x = as.matrix(genetic_distance_dist), FUN = function(x) nj(as.dist(x)))
plot(boot_phylo, main = "Bootstrapped Phylogenetic Tree")
