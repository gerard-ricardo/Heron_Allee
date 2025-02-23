


dist_matrix <- bitwise.dist(data_gl_filtered_adult)
dist_matrix_mat <- as.matrix(dist_matrix)

dist_matrix_mat[is.na(dist_matrix_mat)] <- max(dist_matrix_mat, na.rm = TRUE)
dist_matrix_mat[is.nan(dist_matrix_mat)] <- max(dist_matrix_mat, na.rm = TRUE)
dist_matrix_mat[is.infinite(dist_matrix_mat)] <- max(dist_matrix_mat, na.rm = TRUE)

data_genind_adult@pop <- factor(rep("population1", nrow(data_genind_adult@tab))) #combine all
dist_matrix_mat <- gd.smouse(data_genind_adult, verbose = TRUE)  #Smouse and Peakall (1999) is a method used to quantify the
dist_matrix_mat[is.na(dist_matrix_mat)] <- max(dist_matrix_mat, na.rm = TRUE)
dist_matrix_mat[is.nan(dist_matrix_mat)] <- max(dist_matrix_mat, na.rm = TRUE)
dist_matrix_mat[is.infinite(dist_matrix_mat)] <- max(dist_matrix_mat, na.rm = TRUE)

pheatmap(dist_matrix_mat, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         display_numbers = FALSE, 
         main = "Distance Matrix Heatmap", 
         color = colorRampPalette(c("white", "blue", "red"))(50)) 

hc <- hclust(as.dist(dist_matrix_mat), method = "ward.D2")

dend <- as.dendrogram(hc)
dend <- rotate(dend, 1:nleaves(dend))
dend <- set(dend, "labels", rownames(dist_matrix_mat))
plot(dend, horiz = TRUE, cex = 0.8)
mtext("Distance", side = 1, line = 2)

D <- dist(as.matrix(data_gl_filtered),upper=TRUE,diag=TRUE)
gl.plot.heatmap(D)




genetic_distance_dist <- dist_matrix

phylo_tree <- nj(genetic_distance_dist) # Neighbor-joining method
plot(phylo_tree, main = "Phylogenetic Tree of SNP Data")
print(phylo_tree$edge.length)  #branch lengths


library(boot)
boot_phylo <- boot.phylo(phy = phylo_tree, x = as.matrix(genetic_distance_dist), FUN = function(x) nj(as.dist(x)))
plot(boot_phylo, main = "Bootstrapped Phylogenetic Tree")
