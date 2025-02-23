
library(dartR)
library(dartR.popgen)
library(PopGenReport)
library(adegenet)
library(tictoc)
library(HardyWeinberg)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(stringr)
library(ggrepel)
library(hierfstat)
library(ape)
library(poppr)
library(pegas)
library(dbscan)
library(sp)
library(rgdal)
library(clustertend)
library(cluster)
library(plotly)
library(pheatmap)
library(dendextend)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2") # set theme in code

load("./Rdata/2022_acro_gl.RData")  #data_gl


data_gl$ind.names
sp_indices <- grepl("^sp", data_gl$ind.names) # "^at" indicates strings starting with 'at'
data_gl <- data_gl[sp_indices, ]

summary(data_gl$other$loc.metrics$coverage)
data_gl$other$loc.metrics
data_gl$other$loc.metrics$coverage <- data_gl$other$loc.metrics$AvgCountRef + data_gl$other$loc.metrics$AvgCountSnp
median(data_gl$other$loc.metrics$coverage) #  28.6
min((data_gl$other$loc.metrics$coverage)) # 5.375
max((data_gl$other$loc.metrics$coverage)) #974.6452
sd(data_gl$other$loc.metrics$coverage) / sqrt(1996) #1.858626
hist(data_gl$other$loc.metrics$coverage)


(data_gl_filtered <- data_gl)

gl.report.secondaries(data_gl_filtered)
(data_gl_filtered <- gl.filter.secondaries(data_gl_filtered, method="random", verbose = 3)) #remove loci fragment that shared SNPs. Only keep 1

gl.report.rdepth(data_gl_filtered)
data_gl_filtered <- gl.filter.rdepth(data_gl_filtered,  lower = 10, v = 3) # filter by loci callrate

gl.report.reproducibility(data_gl_filtered )
data_gl_filtered <- gl.filter.reproducibility(data_gl_filtered, t=0.90, v=3) #filter out loci with limited reproducibility

gl.report.callrate(data_gl_filtered, method = "loc") 
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "loc", threshold = 0.4, v = 3) # filter by loci callrate

list.match <- data_gl_filtered$loc.names[
  which(data_gl_filtered$other$loc.metrics$OneRatioSnp > 0.01 & 
          data_gl_filtered$other$loc.metrics$OneRatioSnp < 0.99 & 
          data_gl_filtered$other$loc.metrics$OneRatioRef < 0.99 & 
          data_gl_filtered$other$loc.metrics$OneRatioRef > 0.01 & 
          data_gl_filtered$other$loc.metrics$coverage > 4)
]
data_gl_filtered <- data_gl_filtered[, match(list.match, data_gl_filtered$loc.names)]

data_gl_filtered <- gl.filter.monomorphs(data_gl_filtered, v=3) #remove monomorphic loci (loci with 1 fixed allele across the entire dataset (no differences) )


gl.report.callrate(data_gl_filtered, method = "ind") 
pre_filt_ind <- data_gl_filtered@ind.names
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "ind", threshold = 0.6, v = 3) # filter by ind callrate
filt_ind <- data_gl_filtered@ind.names
(lost_ind <- setdiff(pre_filt_ind, filt_ind))

data_gl_filtered <- gl.recalc.metrics(data_gl_filtered, v = 3) # recalculate loci metrics


data_gl_filtered <- gl.reassign.pop(data_gl_filtered, as.pop = "genotype")
data_gl_filtered

data_genind <- gl2gi(data_gl_filtered)

adults_indices <- which(data_gl_filtered@other$ind.metrics$stage == "adults")
data_gl_filtered_adult <- data_gl_filtered[adults_indices, ]
data_gl_filtered_adult@other$ind.metrics$stage <- droplevels(data_gl_filtered_adult@other$ind.metrics$stage)
data_gl_filtered_adult_unique = data_gl_filtered[which(data_gl_filtered@other$ind.metrics$rep == "1"),]

data_genind_adult <- gl2gi(data_gl_filtered_adult)
data_genind_ad_unique <- gl2gi(data_gl_filtered_adult_unique)

mat_0_1_2_coded = data_genind_adult$tab
mat_0_1_2_coded_char <- as.character(mat_0_1_2_coded)
mat_0_1_2_coded_char[grepl("^2$", mat_0_1_2_coded_char)] <- "1"
mat_0_1_coded <- matrix(as.numeric(mat_0_1_2_coded_char), nrow = nrow(mat_0_1_2_coded), ncol = ncol(mat_0_1_2_coded))




pca_data <- tab(data_genind_ad_unique, freq = TRUE, NA.method = "mean") %>% na.omit()
pca <- dudi.pca(pca_data, center = TRUE, scale = FALSE, nf = 2, scannf = FALSE) # Perform PCA
pca_complete <- data.frame(pca$li, pop = data_genind_ad_unique$pop)

(explained_variance <- pca$eig / sum(pca$eig) * 100)
scree_plot <- data.frame(PC = 1:length(explained_variance), Variance = explained_variance)

ggplot(scree_plot, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_line(aes(y = cumsum(Variance)), group = 1, color = "red") +
  geom_point(aes(y = cumsum(Variance)), color = "red") +
  labs(title = "Scree Plot", x = "Principal Component", y = "Percentage of Variance Explained") +
  theme_sleek2()

set.seed(123) # for reproducibility
(hopkins_stat <- hopkins(pca_data, n = nrow(pca_data) - 1))

set.seed(123) # for reproducibility
kmeans_result <- kmeans(pca_data, centers = 2, nstart = 25)  #no clear k, but probbabyl 2 the closest considering low sample size.
individuals_in_cluster <- which(kmeans_result$cluster == 1) #find indiv in each cluster
silhouette_score <- silhouette(kmeans_result$cluster, dist(pca_data))
summary(silhouette_score)
nrow(pca_data)
(wcss_1 <- sum(kmeans(pca_data, centers = 1, nstart = 25)$tot.withinss))
(wcss_2 <- sum(kmeans(pca_data, centers = 2, nstart = 25)$tot.withinss))
(wcss_3 <- sum(kmeans(pca_data, centers = 3, nstart = 25)$tot.withinss))
(wcss_4 <- sum(kmeans(pca_data, centers = 4, nstart = 25)$tot.withinss))


pca_complete$kmeans_cluster <- as.factor(kmeans_result$cluster) #add to pca data



kNNdistplot(pca_data, k = 5)  #k is no of nearest neighbours used
elbow = 12.8 # Place this at the elbow of the line
eps_value <- elbow 
abline(h = elbow, col = "red", lty = 2)  

dbscan_result <- dbscan(pca_data, eps = eps_value, minPts = 5)
cluster_col_name <- paste0("cluster_dbscan")
pca_complete[[cluster_col_name]] <- as.factor(dbscan_result$cluster)  #add dbscan to pca data
plot <- ggplot(pca_complete, aes_string(x = "Axis1", y = "Axis2", color = cluster_col_name)) +
    geom_point(alpha = 0.6) +
    labs(title = paste("PCA Plot with DBSCAN clusters (eps =", eps_value, ")"),
         x = "Principal Component 1",
         y = "Principal Component 2") +
    theme_minimal()
plot
silhouette_score <- silhouette(dbscan_result$cluster, dist(pca_data))
(dbscan_result)
(summary(silhouette_score))

pca_complete <- pca_complete %>%
  mutate(
    stage = ifelse(str_detect(row.names(pca_complete), "\\.a$"), "Adult", "Larva"),  #add stage
    id = rownames(pca_complete),
    new_id = id,
    kmeans_clust = kmeans_cluster,
    dbscan_clust = cluster_dbscan
  )



data1 <- dplyr::arrange(pca_complete, Axis1) # 
pca_complete <- pca_complete %>% mutate(across(c(stage, pop), as.factor))
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
  geom_point(aes(color = factor(dbscan_clust )), shape = 22, 
             size = 3, stroke = 1, alpha = 0.7, position = position_jitter(width = 0.1, height = 0.1)
  ) +
  geom_text_repel(aes(label = new_id), size = 3, max.overlaps = 38, point.padding = 0.5, box.padding = 0.5) +
  scale_color_manual(values = c("1" = "dodgerblue", "2" = "salmon", "3" = "mediumseagreen")) +
  theme_sleek2() +
  labs(
    x = paste0("PCA1 (", round(explained_variance[1], 2), "%)"),
    y = paste0("PCA2 (", round(explained_variance[2], 2), "%)"),
    color = "dbscan_clust ", fill = "Population", shape = "stage"
  ) 
t2






lar_indices <- which(data_gl_filtered@other$ind.metrics$id != c("at.l.2.14", "at.l.2.16"))
data_gl_filtered <- data_gl_filtered[lar_indices, ]

pca = gl.pcoa(data_gl_filtered)
gl.pcoa.plot(glPca = pca, data_gl_filtered)

pca_data <- tab(data_gl_filtered, freq = TRUE, NA.method = "mean") %>% na.omit() # Convert to tabular format and omit NAs
pca <- dudi.pca(pca_data, center = TRUE, scale = FALSE, nf = 2, scannf = FALSE) # Perform PCA
pca_complete2 <- data.frame(pca$li, pop = data_gl_filtered$pop) # Combine PCA results with population data

(explained_variance <- pca$eig / sum(pca$eig) * 100)
scree_plot <- data.frame(PC = 1:length(explained_variance), Variance = explained_variance)

ggplot(scree_plot, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_line(aes(y = cumsum(Variance)), group = 1, color = "red") +
  geom_point(aes(y = cumsum(Variance)), color = "red") +
  labs(title = "Scree Plot", x = "Principal Component", y = "Percentage of Variance Explained") +
  theme_sleek2()

set.seed(123) # for reproducibility
(hopkins_stat <- hopkins(pca_data, n = nrow(pca_data) - 1))

set.seed(123) # for reproducibility
kmeans_result <- kmeans(pca_data, centers = 1, nstart = 25)
individuals_in_cluster3 <- which(kmeans_result$cluster == 3) #find indiv in each cluster
silhouette_score <- silhouette(kmeans_result$cluster, dist(pca_data))
summary(silhouette_score)
pca_complete2$cluster <- as.factor(kmeans_result$cluster)

pca_complete2 <- pca_complete2 %>%
  mutate(
    stage = ifelse(str_detect(row.names(pca_complete2), "\\.a$"), "Adu", "Lar"),
    mum_id = ifelse(stage == "Lar", as.character(cumsum(stage == "Lar")), 
                    str_extract(row.names(pca_complete2), "(?<=\\.)\\d+$")),
    geno = str_extract(row.names(pca_complete2), "(?<=at)\\d+"),
    id = ifelse(!is.na(mum_id), mum_id, geno),
    new_id = paste0(stage, id)
  )

my_palette <- c(
  "dodgerblue", "firebrick", "mediumseagreen", "orchid", "darkorange", "gold",
  "skyblue", "sandybrown", "palevioletred", "mediumturquoise", "khaki",
  "darkslategray", "plum", "lightslategray", "limegreen", "cornflowerblue",
  "tomato", 'red'
)

t2 <- ggplot(pca_complete2, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(fill = pop, shape = stage, color = ifelse(grepl("Lar", stage), "red", "black")),
             size = 3, stroke = 1, alpha = 0.7, position = position_jitter(width = 0.1, height = 0.1)) +
  geom_text_repel(aes(label = new_id), size = 3, max.overlaps = 38, point.padding = 0.5, box.padding = 0.5) +
  scale_fill_manual(values = my_palette) +
  scale_color_manual(values = c("1" = "dodgerblue", "2" = "salmon", "3" = "mediumseagreen", "red" = "red", "black" = "black")) +
  scale_shape_manual(values = c("Adu" = 22, "Lar" = 21)) + # Set shapes: squares for adults and circles for larvae
  theme_sleek2() +
  labs(
    x = paste0("PCA1 (", round(explained_variance[1], 2), "%)"),
    y = paste0("PCA2 (", round(explained_variance[2], 2), "%)"),
    color = "cluster", fill = "Population", shape = "stage") # Add labels to the axes and legend
t2






