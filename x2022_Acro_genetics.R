## acro orpheus

# load libraries ----------------------------------------------------------
#install.packages('dartr')
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

# data_gl <- gl.read.dart(filename = "./data/Report_DAc23-7804_SNP_2 - Copy.csv", ind.metafile = "./data/meta_acro_ordered.csv", topskip = 6)
# data_gl <- gl.reassign.pop(data_gl, as.pop = "stage")
# # #recalculate metrics
# data_gl <- gl.recalc.metrics(data_gl, v = 3) # recalculate loci metrics
# save(data_gl, file = file.path("./Rdata", "2022_acro_gl.RData"))
load("./Rdata/2022_acro_gl.RData")  #data_gl


#Extrct tenuis from all
data_gl$ind.names
at_indices <- grepl("^at", data_gl$ind.names) # "^at" indicates strings starting with 'at'
data_gl <- data_gl[at_indices, ]

# calculate coverage metrics - mean number of reads that cover reference (30 good). Inc depth/reads will beter this.
summary(data_gl$other$loc.metrics$coverage)
data_gl$other$loc.metrics
data_gl$other$loc.metrics$coverage <- data_gl$other$loc.metrics$AvgCountRef + data_gl$other$loc.metrics$AvgCountSnp
median(data_gl$other$loc.metrics$coverage) #  28.6
min((data_gl$other$loc.metrics$coverage)) # 5.375
max((data_gl$other$loc.metrics$coverage)) #974.6452
sd(data_gl$other$loc.metrics$coverage) / sqrt(1996) #1.858626
hist(data_gl$other$loc.metrics$coverage)

# data filtering ----------------------------------------------------------

data_gl_filtered <- data_gl
#ind = 45, loc = 1996
#note already low loci

##secondaries
gl.report.secondaries(data_gl_filtered)
data_gl_filtered <- gl.filter.secondaries(data_gl_filtered, method="random", verbose = 3) #remove loci fragment that shared SNPs. Only keep 1
#ind = 45, loc = 1607

#rdepth
gl.report.rdepth(data_gl_filtered)
#med =  28.6  Generally 10 is considered min
data_gl_filtered <- gl.filter.rdepth(data_gl_filtered,  lower = 10, v = 3) # filter by loci callrate
# 45 ind, 1003    loc

##reproducibility 
gl.report.reproducibility(data_gl_filtered )
data_gl_filtered <- gl.filter.reproducibility(data_gl_filtered, t=0.95, v=3) #filter out loci with limited reproducibility
#Platy at 95%: ind = 45, loci = 95 
#this removed alot

# callrate loci (non missing data)
gl.report.callrate(data_gl_filtered, method = "loc") 
#PLaty = 40%
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "loc", threshold = 0.7, v = 3) # filter by loci callrate
##At 70%: ind = 45, loci = 101

#Minor Allele Frequency (MAF) and Coverage Filter:
list.match <- data_gl_filtered$loc.names[
  which(data_gl_filtered$other$loc.metrics$OneRatioSnp > 0.01 & 
          data_gl_filtered$other$loc.metrics$OneRatioSnp < 0.99 & 
          data_gl_filtered$other$loc.metrics$OneRatioRef < 0.99 & 
          data_gl_filtered$other$loc.metrics$OneRatioRef > 0.01 & 
          data_gl_filtered$other$loc.metrics$coverage > 4)
]
data_gl_filtered <- data_gl_filtered[, match(list.match, data_gl_filtered$loc.names)]
#ind = 45, loci = 95

#filter monomorphs
data_gl_filtered <- gl.filter.monomorphs(data_gl_filtered, v=3) #remove monomorphic loci (loci with 1 fixed allele across the entire dataset (no differences) )
#ind = 45, loci = 95


## call rate ind (non missing data). low could indicate poor extract or reference genome or contamination.
#individuals
gl.report.callrate(data_gl_filtered, method = "ind") 
#~86%. Lowest is 0.67
pre_filt_ind <- data_gl_filtered@ind.names
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "ind", threshold = 0.58, v = 3) # filter by ind callrate
filt_ind <- data_gl_filtered@ind.names
(lost_ind <- setdiff(pre_filt_ind, filt_ind))
#none lost  #ind = 45, loci = 95

data_gl_filtered <- gl.recalc.metrics(data_gl_filtered, v = 3) # recalculate loci metrics

# population filtering and objects ----------------------------------------

# look into genotype as population
data_gl_filtered <- gl.reassign.pop(data_gl_filtered, as.pop = "genotype")
data_gl_filtered

# Convert GENIND OBJECT all indiv
data_genind <- gl2gi(data_gl_filtered)
#genind object are 2-col (ref/var) loci format, where counts =  numbers of each allele i.e 2/0 means two reference. 

# Filter out eggs and larvae to keep only adults
adults_indices <- which(data_gl_filtered@other$ind.metrics$stage == "adult")
data_gl_filtered_adult <- data_gl_filtered[adults_indices, ]
data_gl_filtered_adult@other$ind.metrics$stage <- droplevels(data_gl_filtered_adult@other$ind.metrics$stage)

# Convert genind adults only
data_genind_adult <- gl2gi(data_gl_filtered_adult)

#create 0_1 coded df
mat_0_1_2_coded = data_genind_adult$tab
mat_0_1_2_coded_char <- as.character(mat_0_1_2_coded)
mat_0_1_2_coded_char[grepl("^2$", mat_0_1_2_coded_char)] <- "1"
mat_0_1_coded <- matrix(as.numeric(mat_0_1_2_coded_char), nrow = nrow(mat_0_1_2_coded), ncol = ncol(mat_0_1_2_coded))


# PCA ---------------------------------------------------------------------
#   adult only ---------------------------------------------------------------------

#quick plot
pca = gl.pcoa(data_gl_filtered_adult)
gl.pcoa.plot(glPca = pca, data_gl_filtered_adult)

# PCA Analysis
#pca_data <- tab(data_gl_filtered_adult, freq = TRUE, NA.method = "mean") %>% na.omit() # Convert to tabular format and omit NAs
pca_data <- tab(data_genind_adult, freq = TRUE, NA.method = "mean") %>% na.omit()

pca <- dudi.pca(pca_data, center = TRUE, scale = FALSE, nf = 2, scannf = FALSE) # Perform PCA
#pca_complete <- data.frame(pca$li, pop = data_gl_filtered_adult$pop) # Combine PCA results with population data
pca_complete <- data.frame(pca$li, pop = data_genind_adult$pop)

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
pca_complete$cluster <- as.factor(kmeans_result$cluster)
#PD: cluster 3 is quite strong, others poor to mod. 

# DBSCAN clustering
# Find the appropriate eps value using kNNdistplot
kNNdistplot(pca_data, k = 5)
elbow = 0 # Place this at the elbow of the line
abline(h = elbow, col = "red", lty = 2)  
library(dbscan)
# Function to perform DBSCAN clustering and plot results
perform_dbscan <- function(pca_data, pca_complete, eps_value, min_pts = 5) {
  dbscan_result <- dbscan(pca_data, eps = eps_value, minPts = min_pts)
  cluster_col_name <- paste0("cluster_dbscan_", eps_value)
  pca_complete[[cluster_col_name]] <- as.factor(dbscan_result$cluster)
  plot <- ggplot(pca_complete, aes_string(x = "Axis1", y = "Axis2", color = cluster_col_name)) +
    geom_point(alpha = 0.6) +
    labs(title = paste("PCA Plot with DBSCAN clusters (eps =", eps_value, ")"),
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
#No clustering


# plotting
#PD
pca_complete <- pca_complete %>%
  mutate(
    stage = ifelse(str_detect(row.names(pca_complete), "\\.a$"), "Adult", "Larva"),  #add stage
    id = rownames(pca_complete),
    new_id = id
  )

#add cluster to meta data of objects
data_gl_filtered_adult@other$ind.metrics = left_join(data_gl_filtered_adult@other$ind.metrics, pca_complete, by  = 'id') %>% 
  dplyr::select(-c(service, plate_location, stage.y)) 
ind_metrics <- data_genind_adult@other$ind.metrics
ind_metrics_updated <- left_join(ind_metrics, pca_complete, by = 'id') %>%
  dplyr::select(-c(service, plate_location, stage.y))
data_genind_adult@other$ind.metrics <- ind_metrics_updated
# subset by group
clusters <- data_genind_adult@other$ind.metrics$cluster
data_genind_adult_subset1 <- data_genind_adult[clusters == "1", ]
data_genind_adult_subset2 <- data_genind_adult[clusters == "2", ]
data_genind_adult_subset3 <- data_genind_adult[clusters == "3", ]


data1 <- dplyr::arrange(pca_complete, Axis1) # 
pca_complete <- pca_complete %>% mutate(across(c(stage, pop), as.factor))
#pca_complete1 = pca_complete %>% select(cluster, pop) %>% rename(id = pop)  #for bathymetry file 
str(pca_complete)
my_palette <- c(
  "dodgerblue", "firebrick", "mediumseagreen", "orchid", "darkorange", "gold",
  "skyblue", "sandybrown", "palevioletred", "mediumturquoise", "khaki",
  "darkslategray", "plum", "lightslategray", "limegreen", "cornflowerblue",
  "tomato",'red'
)

#color individuals
t2 <- ggplot(pca_complete, aes(x=Axis1, y=Axis2, group=new_id)) +
  geom_point(aes(fill=pop), shape=21, size=4, stroke=1, alpha=0.8) +  # Points
  geom_text(aes(label=pop), vjust=1.5, hjust=0.5, color="black", size=3) +  # Add text labels
  scale_fill_manual(values=my_palette) +
  labs(x = "PCA1", y = "PCA2", color = "Population", fill = "Population") +
  theme_minimal() 
t2

#per cluster
t2 <- ggplot(pca_complete, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(color = factor(cluster)), shape = 22, 
             size = 3, stroke = 1, alpha = 0.7, position = position_jitter(width = 0.1, height = 0.1)
  ) +
  geom_text_repel(aes(label = new_id), size = 3, max.overlaps = 38, point.padding = 0.5, box.padding = 0.5) +
  scale_color_manual(values = c("1" = "dodgerblue", "2" = "salmon", "3" = "mediumseagreen")) +
  stat_ellipse(aes(x = Axis1, y = Axis2, group = cluster, color = cluster), level = 0.95, linetype = 2, size = 1) + # Add ellipses around clusters
  theme_sleek2() +
  labs(
    x = paste0("PCA1 (", round(explained_variance[1], 2), "%)"),
    y = paste0("PCA2 (", round(explained_variance[2], 2), "%)"),
    color = "cluster", fill = "Population", shape = "stage"
  ) 
t2
#ggsave(t2, filename = 'heron_pca_clusters.tiff',  path = "./plots", device = "tiff",  width = 5, height = 5)  #this often works better than pdf


# Convert the ggplot to an interactive plotly plot
t2_interactive <- ggplotly(t2)
t2_interactive



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
pca_complete2$cluster <- as.factor(kmeans_result$cluster)
#PD: cluster 3 is quite strong, others poor to mod. 

pca_complete2 <- pca_complete2 %>%
  mutate(
    stage = ifelse(str_detect(row.names(pca_complete2), "\\.a$"), "Adu", "Lar"),
    mum_id = str_extract(row.names(pca_complete2), "(?<=\\.)\\d+$"),
    geno = str_extract(row.names(pca_complete2), "(?<=at)\\d+"),
    id = ifelse(!is.na(mum_id), mum_id, geno),
    new_id = paste0(stage,  id)
  )

my_palette <- c(
  "dodgerblue", "firebrick", "mediumseagreen", "orchid", "darkorange", "gold",
  "skyblue", "sandybrown", "palevioletred", "mediumturquoise", "khaki",
  "darkslategray", "plum", "lightslategray", "limegreen", "cornflowerblue",
  "tomato", 'red'
)

# Plot with ggrepel for label lines
t2 <- ggplot(pca_complete2, aes(x = Axis1, y = Axis2)) +
  geom_point(aes(fill = pop, shape = stage, color = ifelse(grepl("Lar", stage), "red", "black")),
             size = 3, stroke = 1, alpha = 0.7, position = position_jitter(width = 0.1, height = 0.1)) +
  geom_text_repel(aes(label = new_id), size = 3, max.overlaps = 38, point.padding = 0.5, box.padding = 0.5) +
  stat_ellipse(aes(x = Axis1, y = Axis2, group = cluster, color = cluster), level = 0.95, linetype = 2, size = 1) + # Add ellipses around clusters
  scale_fill_manual(values = my_palette) +
  scale_color_manual(values = c("1" = "dodgerblue", "2" = "salmon", "3" = "mediumseagreen", "red" = "red", "black" = "black")) +
  scale_shape_manual(values = c("Adu" = 22, "Lar" = 21)) + # Set shapes: squares for adults and circles for larvae
  theme_sleek2() +
  labs(
    x = paste0("PCA1 (", round(explained_variance[1], 2), "%)"),
    y = paste0("PCA2 (", round(explained_variance[2], 2), "%)"),
    color = "cluster", fill = "Population", shape = "stage") # Add labels to the axes and legend
t2
#


# CERVUS proceessing ------------------------------------------------------

# Extracting SNP data list from genind object
data1 <- data_genind@tab  # This assumes SNP data is stored in the 'tab' slot of the genind object
head(data1)

data1 = t(data1) %>% data.frame()
#data1$SNP = data_genind@other$loc.metrics$SNP
colnames(data1) <- gsub("\\.", "_", colnames(data1))
data1$rownames <- rownames(data1)
data1 <- data1 %>% mutate(CloneID = sub("-.*", "", rownames), allele = sub(".*\\.(.)$", "\\1", rownames)) %>% 
  select(CloneID, allele, everything()) %>% select(-rownames) %>% arrange(., CloneID)
#data1 %>%  dplyr::select(., CloneID, pd13_l_14_11, pd13_a_1) %>% filter(., CloneID == check )  #

data1_long = data1 %>% tidyr::pivot_longer(-c(CloneID, allele) ,  names_to = "id" ,values_to = "counts") %>% arrange(., CloneID, id) %>% data.frame()  #keep vec.x, add all other columns to factors , add all their values to meas)


nrow(data1_long)  #5648
row_counts = data1_long %>%
  group_by(CloneID, id) %>%
  summarise(row_count = n()) %>%
  ungroup()
min(row_counts$row_count)

# Filter out groups that do not have exactly 2 rows
data1_long = data1_long %>%
  group_by(CloneID, id) %>%
  filter(n() == 2) %>%
  ungroup() %>%
  data.frame()
nrow(data1_long)  #4576
any(data1_long == "", na.rm = TRUE) # This will return TRUE if there are any empty strings in the dataframe

#data1_long = data1_long[129:140,]

# Function to repeat allele with a separator
repeat_with_separator <- function(allele, counts, sep = "_") {
  if (is.na(allele) | is.na(counts)) {
    return("NA")
  } else {
    return(paste(rep(allele, counts), collapse = sep))
  }
}

# Apply the function to the dataframe
data1_long <- data1_long %>%
  mutate(allele_multiplied = mapply(repeat_with_separator, allele, counts, sep = ",")) # you can change sep to "_" if you want underscores
any(data1_long == "", na.rm = TRUE) # This will return TRUE if there are any empty strings in the dataframe

# Display the intermediate result
print(data1_long)

# Group and summarise the data
summary_data <- data1_long %>%
  group_by(CloneID, id) %>%
  summarise(
    id = first(id),
    base = paste(allele_multiplied, collapse = ",")
  ) %>%
  arrange(id) %>%
  data.frame()

# Remove leading or trailing commas in summary_data$base
summary_data$base <- gsub("^,+|,+$", "", summary_data$base)


# Separate the 'base' column into 'base1' and 'base2'
summary_data <- summary_data %>%
  separate(base, into = c("a", "b"), sep = ",", extra = "drop", fill = "right")

data1_long = summary_data %>% tidyr::pivot_longer(-c(CloneID, id) ,  names_to = "base" ,values_to = "code") %>% arrange(., CloneID, id) %>% data.frame()  #keep vec.x, add all other columns to factors , add all their values to meas)
data1_long$LocusID = paste0(data1_long$CloneID, data1_long$base)
data1_long = data1_long %>% select(-c(CloneID , base))
data1_long <- data1_long %>%
  mutate(code = ifelse(code == "NA", "*", code))
data_wide <- data1_long %>% tidyr::pivot_wider(names_from = LocusID, values_from = code, names_prefix = "X") %>% 
  data.frame()#year goes to columns, their areas go as the values, area is the prefix

#new
write.csv(data_wide, row.names = FALSE,
          file = file.path("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects/3_Heron_Platy_ms/Cervus",
                           "acro_ten_map_letters_code2_2.csv"))

### Offspring file
indices_with_l <- grep("l", data_wide$id)
labels_with_l <- data_wide$id[indices_with_l]
df1 <- data.frame(offspring = labels_with_l)
dam <- sub("_.*", "", df1$offspring)
#dam_counts <- table(df1$dam)
#known_dam <- rep(names(dam_counts), dam_counts)
# known_dam <- gsub("_(\\d)$", "_0\\1", known_dam)
#df1$known_dam <- known_dam
# df1$Offspring <- gsub("pd(\\d)_", "pd0\\1_", df1$Offspring)
# df1$known_dam <- gsub("pd(\\d)_", "pd0\\1_", df1$known_dam)

indices_with_a <- grep("_a", data_wide$id)
labels_with_a <- data_wide$id[indices_with_a]
first_matching_label <- sapply(dam, function(dam) {
  matching_labels <- grep(dam, labels_with_a, value = TRUE)
  if (length(matching_labels) > 0) {
    return(matching_labels[1])
  } else {
    return(NA)
  }
})
#df1$known_dam  = first_matching_label
#find_error_dam = which(df1$known_dam == 'pd5_a_1')
#df1[find_error_dam, 2] = 'pd13_a_1'
len = length(unique(labels_with_a))
not_known_dam <- labels_with_a
length(not_known_dam)
cands <- rep(labels_with_a, length(df1$offspring)) %>% sort(.)
cands_df <- matrix(cands, nrow = length(df1$offspring), ncol = len, byrow = FALSE) %>% data.frame()
colnames(cands_df) <- rep("candidate", len)
offspring_df <- cbind(df1, cands_df)
offspring_df

write.csv(offspring_df, row.names = FALSE, file = file.path("C:/Users/gerar/OneDrive/1_Work/4_Writing/1_Allee_effects/3_Heron_Platy_ms/Cervus", 
                                                            "offspring_acro_ten.csv"))





