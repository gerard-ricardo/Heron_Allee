
load('./Rdata/2022 spawning times.RData') #data1
load("./Rdata/2022_Heron_dist_mat.RData") #lst

str(data1) #check data type is correct
data1 = data1[complete.cases(data1), ]  #make sure import matches NA type
data1

library(ggplot2)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2")  #set theme in code


data2 = split(data1, data1$date)

data3.1 = data2$`13_12`
time_mat <- dist(data3.1[,3], diag = T)
time_mat

data4.1 = data2$`14_12`
time_mat1 <- dist(data4.1[,3], diag = T)
time_mat1

dist_mat = lst$dist_mat
dist_mat1 = lst$dist_mat1

library(vegan)
(mantel_result <- mantel(time_mat, dist_mat))

(mantel_result1 <- mantel(time_mat1, dist_mat1))




meta2 = data_gl_filtered_adult@other$ind.metrics
meta2 <- na.omit(meta2, subset = "spawn_time") # Remove rows where spawn_time is NA
meta2$spawn_time_clean = gsub("[^0-9]", "", meta2$spawn_time) %>% as.numeric()# Replace non-numerals with empty string
str(meta2)
meta2$spawn_time_cat = ifelse(meta2$spawn_time_clean <45, 'early', 'late')
meta2$dummy = ifelse(meta2$spawn_time_cat == 'early', 0, 1)
clus_spawn = meta2 %>%  dplyr::select(id, cluster, spawn_time_clean, spawn_time_cat, rep_id, dummy)
clus_spawn = subset(clus_spawn, rep_id == 2)  #used rep '2' because it had an extra
clus_spawn = subset(clus_spawn, id != 'pd12.a.2')  #remove 12. 
clus_spawn$cluster <- droplevels(clus_spawn$cluster)   #drop unused levels




(table <- table(clus_spawn$cluster, clus_spawn$spawn_time_cat))  # Create a contingency table
clus_spawn[which(clus_spawn$cluster == '3' & clus_spawn$spawn_time_cat == 'late'),]
chisq.test(table)  # Perform the chi-squared test
(fisher.test(table))  # Perform Fisher's Exact test
ggplot(clus_spawn, aes(x = spawn_time_clean, fill = factor(cluster))) +  # Set spawn_time_clean as x, and fill colour by cluster
  geom_density(alpha = 0.5) +  # Adjust transparency with alpha if needed
  labs(x = "Spawn Time Clean", y = "Density", fill = "Cluster") +  # Label the axes and legend
  ggtitle("Density Plot of Spawn Time by Cluster") +  # Add a title
  theme_minimal()  # Use a minimal theme


str(clus_spawn)
clus_spawn$obs <- factor(formatC(1:nrow(clus_spawn), flag = "0", width = 3)) # unique tank ID for later on

library(glmmTMB)
md1 <- glm(dummy ~ cluster, family = "binomial", data = clus_spawn)
summary(md1)

