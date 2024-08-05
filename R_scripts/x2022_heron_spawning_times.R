#2022 Heron spawning time pairwise differences i.eis spawning times related to spatial distances

# 1 Import data -----------------------------------------------------------
#read.excel <- function(header=TRUE,...) {read.table("clipboard", sep= "\t", dec= ",", header=header, na.strings=c("",".","na"),...)}
#data1=read.excel() #read clipboard from excel
#save(data1, file = file.path("./Rdata", "2022 spawning times.RData"))
load('./Rdata/2022 spawning times.RData') #data1
load("./Rdata/2022_Heron_dist_mat.RData") #lst

# 2 Labeling and wrangling -----------------------------------------------
str(data1) #check data type is correct
# data1$raw.x <- as.numeric(as.character(data1$raw.x))
# data1$time <- as.factor(as.character(data1$time))
# data1$tank <- as.factor(as.character(data1$tank))
data1 = data1[complete.cases(data1), ]  #make sure import matches NA type
data1

# 3 Data exploration ------------------------------------------------------
##Visualize data - plot data split at every factor
library(ggplot2)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2")  #set theme in code
# p0 = ggplot()+
#   geom_point(data1, mapping = aes(x = raw.x, y = suc/tot),position = position_jitter(width = .02, height = .02), 
#              alpha = 0.50,size = data1$tot/max(data1$tot)*3 )
# p0 = p0 + facet_wrap(~time+spec)+scale_x_log10(name ="XXXX")#+geom_smooth(data1, mapping = aes(x = raw.x, y = suc/tot))
# #p0= p0+ scale_y_continuous( limits = c(0, 1)) 
# p0

#data1 %>%  group_by(date) %>%  hist(.$tfll)

# pairwise differences ----------------------------------------------------
data2 = split(data1, data1$date)
#as.numeric(dist(data2$`13_12`$tfll))
#data2$`13_12`

data3.1 = data2$`13_12`
time_mat <- dist(data3.1[,3], diag = T)
time_mat

data4.1 = data2$`14_12`
time_mat1 <- dist(data4.1[,3], diag = T)
time_mat1

dist_mat = lst$dist_mat
dist_mat1 = lst$dist_mat1

# Perform the Mantel test for the 13th
library(vegan)
(mantel_result <- mantel(time_mat, dist_mat))
 #strong relations ship (r = ~1) but no sig. effect of distance on timing

# Perform the Mantel test for the 14th
(mantel_result1 <- mantel(time_mat1, dist_mat1))
 #almost no relationship and no effect of distance on timing


# library(dplyr)
# data1 %>% 
#   group_split(date) %>% 
#   lapply(function(x) {
#     combinations <- expand.grid(ID1 = unique(x$ID), ID2 = unique(x$ID))   #crates all combination
#     left_join(combinations, x, by = c("ID1" = "ID")) %>%
#       rename(time1 = tfll) %>%
#       left_join(x, by = c("ID2" = "ID")) %>%
#       rename(time2 = tfll) %>%
#       mutate(time_diff = time2 - time1) %>% 
#       filter(ID1 <= ID2) %>%
#       pivot_wider(id_cols = ID1, names_from = ID2, values_from = time_diff, names_prefix = "") %>% #creates a pairwise table
#       data.frame()
#   })


# cluster vs spawning time ------------------------------------------------
#cat spawning times
meta2 = data_gl_filtered_adult@other$ind.metrics
meta2 <- na.omit(meta2, subset = "spawn_time") # Remove rows where spawn_time is NA
meta2$spawn_time_clean = gsub("[^0-9]", "", meta2$spawn_time) %>% as.numeric()# Replace non-numerals with empty string
str(meta2)
meta2$spawn_time_cat = ifelse(meta2$spawn_time_clean <45, 'early', 'late')
meta2$dummy = ifelse(meta2$spawn_time_cat == 'early', 0, 1)
clus_spawn = meta2 %>%  dplyr::select(id, cluster, spawn_time_clean, spawn_time_cat, rep_id, dummy)
#nrow(clus_spawn)/2
#subset(clus_spawn, rep_id == 1)
clus_spawn = subset(clus_spawn, rep_id == 2)  #used rep '2' because it had an extra
clus_spawn = subset(clus_spawn, id != 'pd12.a.2')  #remove 12. 
clus_spawn$cluster <- droplevels(clus_spawn$cluster)   #drop unused levels




#analysis
(table <- table(clus_spawn$cluster, clus_spawn$spawn_time_cat))  # Create a contingency table
clus_spawn[which(clus_spawn$cluster == '3' & clus_spawn$spawn_time_cat == 'late'),]
chisq.test(table)  # Perform the chi-squared test
(fisher.test(table))  # Perform Fisher's Exact test
ggplot(clus_spawn, aes(x = spawn_time_clean, fill = factor(cluster))) +  # Set spawn_time_clean as x, and fill colour by cluster
  geom_density(alpha = 0.5) +  # Adjust transparency with alpha if needed
  labs(x = "Spawn Time Clean", y = "Density", fill = "Cluster") +  # Label the axes and legend
  ggtitle("Density Plot of Spawn Time by Cluster") +  # Add a title
  theme_minimal()  # Use a minimal theme


#try binary glm (sample size probably too small)
str(clus_spawn)
clus_spawn$obs <- factor(formatC(1:nrow(clus_spawn), flag = "0", width = 3)) # unique tank ID for later on

library(glmmTMB)
md1 <- glm(dummy ~ cluster, family = "binomial", data = clus_spawn)
summary(md1)

