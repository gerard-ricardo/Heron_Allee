#2022 Heron spawning time pairwise differences i.eis spawning times related to distances

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




