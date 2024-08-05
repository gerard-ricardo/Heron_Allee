######### (Heron adult colonies coords)#########################

# 1. Load Libraries ------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(tidyr)
library(rgdal)
library(dplyr)
library(spatstat)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2") # set theme in code

# 1 Import data -----------------------------------------------------------

# read.excel <- function(header=TRUE,...) {read.table("clipboard",sep="\t",header=header,...)}
# data1=read.excel() #read clipboard from excel
#
# save(data1, file = file.path("./Rdata", "2022_Heron.RData"))
load("./Rdata/2022_Heron.RData")

# note: this needs to rotate and a correct envelope used.

# 2 Labeling and wrangling -----------------------------------------------
str(data1) # check data type is correct
data1$y1 <- as.numeric(as.character(data1$latitude))
data1$x1 <- as.numeric(as.character(data1$longitude))
data1
data1 <- data1[complete.cases(data1), ] # make sure import matches NA type
data1 = data1 %>% arrange(desc)

pca_complete_clus = pca_complete %>% dplyr::select(c(pop, cluster))
pca_complete_clus_unique <- pca_complete_clus %>%
  group_by(pop) %>%
  slice(1) %>%
  ungroup()
data1$pop = data1$desc
data1 = left_join(data1, pca_complete_clus_unique)

# 3 Data exploration ------------------------------------------------------
#### spatial conversion
range(data1$x1)
cord.dec <- SpatialPoints(cbind(data1$x1, -data1$y1), proj4string = CRS("+proj=longlat")) # convert to spatial points
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:32748")) %>%
  data.frame() %>%
  mutate(., x = coords.x1, y = coords.x2) # convert to northing &   easting
data1$x <- cord.UTM$x
data1$y <- cord.UTM$y

## Visualize data - plot data split at every factor

p0 <- ggplot() +
  geom_point(data1, mapping = aes(x = x, y = y), position = position_jitter(width = .02, height = .02), alpha = 0.50, size = 3)
p0 + geom_text(data1, mapping = aes(x = x, y = y, label = desc))


# rotating to north -------------------------------------------------------

#http://rstudio-pubs-static.s3.amazonaws.com/9694_61d30d86c193466ab38c7bea58221e35.html
rotatefun <- function(angle) {
  boxpts1 <- SpatialPoints(t(bbox(cord.dec)), proj4string = CRS(proj4string(cord.dec)))
  # convert to lat-long
  # boxLL1 = bbox(spTransform(boxpts1, CRS("+init=epsg:32748")))
  boxLL2 <- boxpts1 %>%
    data.frame() %>%
    t() # wrangle to useable format
  # find the centre
  llc1 <- apply(boxLL2, 1, mean)
  names(llc1) <- c("x", "y")
  # construct the proj4 string
  prj1 <- paste0(
    "+proj=omerc +lat_0=", llc1[2], " +lonc=", llc1[1], " +alpha=",
    angle, " +gamma=0.0 +k=1.000000 +x_0=0.000 +y_0=0.000 +ellps=WGS84 +units=m "
  )
  # return as a CRS:
  ddd <- CRS(prj1)
  regionR1 <- spTransform(cord.dec, ddd)
  data2 <- regionR1 %>%
    data.frame() %>%
    mutate(., x = coords.x1, y = coords.x2)
  # data3 = data2 %>% t() #%>% mutate(.,x = coords.x1, y = coords.x2)
  return(data2)
}

data3 <- rotatefun(angle = 50)
data3$id <- data1$desc
data3$cluster = data1$cluster

# p0 = ggplot()+geom_point(out, mapping = aes(x = x, y = y),position = position_jitter(width = .02, height = .02), alpha = 0.50,size = 3 )
# p0

plot(data3$x, data3$y, xlim = c(60, -60), ylim = c(60, -60))
text(data3$x, data3$y, data3$id, pos = 1)
# note data are transposed

############################
rangex <- range(data3$x) + cbind(-2, 2) 
rangey <- range(data3$y) + cbind(-2, 2)   #plus buffer?
mypattern <- ppp(data3$x, data3$y, c(rangex), c(rangey), marks = data3$id) # imports using subset as x and y, then ranges of x and y
plot(mypattern)
sum_box = summary(mypattern)  
sum_box$intensity  #0.003591173 points per square unit

1 / 0.004113922  #243 corals er
nearne <- nndist(mypattern) # Computes the distance from each point to its nearest neighbour
quantile(nndist(mypattern))

p1 <- ggplot() +
  geom_density(aes(nearne), alpha = 0.3, color = "steelblue", fill = "steelblue") +
  tidybayes::stat_pointinterval(aes(y = 0.00, x = nearne), .width = c(.66, .95)) #+facet_wrap(~contrast+time, nrow = 3, ncol = 2)+
p1 <- p1 + coord_cartesian(ylim = c(0.0, 0.15))
p1 <- p1 + scale_x_continuous(name = "Nearest neighbour distance (m)")
p1 <- p1 + scale_y_continuous(name = "Frequency")
p1

#ggsave(p1, filename = 'heron_nn_dist.tiff',  path = "./plots", device = "tiff",  width = 8, height = 5)  #this often works better than pdf

G <- Gest(mypattern)
plot(G, cbind(km, rs, han) ~ r, main = "Nearest neighbor distance distribution")

# pairwise distance
pd1 <- pairdist(mypattern)
colnames(pd1) <- data3$id
rownames(pd1) <- data3$id

# get distance matrix for colonies spawning on 13th
data5 <- subset(data3, id %in% c("4", "13", "14", "15")) # remove factor treatment level. Use '%in%' to keep.
dist_mat <- dist(data5[, 3:4], diag = T)
dist_mat

# get distance matrix for 14th
data6 <- subset(data3, id %in% c("3", "5", "9", "13", "14", "15")) # remove factor treatment level. Use '%in%' to keep.
dist_mat1 <- dist(data6[, 3:4], diag = T)
dist_mat1

lst = list(dist_mat = dist_mat, dist_mat1 = dist_mat1)
#save(lst, file = file.path("./Rdata", "2022_Heron_dist_mat.RData"))

clarkevans(mypattern) # clarke-evans (A value R>1 suggests ordering, while R<1 suggests clustering.)
plot(Kest(mypattern, correction = c("best"))) # Ripley’s K-function. Above = clustered, below = regular.
plot(envelope(mypattern, fun = Kest, nsim = 780, nrank = 20)) # Envelopes of K-function: This is a hypothesis test
plot(envelope(mypattern, Kest, correction = "Ripley", verbose = F))
plot(envelope(mypattern, Lest, correction = "Ripley", verbose = F)) # this is often preferred over K test (transformed K)
plot(density(mypattern), main = "Density") # kernel smoother of point density:
p2

#### distance parallel to shore##################
df <- data.frame(id = c("5", "9", "15", "13", "14"), dist = c(158, 165, 148, 176, 159))
mean(df$dist)
sd(df$dist)
#########################


# split by clusters -------------------------------------------------------
## group 1
data4 = data3 %>% subset(., cluster%in% '1') 
rangex <- range(data4$x) + cbind(-2, 2) 
rangey <- range(data4$y) + cbind(-2, 2)   #plus buffer?
mypattern <- ppp(data4$x, data4$y, c(rangex), c(rangey), marks = data4$id) # imports using subset as x and y, then ranges of x and y
plot(mypattern)
sum_box = summary(mypattern)  
sum_box$intensity  #0.002426236 points per square unit

#1 / sum_box$intensity  #412 corals er
nearne <- nndist(mypattern) # Computes the distance from each point to its nearest neighbour
quantile(nndist(mypattern))
#0%       25%       50%       75%      100% 
#5.258298 12.430135 14.437074 16.053784 31.39961
p1 <- ggplot() +
  geom_density(aes(nearne), alpha = 0.3, color = "steelblue", fill = "steelblue") +
  tidybayes::stat_pointinterval(aes(y = 0.00, x = nearne), .width = c(.66, .95)) #+facet_wrap(~contrast+time, nrow = 3, ncol = 2)+
p1 <- p1 + coord_cartesian(ylim = c(0.0, 0.15))
p1 <- p1 + scale_x_continuous(name = "Nearest neighbour distance (m)")
p1 <- p1 + scale_y_continuous(name = "Frequency")
p1

## group 3
data5 = data3 %>% subset(., cluster%in% '3') 
rangex <- range(data5$x) + cbind(-2, 2) 
rangey <- range(data5$y) + cbind(-2, 2)   #plus buffer?
mypattern <- ppp(data5$x, data5$y, c(rangex), c(rangey), marks = data5$id) # imports using subset as x and y, then ranges of x and y
plot(mypattern)
sum_box = summary(mypattern)  
sum_box$intensity  #0.003377019 points per square unit

#1 / sum_box$intensity  #412 corals er
nearne <- nndist(mypattern) # Computes the distance from each point to its nearest neighbour
quantile(nndist(mypattern), c(0.025, 0.5, 0.975))
#0%       25%       50%       75%      100% 
#5.258298 12.430135 14.437074 16.053784 31.39961
p1 <- ggplot() +
  geom_density(aes(nearne), alpha = 0.3, color = "steelblue", fill = "steelblue") +
  tidybayes::stat_pointinterval(aes(y = 0.00, x = nearne), .width = c(.66, .95)) #+facet_wrap(~contrast+time, nrow = 3, ncol = 2)+
p1 <- p1 + coord_cartesian(ylim = c(0.0, 0.15))
p1 <- p1 + scale_x_continuous(name = "Nearest neighbour distance (m)")
p1 <- p1 + scale_y_continuous(name = "Frequency")
p1
