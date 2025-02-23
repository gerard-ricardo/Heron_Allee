
library(tidyverse)
library(ggplot2)
library(tidyr)
library(rgdal)
library(dplyr)
library(spatstat)
library(tidybayes)
library(viridis)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2") # set theme in code


load("./Rdata/2022_Heron.RData")



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

range(data1$x1)
cord.dec <- SpatialPoints(cbind(data1$x1, -data1$y1), proj4string = CRS("+proj=longlat")) # convert to spatial points
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:32748")) %>%
  data.frame() %>%
  mutate(., x = coords.x1, y = coords.x2) # convert to northing &   easting
data1$x <- cord.UTM$x
data1$y <- cord.UTM$y


p0 <- ggplot() +
  geom_point(data1, mapping = aes(x = x, y = y), position = position_jitter(width = .02, height = .02), alpha = 0.50, size = 3)
p0 + geom_text(data1, mapping = aes(x = x, y = y, label = desc))


rotatefun <- function(angle) {
  boxpts1 <- SpatialPoints(t(bbox(cord.dec)), proj4string = CRS(proj4string(cord.dec)))
  boxLL2 <- boxpts1 %>%
    data.frame() %>%
    t() # wrangle to useable format
  llc1 <- apply(boxLL2, 1, mean)
  names(llc1) <- c("x", "y")
  prj1 <- paste0(
    "+proj=omerc +lat_0=", llc1[2], " +lonc=", llc1[1], " +alpha=",
    angle, " +gamma=0.0 +k=1.000000 +x_0=0.000 +y_0=0.000 +ellps=WGS84 +units=m "
  )
  ddd <- CRS(prj1)
  regionR1 <- spTransform(cord.dec, ddd)
  data2 <- regionR1 %>%
    data.frame() %>%
    mutate(., x = coords.x1, y = coords.x2)
  return(data2)
}

data3 <- rotatefun(angle = 50)
data3$id <- data1$desc
data3$cluster = data1$cluster


plot(data3$x, data3$y, xlim = c(60, -60), ylim = c(60, -60))
text(data3$x, data3$y, data3$id, pos = 1)

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
  tidybayes::stat_pointinterval(aes(y = 0.00, x = nearne), .width = c(.66, .95)) +
  coord_cartesian(ylim = c(0.0, 0.115)) +
  scale_x_continuous(name = "Intercolonial distance (m)") +
  scale_y_continuous(name = "Frequency") +
  theme_sleek3()+
  scale_fill_viridis_c(option = "viridis", direction = -1, name = "Intercolonial Distance")
p1



load("./Rdata/heron_intercol_all.RData")  #p1




data3 <- rotatefun(angle = 50)
data3$id <- data1$desc
data3$cluster = data1$cluster


plot(data3$x, data3$y, xlim = c(60, -60), ylim = c(60, -60))
text(data3$x, data3$y, data3$id, pos = 1)


G <- Gest(mypattern)
plot(G, cbind(km, rs, han) ~ r, main = "Nearest neighbor distance distribution")

pd1 <- pairdist(mypattern)
colnames(pd1) <- data3$id
rownames(pd1) <- data3$id

data5 <- subset(data3, id %in% c("4", "13", "14", "15")) # remove factor treatment level. Use '%in%' to keep.
dist_mat <- dist(data5[, 3:4], diag = T)
dist_mat

data6 <- subset(data3, id %in% c("3", "5", "9", "13", "14", "15")) # remove factor treatment level. Use '%in%' to keep.
dist_mat1 <- dist(data6[, 3:4], diag = T)
dist_mat1

lst = list(dist_mat = dist_mat, dist_mat1 = dist_mat1)

clarkevans(mypattern) # clarke-evans (A value R>1 suggests ordering, while R<1 suggests clustering.)
plot(density(mypattern), main = "")

density_df <- as.data.frame(as.table(as.matrix(density(mypattern))))

density_df <- density_df[, c("Var2", "Var1", "Freq")]
names(density_df) <- c("x", "y", "value")

convert_to_numeric <- function(column) {
  as.numeric(factor(column, levels = unique(column)))
}

density_df$x <- convert_to_numeric(density_df$x)
density_df$y <- convert_to_numeric(density_df$y)
head(density_df)
tail(density_df)

density_df$x <- max(density_df$x) - density_df$x + 1
density_df$y <- max(density_df$y) - density_df$y + 1

aspect_ratio <- 3.292593





theme_minimal_no_grid <- function(base_size = 11, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.ticks.length = unit(0.25, "cm"),  # Set the length of axis ticks
      axis.ticks = element_line(color = "grey20"),  # Keep axis ticks visible
      panel.background = element_rect(fill = "white", color = NA),  # White background for the plot panel
      plot.background = element_rect(fill = "white", color = NA)    # White background for the entire plot
    )
}



p2 <- ggplot(density_df, aes(x, y, fill=value)) +
  geom_raster() +
  scale_fill_viridis_c(option = "viridis", name = "Density") +  # Add "Density" as the legend title
  coord_fixed(ratio = aspect_ratio) +  # Adjust the aspect ratio for the rotated plot
  theme_minimal_no_grid() +
  scale_x_continuous(
    breaks = c(0, 64, 128),  # Original x values
    labels = c(0, 20, 40)    # New x labels
  ) +
  labs(x = "Cross-shore",  y = "Alongshore",  fill = "Density")
p2

ggsave(p2, filename = 'adult_density.pdf',  path = "./plots", device = 'pdf',  width = 6, height = 5.5)  #

load("./Rdata/heron_density.RData")



df <- data.frame(id = c("5", "9", "15", "13", "14"), dist = c(158, 165, 148, 176, 159))
mean(df$dist)
sd(df$dist)



data3 %>% filter(id %in% c('1', '15', '12'))   %>%
  mutate(statement = paste("ID", id, "corresponds to cluster", cluster)) %>%
  pull(statement) %>% walk(print)  # Using walk from purrr to print each statement

data4 = data3 %>% subset(., cluster %in% '3') 
rangex <- range(data4$x) + cbind(-2, 2) 
rangey <- range(data4$y) + cbind(-2, 2)   #plus buffer?
mypattern <- ppp(data4$x, data4$y, c(rangex), c(rangey), marks = data4$id) # imports using subset as x and y, then ranges of x and y
plot(mypattern)
sum_box = summary(mypattern)  
sum_box$intensity  
nearne <- nndist(mypattern) # Computes the distance from each point to its nearest neighbour
quantile(nndist(mypattern))
p1 <- ggplot() +
  geom_density(aes(nearne), alpha = 0.3, color = "steelblue", fill = "steelblue") +
  tidybayes::stat_pointinterval(aes(y = 0.00, x = nearne), .width = c(.66, .95)) #+facet_wrap(~contrast+time, nrow = 3, ncol = 2)+
p1 <- p1 + coord_cartesian(ylim = c(0.0, 0.15))
p1 <- p1 + scale_x_continuous(name = "Nearest neighbour distance (m)")
p1 <- p1 + scale_y_continuous(name = "Frequency")
p1

data5 = data3 %>% subset(., cluster%in% '2') 
rangex <- range(data5$x) + cbind(-2, 2) 
rangey <- range(data5$y) + cbind(-2, 2)   #plus buffer?
mypattern <- ppp(data5$x, data5$y, c(rangex), c(rangey), marks = data5$id) # imports using subset as x and y, then ranges of x and y
plot(mypattern)
sum_box = summary(mypattern)  
sum_box$intensity  #0.003377019 points per square unit

nearne <- nndist(mypattern) # Computes the distance from each point to its nearest neighbour
quantile(nndist(mypattern), c(0.025, 0.5, 0.975))

p1 <- ggplot() +
  geom_density(aes(nearne), alpha = 0.3, color = "steelblue", fill = "steelblue") +
  tidybayes::stat_pointinterval(aes(y = 0.00, x = nearne), .width = c(.66, .95)) #+facet_wrap(~contrast+time, nrow = 3, ncol = 2)+
p1 <- p1 + coord_cartesian(ylim = c(0.0, 0.15))
p1 <- p1 + scale_x_continuous(name = "Nearest neighbour distance (m)")
p1 <- p1 + scale_y_continuous(name = "Frequency")
p1
