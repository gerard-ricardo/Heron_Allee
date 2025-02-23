

library(raster) # Load the raster package
library(dplyr)
library(rgdal) # Load the rgdal package
library(sp) # Load the sp package
library(RColorBrewer) # Load the RColorBrewer package
library(ggsn)
library(ggmap)
library(scales)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(units)
library(sf)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2") # set theme in code
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek1") # set theme in code

load("./Rdata/2022_Heron.RData") #data1
data1 <- data1 %>% dplyr::rename(lat = latitude, lon = longitude) %>% dplyr::select(c(lat, lon, desc))


map <- get_googlemap(center = c(151.9233, -23.4544), zoom = 19, maptype = "satellite")
bbox <- attr(map, "bb") %>% data.frame()
bbox$ll.lon

p3 <- ggmap(map) +
  geom_point(data = data1, aes(x = lon, y = lat, col = desc), alpha = 0.5, size = 2) +
  geom_text(data = data1, aes(x = lon, y = lat, label = desc), vjust = -1, hjust = 0.5, size = 3, color = "white") + 
  labs(x = "Longitude", y = "Latitude") +
  ggsn::scalebar(
    x.min = bbox$ll.lon,  # Lower-left longitude
    x.max = bbox$ur.lon,  # Upper-right longitude
    y.min = bbox$ll.lat,  # Lower-left latitude
    y.max = bbox$ur.lat,  # Upper-right latitude
    dist = 20,  # Distance represented by the scale bar in metres (now 20m)
    dist_unit = "m",  # Units of the scale bar
    transform = TRUE,  # Convert to planar coordinates for accuracy
    model = "WGS84",  # Coordinate reference system
    location = "bottomright",  # Location of the scale bar
    st.dist = 0.5,  # Distance between the scale bar and text, increased to move it higher
    st.size = 3,  # Text size
    st.color = "white",  # Text color
    height = 0.03,  # Scale bar height
    box.fill = c("white", "white"),  # Fill colors
    box.color = "black"  # Border color
  )+
  scale_y_continuous(limits = c(-23.45514, bbox$ur.lat), expand = c(0, 0))
p3

load("./Rdata/heron_adult_site.RData") #p3



load("./Rdata/Heron_bath_1m_site.RData")  #cropped_raster
load("./Rdata/2022_Heron.RData") #data1

plot(cropped_raster) # Plot the raster data

values_tif <- getValues(cropped_raster)  #itff = 6 218 852 , 59 050 608
hist(values_tif, na.rm = T)

window_size <- 3
smoothed_raster <- focal(cropped_raster, w = matrix(1, window_size, window_size), fun = mean, na.rm = TRUE)

raster_df <- as.data.frame(rasterToPoints(smoothed_raster), stringsAsFactors = FALSE)
colnames(raster_df) <- c("lon", "lat", "value")  # Rename columns to meaningful names
blues_palette <- brewer.pal(9, "YlGnBu")[3:9]  # Get a Blues palette with 9 colors
custom_palette <- colorRampPalette(blues_palette)
ggplot(raster_df, aes(x = lon, y = lat, fill = value)) +
  geom_raster() +  # or use geom_tile() as an alternative
  scale_fill_gradientn(colours = custom_palette(5), name = "Depth (m)") +
  theme_sleek2() +
  labs(x = "Longitude", y = "Latitude") +
  coord_fixed()  # Ensures that the aspect ratio is maintained








bathy_plot <- ggplot(raster_df, aes(x = lon, y = lat, fill = value)) +
  geom_raster() +
  theme_sleek2() +
  scale_fill_gradientn(colours = custom_palette(100), name = "Depth (m)") +
  geom_point(data = data_genind_adult_unique@other$ind.metrics, mapping = aes(x = lon, y = lat), shape = 21, size = 3, 
             fill = data_genind_adult_unique@other$ind.metrics$cluster_colour, colour = 'grey10', alpha = 0.5) + 
  scale_x_continuous(labels = label_number(accuracy = 0.0001), breaks = scales::breaks_extended(n = 4)) +   
  scale_y_continuous(labels = label_number(accuracy = 0.0001)) +  
  labs(x = "Longitude", y = "Latitude") +
  coord_fixed(expand = FALSE) +  
  theme(
    legend.position = c(0.88, 0.75),
    legend.title = element_text(size = 10),  # Set legend title font size
    axis.text.x = element_text(size = 9, margin = margin(t = -15)),  # Move x-axis tick labels inside
    axis.text.y = element_text(size = 9, margin = margin(r = -40)),  # Move y-axis tick labels inside
    axis.title.x = element_text(margin = margin(t = 10)),  # Move x-axis label further down
    axis.title.y = element_text(margin = margin(r = 10)))   # Move y-axis label further left

bathy_plot


points_df <- data_genind_adult_unique@other$ind.metrics %>% dplyr::select(lon, lat, new_id, cluster_colour)
find_nearest <- function(point, raster_df) {
  distances <- sqrt((raster_df$lon - point["lon"])^2 + (raster_df$lat - point["lat"])^2)  # Calculate distances
  nearest_index <- which.min(distances)  # Find the index of the nearest point
  return(c(nearest_index, distances[nearest_index]))
}
nearest_indices <- apply(points_df[, c("lon", "lat")], 1, find_nearest, raster_df = raster_df)
matched_values <- points_df %>%
  mutate(
    nearest_index = nearest_indices[1, ],  # Add the index of the nearest raster point
    distance = nearest_indices[2, ],  # Add the distance to the nearest raster point
    value = raster_df$value[nearest_indices[1, ]]  # Add the value of the nearest raster point
  )

filtered_df <- matched_values %>% filter(cluster_colour != "mediumseagreen") # Filter out 'mediumseagreen'
unique_colours <- unique(filtered_df$cluster_colour)
(t_test_result <- t.test(value ~ cluster_colour, data = filtered_df))










library(ggmap)
library(ggplot2)
library(sf)
library(dplyr)
library(patchwork)

data1 <- data.frame(
  ID = c('305_1', '10_1', '24_1', '11_1', '25_1', '306_1'),   #id_TIMEPOINT
  lat1 = c(-23.450642, -23.450656, -23.450422, -23.462564, -23.462383, -23.462297),
  lon1 = c(151.916614, 151.916667, 151.916719, 151.928894, 	151.928697, 151.928564),
  second = c('305_7', '10_2', '24_2', '11_8', '25_8', '306_8'), 
  lat2 = c(-23.453875, -23.454042, -23.453064, -23.449714, -23.447156, -23.449331),
  lon2 = c(151.922756, 151.922786, 151.921297, 151.912953, 151.9089, 151.912317)
)

load("./Rdata/heron_container_tracks.RData") #data2
data2
data2 <- data2[complete.cases(data2), ] # make sure import matches NA type
data2$time = data2$time_min * 60



data2 <- data2 %>%
  filter(!is.na(ID)) %>% # Remove rows with NA in ID
  group_by(date, ID) %>% # Group by ID
  filter(time_min == max(time_min, na.rm = TRUE)) %>% # Filter for max time_min in each group
  ungroup() %>% data.frame() # Ungroup for clean output

data2 = data2[!(data2$date == '13th' & data2$ID == '14'), ] #remove container with broken filter
data2 = data2[!(data2$date == '14th' & data2$ID == '5'), ] #remove containers that was mashed



points_x <- st_as_sf(data1, coords = c("lon1", "lat1"), crs = 4326)
points_y <- st_as_sf(data1, coords = c("lon2", "lat2"), crs = 4326)
points_x1 <- st_as_sf(data2, coords = c("lon1", "lat1"), crs = 4326)
points_y1 <- st_as_sf(data2, coords = c("lon2", "lat2"), crs = 4326)

points_x_proj <- st_transform(points_x, crs = 32756)
points_y_proj <- st_transform(points_y, crs = 32756)
points_x_proj1 <- st_transform(points_x1, crs = 32756)
points_y_proj1 <- st_transform(points_y1, crs = 32756)

data1$dist <- st_distance(points_x_proj, points_y_proj, by_element = TRUE)
data1$dist_m <- as.numeric(data1$dist)
data2$dist <- st_distance(points_x_proj1, points_y_proj1, by_element = TRUE)
data2$dist_m <- as.numeric(data2$dist)

data1$time <- c(3600, 3600, 3000, 4200, 4200, 4200)

data1$speed_m_s <- data1$dist_m / data1$time 
(data2$speed_m_s <- data2$dist_m / data2$time)


map2 <- get_googlemap(center = c(151.9225, -23.4545), zoom = 17, maptype = "satellite")

p2 <- ggmap(map2) +
  geom_segment(data = data2, mapping = aes(x = lon1, y = lat1, xend = lon2, yend = lat2, size = speed_m_s, color = date), 
               size = 1, alpha = 0.9, arrow = arrow(length = unit(0.3, "cm"), type = "closed")) +
  labs(x = "Longitude", y = "Latitude", size = "Speed (m/s)") +
  theme_minimal() +
  scale_size_continuous(range = c(0.5, 2))+
  theme(legend.position = c(0.85, 0.85), legend.background = element_rect(fill = "white", colour = "black")  # White background with a border
  ) 
p2

map <- get_googlemap(center = c(151.9195, -23.452), zoom = 15, maptype = "satellite")
p1 <- ggmap(map) +
  geom_segment(data = data1, mapping = aes(x = lon1, y = lat1, xend = lon2, yend = lat2, size = speed_m_s), 
               color = "red", alpha = 0.4, arrow = arrow(length = unit(0.3, "cm"), type = "closed")) +
  geom_point(mapping = aes(x = 151.923240, y = -23.454442), colour = "grey", size = 4) + # Add the dark blue point
  labs(x = "Longitude", y = "Latitude", size = "Speed (m/s)") +
  scale_size_continuous(range = c(0.5, 2)) +
  theme_minimal() +
  scale_size_continuous(range = c(0.5, 2)) +
  theme(legend.position = c(0.85, 0.80), legend.background = element_rect(fill = "white", colour = "black")  # White background with a border
  ) 
p1


panel_plot <- p2 | p1 
panel_plot



