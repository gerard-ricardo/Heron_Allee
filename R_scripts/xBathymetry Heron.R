#Bathymetry data Heron


# 1. Load Libraries ------------------------------------------------------
library(raster) # Load the raster package
library(dplyr)
library(rgdal) # Load the rgdal package
library(sp) # Load the sp package
library(RColorBrewer) # Load the RColorBrewer package
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2") # set theme in code


# # 1 Import data -----------------------------------------------------------
# #file_path <- "C:/Users/gerar/OneDrive/1_Work/3_Results/0 Data for experiental design/Capricorn Bunker/heron_1m_1.asc" # Specify the file path
# #file_path <- "C:/Users/gerar/OneDrive/1_Work/3_Results/0 Data for experiental design/Capricorn Bunker/Heron_bath.tif" # Specify the file path
# file_path <- "C:/Users/gerar/OneDrive/1_Work/3_Results/0 Data for experiental design/Capricorn Bunker/BATHY_1M_depth_detail/HRCASI_full_zpl_pair13_dist_BATHY_1M_depth_detail.tif" # Specify the file path
# #tif is much faster
# #this is the extent of the original file 'HRCASI_full_zpl_pair13_dist_BATHY_1M_depth_detail.tif'
# raster_data <- raster(file_path) # Read the .asc file into a raster object
# 
# 
# # 3. Data wrangling and exploration ----------------------------------------------------
# plot(raster_data) # Plot the raster data
# raster_data[raster_data == -999] <- NA # Replace -999 values with NA
# plot(raster_data) # Plot the raster data
# 
# raster_data@extent  
# #data in northing easting  form
# #Latitude Start: -23.432317 * Longitude Start: 151.898827 * Latitude End: -23.449153 * Longitude End: 151.924951
# 
# # Define the latitude and longitude coordinates
# #-23.453873° to  -23.454934°, 151.922649° to  151.923914°
#
# target_crs <- CRS("+proj=longlat +datum=WGS84")
# raster_data_latlon <- projectRaster(raster_data, crs = target_crs)
# raster_data = raster_data_latlon
# plot(raster_data, main = "Raster Data in Geographic Coordinates")
#save(raster_data, file = file.path("./Rdata", "Heron_bath_1m.RData"))
#load("./Rdata/Heron_bath_1m.RData")  #raster_data     - MOVED TO SCRATCH = TOO BIG
#
# Create a SpatialPoints object for the lat/lon coordinates
#lat <- c(-23.453873, -23.454934)
#lon <- c(151.9226, 151.923914)
#coordinates <- data.frame(lon, lat)
# coordinates_sp <- SpatialPoints(coordinates, proj4string = CRS("+proj=longlat +datum=WGS84"))
# # Define the projection (CRS) for UTM
# crs_utm <- CRS("+proj=utm +zone=56 +south +datum=WGS84 +units=m +no_defs")
# # Transform the coordinates to UTM
# coordinates_utm <- spTransform(coordinates_sp, crs_utm)
# # Extract the UTM coordinates
# (utm_coords <- coordinates(coordinates_utm)) 
# utm_coords = data.frame(utm_coords)
#
#extract AOI
#extent_to_crop <- extent(coordinates$lon[1], coordinates$lon[2], coordinates$lat[2], coordinates$lat[1]) # Define the extent
#cropped_raster <- crop(raster_data, extent_to_crop) # Crop the raster to the specified extent
#save(cropped_raster, file = file.path("./Rdata", "Heron_bath_1m_site.RData"))
load("./Rdata/Heron_bath_1m_site.RData")  #cropped_raster


# ## plot in base r
plot(cropped_raster) # Plot the raster data

#check values
values_tif <- getValues(cropped_raster)  #itff = 6 218 852 , 59 050 608
hist(values_tif, na.rm = T)

# Apply focal filter for smoothing
# Define the size of the moving window (3x3 window in this example)
window_size <- 3
# Apply the focal filter with mean function
smoothed_raster <- focal(cropped_raster, w = matrix(1, window_size, window_size), fun = mean, na.rm = TRUE)
# plot(smoothed_raster) # Plot the smoothed raster data
# #check values
# values_tif2 <- getValues(smoothed_raster)  #itff = 6 218 852 , 59 050 608 
# hist(values_tif2, na.rm = T)
# #custom_palette <- colorRampPalette(rev(c("royalblue3", "royalblue2",  "royalblue1", 'steelblue3', 'steelblue2', "steelblue1", 'paleturquoise2', "paleturquoise1", "lightcyan")))
# 
# # Assuming smoothed_raster is already created and processed
# 
# # Define a ColorBrewer palette (GnBu) and create a gradient palette
# gnbu_palette <- brewer.pal(5, "GnBu")
# custom_palette <- colorRampPalette(gnbu_palette)
# 
# # Plot the smoothed raster data without the default legend
# plot(smoothed_raster, col = custom_palette(100), legend = FALSE)
# # Define the range of values for the legend
# legend_values <- seq(minValue(smoothed_raster), maxValue(smoothed_raster), length.out = 100)
# # Define the position and size of the custom legend
# vz = matrix(legend_values, nrow = 1)
# vx = c(185, 195)  # x position of the legend
# vy = seq(-7.8, 10, length.out = 100)  # y position of the legend

## geom_raster
# Convert the smoothed raster to a data frame for ggplot
raster_df <- as.data.frame(rasterToPoints(smoothed_raster), stringsAsFactors = FALSE)
colnames(raster_df) <- c("lon", "lat", "value")  # Rename columns to meaningful names
# Plot the raster data using ggplot2
blues_palette <- brewer.pal(9, "YlGnBu")[3:9]  # Get a Blues palette with 9 colors
custom_palette <- colorRampPalette(blues_palette)
ggplot(raster_df, aes(x = lon, y = lat, fill = value)) +
  geom_raster() +  # or use geom_tile() as an alternative
  scale_fill_gradientn(colours = custom_palette(5), name = "Depth (m)") +
  theme_sleek2() +
  labs(x = "Longitude", y = "Latitude") +
  coord_fixed()  # Ensures that the aspect ratio is maintained





# add adults --------------------------------------------------------------
load("./Rdata/2022_Heron.RData") #data1

# note: this needs to rotate and a correct envelope used.

# 2 Labelling and wrangling -----------------------------------------------
# str(data1) # check data type is correct
# data1$y1 <- as.numeric(as.character(data1$latitude))
# data1$x1 <- as.numeric(as.character(data1$longitude))
# data1
# data1 <- data1[complete.cases(data1), ] # make sure import matches NA type
# data2 = data.frame(lat = data1$y1, lon = data1$x1, id = data1$desc)
# data3 = left_join(data2, data_gl_filtered_adult@other$ind.metrics, by = 'id')

meta1 = data_gl_filtered_adult@other$ind.metrics  #note PCA with cluster needs to run first
cluster_colours <- c("dodgerblue", "salmon", "mediumseagreen")
point_colours <- cluster_colours[meta1$cluster]
# points(meta1$lon, meta1$lat, pch = 21, col = "black", bg = point_colours, cex = 1.5)
# 
# 
# # Add the custom legend
# par(new = TRUE)
# plot(1, xlab = "Latitude", ylab = "Longitude", axes = FALSE, type = "n", xlim = c(-200, 200), ylim = c(-20, 20))
# image(vx, vy, vz, col = rev(custom_palette(100)), axes = FALSE, xlab = "", ylab = "", add = TRUE)
# polygon(c(185, 195, 195, 185), c(-7.8, -7.8, 10, 10))
# minValue(smoothed_raster)
# maxValue(smoothed_raster)
# axis(4, at = c(5.4, 0.4, -4.6), labels = c(5, 10, 15), las = 1)
# text(x = 185, y = 12, labels = 'Depth (m)')

# library(gridGraphics)
# library(grid)
# grid.echo()
# a <- grid.grab()
# grid.newpage()
# p5 = grid.draw(a)
# # Step 1: Create and capture the base plot
# recorded_plot <- recordPlot()
# 
# # Step 2: Replay the recorded plot and capture it as a grob
# grid.newpage()  # Clear the page
# replayPlot(recorded_plot)
# plot_grob <- grid.grab()

library(scales)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek1") # set theme in code

p5 <- ggplot(raster_df, aes(x = lon, y = lat, fill = value)) +
  geom_raster() +
  theme_sleek1() +
  scale_fill_gradientn(colours = custom_palette(100), name = "Depth (m)") +
  geom_point(data = meta1, mapping = aes(x = lon, y = lat), 
             shape = 21, size = 3, fill = cluster_colours[meta1$cluster], colour = 'grey10', alpha = 0.5) +
  scale_x_continuous(labels = label_number(accuracy = 0.0001)) +  # Set longitude to 4 decimal places
  scale_y_continuous(labels = label_number(accuracy = 0.0001)) +  # Set latitude to 4 decimal places
  labs(x = "Longitude", y = "Latitude") +
  coord_fixed() +
  theme(legend.position = c(0.85, 0.7))

p5

save(p5, file = file.path("./Rdata", "bath_cluster.RData"))
load("./Rdata/bath_cluster.RData")  #p5


# google maps -------------------------------------------------------------
data1 <- data1 %>% dplyr::rename(lat = latitude, lon = longitude) %>% dplyr::select(c(lat, lon, desc))

library(ggsn)
library(ggmap)
# p3 <- ggmap(get_googlemap(center = c(151.9234, -23.4544), zoom = 19,  maptype = "satellite")) +  #centre coordinat
#   geom_point(data = data1, aes(lon, lat, col = desc), alpha = 0.5, size = 2)+
#   labs(
#     x = "Longitude",  # Label for the x-axis
#     y = "Latitude",  # Label for the y-axis
#   )
# p3

#scale bar
# Get the map object
map <- get_googlemap(center = c(151.9233, -23.4544), zoom = 19, maptype = "satellite")
bbox <- attr(map, "bb") %>% data.frame()
bbox$ll.lon


library(ggsn)
library(ggmap)

# Create the map with points and a scale bar
p3 <- ggmap(map) +
  geom_point(data = data1, aes(x = lon, y = lat, col = desc), alpha = 0.5, size = 2) +
  geom_text(data = data1, aes(x = lon, y = lat, label = desc), vjust = -1, hjust = 0.5, size = 3, color = "white") + 
  labs(x = "Longitude", y = "Latitude") +
  # geom_rect(aes(xmin = bbox$ll.lon, xmax = bbox$ur.lon,  ymin = bbox$ll.lat , ymax = bbox$ll.lat- 0.05),fill = "black") +
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

#save(p3, file = file.path("./Rdata", "heron_adult_site.RData"))
load("./Rdata/heron_adult_site.RData") #p3



