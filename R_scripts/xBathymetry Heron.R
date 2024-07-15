#Bathymetry data Heron


# 1. Load Libraries ------------------------------------------------------
library(raster) # Load the raster package
library(dplyr)
library(rgdal) # Load the rgdal package
library(sp) # Load the sp package
library(RColorBrewer) # Load the RColorBrewer package




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
plot(cropped_raster) # Plot the raster data

#check values
values_tif <- getValues(cropped_raster)  #itff = 6 218 852 , 59 050 608 
hist(values_tif, na.rm = T)

# Apply focal filter for smoothing
# Define the size of the moving window (3x3 window in this example)
window_size <- 3
# Apply the focal filter with mean function
smoothed_raster <- focal(cropped_raster, w = matrix(1, window_size, window_size), fun = mean, na.rm = TRUE)
plot(smoothed_raster) # Plot the smoothed raster data
#check values
values_tif2 <- getValues(smoothed_raster)  #itff = 6 218 852 , 59 050 608 
hist(values_tif2, na.rm = T)
#custom_palette <- colorRampPalette(rev(c("royalblue3", "royalblue2",  "royalblue1", 'steelblue3', 'steelblue2', "steelblue1", 'paleturquoise2', "paleturquoise1", "lightcyan")))


# Assuming smoothed_raster is already created and processed

# Define a ColorBrewer palette (GnBu) and create a gradient palette
gnbu_palette <- brewer.pal(5, "GnBu")
custom_palette <- colorRampPalette(gnbu_palette)

# Plot the smoothed raster data without the default legend
plot(smoothed_raster, col = custom_palette(100), legend = FALSE)
# Define the range of values for the legend
legend_values <- seq(minValue(smoothed_raster), maxValue(smoothed_raster), length.out = 100)
# Define the position and size of the custom legend
vz = matrix(legend_values, nrow = 1)
vx = c(185, 195)  # x position of the legend
vy = seq(-7.8, 10, length.out = 100)  # y position of the legend






# add adults --------------------------------------------------------------
load("./Rdata/2022_Heron.RData")

# note: this needs to rotate and a correct envelope used.

# 2 Labelling and wrangling -----------------------------------------------
str(data1) # check data type is correct
data1$y1 <- as.numeric(as.character(data1$latitude))
data1$x1 <- as.numeric(as.character(data1$longitude))
data1
data1 <- data1[complete.cases(data1), ] # make sure import matches NA type
data2 = data.frame(lat = data1$y1, lon = data1$x1, id = data1$desc)
data3 = left_join(data2, pca_complete1, by = 'id')
cluster_colours <- c("red", "blue", "green")
point_colours <- cluster_colours[data3$Cluster]
points(data3$lon, data3$lat, pch = 21, col = "black", bg = point_colours)


# Add the custom legend
par(new = TRUE)
plot(1, xlab = "", ylab = "", axes = FALSE, type = "n", xlim = c(-200, 200), ylim = c(-20, 20))
image(vx, vy, vz, col = rev(custom_palette(100)), axes = FALSE, xlab = "", ylab = "", add = TRUE)
polygon(c(185, 195, 195, 185), c(-7.8, -7.8, 10, 10))
minValue(smoothed_raster)
maxValue(smoothed_raster)
axis(4, at = c(5.4, 0.4, -4.6), labels = c(5, 10, 15), las = 1)




