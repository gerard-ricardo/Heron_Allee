


library(plyr)
library(dataaimsr)
library(tidyverse)
library(ggmap)


library(ggplot2)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2")  #set theme in code
data1 = aims_filter_values("weather", filter_name = "series")  #filter for weatehr

heron_id = dplyr::filter(data1, grepl('Heron', series))  #get all xxx   IDs
heron_wind_id <- filter(heron_id, grepl('Wind', series)) %>% filter(., grepl('Speed', series))  #186 looks good


heron_wind <- aims_data("weather", api_key = my_api_key,
                        filters = list(series_id = heron_wind_id$series_id[6],
                                       from_date = "2022-12-13",
                                       thru_date = "2022-12-16"))  

heron_wind$knots = heron_wind$qc_val/1.852   #change to knots
tz = "Australia/Brisbane"
heron_wind$time = as.POSIXct(heron_wind$time, tz = tz)

heron_wind$t_date = substr(heron_wind$time, start = 1, stop = 10)  #seperates date
heron_wind$t_time = substr(heron_wind$time, start = 12, stop = 20)  #seperate time
heron_wind_spawn_n = subset(heron_wind, t_date == '2022-12-13' & t_time > '18:00:00' & t_time < '23:00:00')
heron_wind_spawn_n = heron_wind_spawn_n %>% select(time, t_date, t_time, knots)

(night_ave = median(heron_wind_spawn_n$knots))
sd(heron_wind_spawn_n$knots)

heron_wind_spawn_n$time <- as.POSIXct(heron_wind_spawn_n$time, tz = "UTC")

xmin <- as.POSIXct("2022-12-13 18:00:00", tz = "Australia/Brisbane")
xmax <- as.POSIXct("2022-12-13 23:00:00", tz = "Australia/Brisbane")



p4 = ggplot() + 
  geom_point(heron_wind, mapping = aes(time, knots), col = 'steelblue', alpha = 0.5) +
  annotate("rect", xmin = as.POSIXct("2022-12-13 18:00:00", tz = tz), xmax = as.POSIXct("2022-12-13 23:00:00", tz = tz), 
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "grey") + # add grey rectangle
  annotate("rect", xmin = as.POSIXct("2022-12-14 18:00:00", tz = tz), xmax = as.POSIXct("2022-12-14 23:00:00", tz = tz), 
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "grey") + # add grey rectangle
  annotate("rect", xmin = as.POSIXct("2022-12-15 18:00:00", tz = tz), xmax = as.POSIXct("2022-12-15 23:00:00", tz = tz), 
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "grey") +
  labs(x = "Date", y = "Knots")+ # add another rectangle for the next day
scale_y_continuous(name = "Knots", sec.axis = sec_axis(~ . * 1.852, name = "km/h")) # Convert knots to km/h
p4

ggsave(p4, filename = 'fig1_p4.pdf',  path = "./plots", device = 'pdf',  width = 7, height = 3)  #


ggplot() + 
  geom_point(heron_wind_spawn_n, mapping = aes(time, knots), col = 'steelblue', alpha = 0.5) +
  annotate("rect", xmin = as.POSIXct("2022-12-13 18:00:00", tz = tz), xmax = as.POSIXct("2022-12-13 23:00:00", tz = tz), 
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "grey")  # add grey rectangle










heron.id.dir <- filter(heron.id,
                           grepl('Wind', series)) %>% filter(.,
                                                             grepl('Direction', series))  #


heron.wind.dir <- aims_data("weather", api_key = my_api_key,
                           filters = list(series_id = 188,
                                          from_date = "2022-12-13",
                                          thru_date = "2022-12-15"))  #last 10 years



tz = "Australia/Brisbane"
heron.wind.dir$time = as.POSIXct(heron.wind.dir$time, tz = tz)

heron.wind.dir$t.date = substr(heron.wind.dir$time, start = 1, stop = 10)  #seperates date
heron.wind.dir$t.time = substr(heron.wind.dir$time, start = 12, stop = 20)  #seperate time
heron.wind.spawn.d = subset(heron.wind.dir, t.date %in% '2022-12-15')  #subset 14th
tail(heron.wind.spawn.d)
heron.windd.spawn.n1 = subset(heron.wind.spawn.d, t.time > '17:30:00')  #remove factor treatment level. Use '%in%' to keep.
heron.windd.spawn.n1 = subset(heron.windd.spawn.n1, t.time < '23:30:00')  #remove factor treatment level. Use '%in%' to keep.
median(heron.windd.spawn.n1$qc_val, na.rm=T)
heron.windd.spawn.n1 %>% select(qc_val )
str(heron.windd.spawn.n1)

breaks = 20
ff = hist(heron.windd.spawn.n1$qc_val, breaks = breaks) #breaks it into density
df1 = data.frame(breaks = ff$breaks, counts = ff$counts) #bin counts

library(ggplot2)
ggplot(df1, aes(x = breaks, y = counts)) +
  coord_polar(theta = "x", start = -0.15) +
  geom_bar(stat = "identity") + theme_light()+
  scale_x_continuous(breaks = seq(0, 360, 60))

library(openair)
pollutionRose(heron.windd.spawn.n1, pollutant = "qc_val")

pollutionRose(mydata, pollutant = "pm10", type = "year", statistic = "prop.mean")


set.seed(123) # For reproducibility
wind_data <- data.frame(
  wind_speed = rnorm(100, mean = 10, sd = 3),  # Simulated wind speed data
  wind_dir = sample(0:360, 100, replace = TRUE)  # Simulated wind direction data
)
windRose(wind_data, ws = "wind_speed", wd = "wind_dir", 
         main = "Wind Rose for Heron Island",
         key.title = "Wind Speed (m/s)",
         palette = "YlGnBu")



