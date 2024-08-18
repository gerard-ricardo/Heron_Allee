#2022 Winds at Heron

#setwd("C:/Users/g_ric/OneDrive/1 Work/4 Writing/1 Allee effects/allee effects")

#devtools::install_github("ropensci/dataaimsr", upgrade = 'never')  #

library(plyr)
library(dataaimsr)
library(tidyverse)
library(ggmap)

#see myfun for api

library(ggplot2)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2")  #set theme in code
#####windspeed##############
#Heron
data1 = aims_filter_values("weather", filter_name = "series")  #filter for weatehr
# heron_id = dplyr::filter(data1,
#                             grepl('Heron', series))  #get all heronree  IDs
# heron_wind_id <- filter(heron_id, grepl('Wind', series)) %>% 
#               filter(., grepl('Speed', series))  #159 (10 mins)  or 184 (30 mins)
# heron_wind_id
# 
# heron_wind <- aims_data("weather", api_key = my_api_key,
#                        filters = list(series_id = 186,
#                                       from_date = "2022-01-01",
#                                       thru_date = "2022-12-01"))  

###Heron
heron_id = dplyr::filter(data1, grepl('Heron', series))  #get all xxx   IDs
heron_wind_id <- filter(heron_id, grepl('Wind', series)) %>% filter(., grepl('Speed', series))  #186 looks good


#extract mean wind
heron_wind <- aims_data("weather", api_key = my_api_key,
                        filters = list(series_id = heron_wind_id$series_id[6],
                                       from_date = "2022-12-13",
                                       thru_date = "2022-12-16"))  

heron_wind$knots = heron_wind$qc_val/1.852   #change to knots
#corect time zone
tz = "Australia/Brisbane"
heron_wind$time = as.POSIXct(heron_wind$time, tz = tz)

heron_wind$t_date = substr(heron_wind$time, start = 1, stop = 10)  #seperates date
heron_wind$t_time = substr(heron_wind$time, start = 12, stop = 20)  #seperate time
heron_wind_spawn_n = subset(heron_wind, t_date == '2022-12-13' & t_time > '18:00:00' & t_time < '23:00:00')
heron_wind_spawn_n = heron_wind_spawn_n %>% select(time, t_date, t_time, knots)

(night_ave = median(heron_wind_spawn_n$knots))
sd(heron_wind_spawn_n$knots)

# Convert the time column to POSIXct with the correct time zone for Heron Island (Australia/Brisbane)
heron_wind_spawn_n$time <- as.POSIXct(heron_wind_spawn_n$time, tz = "UTC")

# Set the xmin and xmax with the same time zone
xmin <- as.POSIXct("2022-12-13 18:00:00", tz = "Australia/Brisbane")
xmax <- as.POSIXct("2022-12-13 23:00:00", tz = "Australia/Brisbane")



# Plotting
p4 = ggplot() + 
  geom_point(heron_wind, mapping = aes(time, knots), col = 'steelblue', alpha = 0.5) +
  annotate("rect", xmin = as.POSIXct("2022-12-13 18:00:00", tz = tz), xmax = as.POSIXct("2022-12-13 23:00:00", tz = tz), 
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "grey") + # add grey rectangle
  annotate("rect", xmin = as.POSIXct("2022-12-14 18:00:00", tz = tz), xmax = as.POSIXct("2022-12-14 23:00:00", tz = tz), 
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "grey") + # add grey rectangle
  annotate("rect", xmin = as.POSIXct("2022-12-15 18:00:00", tz = tz), xmax = as.POSIXct("2022-12-15 23:00:00", tz = tz), 
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "grey") +
  labs(x = "Date", y = "Knots") # add another rectangle for the next day
  # geom_vline(xintercept = as.POSIXct("2022-12-14 00:00:00", tz = tz), col = "red", linetype = "dashed") + # vertical line on 14th
  # geom_vline(xintercept = as.POSIXct("2022-12-15 00:00:00", tz = tz), col = "red", linetype = "dashed") + # vertical line on 15th
  # geom_vline(xintercept = as.POSIXct("2022-12-16 00:00:00", tz = tz), col = "red", linetype = "dashed")   # vertical line on 16th


ggsave(p4, filename = 'fig1_p4.pdf',  path = "./plots", device = 'pdf',  width = 7, height = 3)  #


#single day
ggplot() + 
  geom_point(heron_wind_spawn_n, mapping = aes(time, knots), col = 'steelblue', alpha = 0.5) +
  annotate("rect", xmin = as.POSIXct("2022-12-13 18:00:00", tz = tz), xmax = as.POSIXct("2022-12-13 23:00:00", tz = tz), 
           ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "grey")  # add grey rectangle





#####################

# #save(onet.wind, file = "onetwind10yr.Rdata")
# load('onetwind10yr.Rdata')
# 
# onet.wind$knots = onet.wind$qc_val/1.852   #change to knots
# 
# #ggplot()+geom_point(onet.wind, mapping = aes(time, knots), col = 'steelblue', alpha = 0.1) # wind across time
# 
# #filter the dates coarsely night 6pm to 6am
# onet.wind$t.time = substr(onet.wind$time, start = 12, stop = 20)  #extracts first letters
# #convert time to proportion
# onet.wind$j.time = sapply(strsplit(onet.wind$t.time,":"),
#                           function(x) {
#                             x <- as.numeric(x)
#                             (x[1]+x[2]/60 + x[3]/60)/24
#                           })  #split by : and then convert to prop
# 
# 
# ot.w.day = onet.wind %>% filter(j.time > 0.25, j.time < 0.75)  #filter
# ot.w.day$day = 'day'
# ot.w.night = anti_join(onet.wind, ot.w.day, by = 'j.time')
# ot.w.night$day = 'night'
# ot.w = rbind(ot.w.day, ot.w.night)
# 
# #does night have less wind than day?
# ot.w %>%  group_by(day) %>%  summarise(quantile(.$knots, probs = c(.01, 0.25, 0.5, 0.75, .99), 'na.rm' = T)) # quartile. median par
# #quantile(ot.w.night$knots, probs = c(.01, 0.25, 0.5, 0.75, .99), 'na.rm' = T) # nightime minutely weaker
# ggplot(ot.w, aes(x=as.factor(day), y=knots)) + 
#   geom_boxplot(fill="slateblue", alpha=0.2) + 
#   xlab("Time of day") 
# 
# #filter for summer (Oct-Dec) and then day/night
# ot.w$t.date = substr(ot.w$time, start = 1, stop = 10)  #extracts first letters
# ot.w$j.month = sapply(strsplit(ot.w$t.date,"-"),
#                       function(x) {
#                         x <- as.numeric(x)
#                         x[2]
#                       })  #split and eztract month
# head(ot.w)
# 
# ot.w.sum = ot.w %>% filter(j.month >= 10, j.time <= 12)  #filter
# tail(ot.w.sum)
# 
# diel = ot.w.sum %>%  split(.$day) 
# 
# diel$day$knots %>%     quantile(., probs = c(.01, 0.25, 0.5, 0.75, .99), 'na.rm' = T) # quartile. median par
# 
# diel$night$knots %>%     quantile(., probs = c(.01, 0.25, 0.5, 0.75, .99), 'na.rm' = T) # small diff less in summer
# 
# 
# ggplot(ot.w.sum, aes(x=as.factor(day), y=knots)) + 
#   geom_boxplot(fill="slateblue", alpha=0.2) + 
#   xlab("Time of day")+theme_sleek1()



#########wind direction##############################

heron.id.dir <- filter(heron.id,
                           grepl('Wind', series)) %>% filter(.,
                                                             grepl('Direction', series))  #


heron.wind.dir <- aims_data("weather", api_key = my_api_key,
                           filters = list(series_id = 188,
                                          from_date = "2022-12-13",
                                          thru_date = "2022-12-15"))  #last 10 years

#save(onet.wind.dir, file = "onetwind10yr.dir.Rdata")
#load('onetwind10yr.Rdata')

#range(onet.wind.dir$time)
#ggplot()+geom_point(onet.wind.dir, mapping = aes(time, qc_val), col = 'steelblue', alpha = 0.1)

#corect time zone
tz = "Australia/Brisbane"
heron.wind.dir$time = as.POSIXct(heron.wind.dir$time, tz = tz)

heron.wind.dir$t.date = substr(heron.wind.dir$time, start = 1, stop = 10)  #seperates date
heron.wind.dir$t.time = substr(heron.wind.dir$time, start = 12, stop = 20)  #seperate time
heron.wind.spawn.d = subset(heron.wind.dir, t.date %in% '2022-12-15')  #subset 14th
tail(heron.wind.spawn.d)
#subset for spawning period
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
#df1$ws = df1
pollutionRose(heron.windd.spawn.n1, pollutant = "qc_val")

pollutionRose(mydata, pollutant = "pm10", type = "year", statistic = "prop.mean")


# Example data frame
set.seed(123) # For reproducibility
wind_data <- data.frame(
  wind_speed = rnorm(100, mean = 10, sd = 3),  # Simulated wind speed data
  wind_dir = sample(0:360, 100, replace = TRUE)  # Simulated wind direction data
)
# Create a wind rose plot
windRose(wind_data, ws = "wind_speed", wd = "wind_dir", 
         main = "Wind Rose for Heron Island",
         key.title = "Wind Speed (m/s)",
         palette = "YlGnBu")



