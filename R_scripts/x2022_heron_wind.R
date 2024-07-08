#2022 Winds at Heron

#setwd("C:/Users/g_ric/OneDrive/1 Work/4 Writing/1 Allee effects/allee effects")

#devtools::install_github("ropensci/dataaimsr", upgrade = 'never')  #

library(plyr)
library(dataaimsr)
library(tidyverse)
library(ggmap)
my_api_key <- "P1SNKfYxyW9dYjtpp05CP3wFCJMhm8CX3TfkgsB9"

library(ggplot2)
source("https://raw.githubusercontent.com/gerard-ricardo/data/master/theme_sleek2")  #set theme in code
#####windspeed##############
# #One Tree
# data1 = aims_filter_values("weather", filter_name = "series")  #filter for weatehr
# one.tree.id = dplyr::filter(data1,
#                             grepl('Tree', series))  #get all onetree  IDs
# onet.wind.id <- filter(one.tree.id,
#                        grepl('Wind', series)) %>% filter(.,
#                                                          grepl('Speed', series))  #159 (10 mins)  or 184 (30 mins)
# onet.wind.id
# 
# onet.wind <- aims_data("weather", api_key = my_api_key,
#                        filters = list(series_id = 184,
#                                       from_date = "2021-01-01",
#                                       thru_date = "2021-12-01"))  #last 10 years (yyyy/mm/dd)

###Heron
heron.id = dplyr::filter(data1, grepl('Heron', series))  #get all xxx   IDs
heron.wind.id <- filter(heron.id,
                        grepl('Wind', series)) %>% filter(.,
                                                          grepl('Speed', series))  #186 looks good
heron.wind.id$series_id

heron.wind <- aims_data("weather", api_key = my_api_key,
                        filters = list(series_id = heron.wind.id$series_id[2],
                                       from_date = "2022-12-10",
                                       thru_date = "2022-12-20"))  #last 10 years

heron.wind$knots = heron.wind$qc_val/1.852   #change to knots

heron.wind$t.date = substr(heron.wind$time, start = 1, stop = 10)  #seperates date
heron.wind$t.time = substr(heron.wind$time, start = 12, stop = 20)  #seperate time
heron.wind.spawn.d = subset(heron.wind, t.date %in% '2022-12-14')  #remove factor treatment level. Use '%in%' to keep.
#heron.wind.spawn.d = heron.wind
day.ave = mean(heron.wind.spawn.d$knots)
heron.wind.spawn.n = subset(heron.wind.spawn.d, t.time > '17:30:00')  #remove factor treatment level. Use '%in%' to keep.
heron.wind.spawn.n = subset(heron.wind.spawn.n, t.time < '23:30:00')  #remove factor treatment level. Use '%in%' to keep.

night.ave = mean(heron.wind.spawn.n$knots)
night.ave
sd(heron.wind.spawn.n$knots)

heron.wind$time
ggplot()+geom_point(heron.wind, mapping = aes(time, knots), col = 'steelblue', alpha = 0.1)

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
                                          from_date = "2022-12-10",
                                          thru_date = "2022-12-20"))  #last 10 years

#save(onet.wind.dir, file = "onetwind10yr.dir.Rdata")
#load('onetwind10yr.Rdata')

#range(onet.wind.dir$time)
#ggplot()+geom_point(onet.wind.dir, mapping = aes(time, qc_val), col = 'steelblue', alpha = 0.1)

heron.wind.dir$t.date = substr(heron.wind.dir$time, start = 1, stop = 10)  #seperates date
heron.wind.dir$t.time = substr(heron.wind.dir$time, start = 12, stop = 20)  #seperate time
heron.wind.spawn.d = subset(heron.wind.dir, t.date %in% '2022-12-14')  #remove factor treatment level. Use '%in%' to keep.
tail(heron.wind.spawn.d)
heron.windd.spawn.n1 = subset(heron.wind.spawn.d, t.time > '17:30:00')  #remove factor treatment level. Use '%in%' to keep.
heron.windd.spawn.n2 = subset(heron.windd.spawn.n1, t.time < '23:30:00')  #remove factor treatment level. Use '%in%' to keep.

mean(heron.windd.spawn.n2$qc_val)

# ff = hist(heron.windd.spawn.n$qc_val, breaks = 20) #breaks it into density
# df1 = data.frame(breaks = ff$breaks[1:18], counts = ff$counts) #bin counts
# 
# library(ggplot2)
# ggplot(df1, aes(x = breaks, y = counts)) +
#   coord_polar(theta = "x", start = -0.15) +
#   geom_bar(stat = "identity") + theme_light()+
#   scale_x_continuous(breaks = seq(0, 360, 60))
# 
# library(openair)
# #df1$ws = df1
# pollutionRose(df1, pollutant = "counts")

