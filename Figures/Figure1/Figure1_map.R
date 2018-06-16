# Original code written by JP Rippe (jpr6mg@gmail.com)
# Updated by Hannah Aichelman March 13, 2018
# This script creates the zoomed in and inset map for figure 1. These images were then put together using Adobe Illustrator.

install.packages("gpclib", type="source")
library(mapproj)
library(dplyr)
library(ggplot2)
library(ggmap)
detach("package:ggmap", unload=TRUE)
library(rgeos)
library(rgdal)
library(maps)
library(mapdata)
library(maptools)

setwd("~/Documents/Manuscripts/AstrangiaPhysiology/Map_Figure")

#clip gshhs_f.b to the broader area that you want to include in your map. gshhs_f.b is the high definition noaa coastline layer
if (!rgeosStatus()) gpclibPermit()
gshhs.f.b <- "gshhs_f.b"
#Crop global map layer to desired extent
sf1 <- getRgshhsMap(gshhs.f.b, xlim = c(-100, -55), ylim = c(20, 55)) %>%
  fortify()

#Read in coordinates of sampling sites
a=read.csv('GPSCoordinates.csv')
# a$reefZone <- c('or','or','ir','or','ir','ir','or','ir') #jp used this for florida map

#plot zoomed in map of virginia and rhode island collection sites
map_zoom <- ggplot() + 
  geom_polygon(data=sf1, aes(x=long, y=lat, group = group), fill = "grey70", color='black', lwd = 0.1)+
  geom_point(data=a[c(1:3)], aes(x=site_long, y=site_lat, shape = Site), size=3, col = 'black') +
  scale_shape_manual(values=c(17,16))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_fixed(ratio=1, xlim = c(-77,-69), ylim = c(35,43))+
  scale_y_continuous(breaks=seq(35,43, 1))
map_zoom

#save ggplot figure
ggsave(file="~/Documents/Manuscripts/AstrangiaPhysiology/Map_Figure/map_zoom.pdf", map_zoom, width = 8, height = 6, units = c("in"), useDingbats=FALSE)

#plot map further away to use for inset
map_distance <- ggplot() + 
  geom_polygon(data=sf1, aes(x=long, y=lat, group = group), fill = "grey70", color='black', lwd = 0.1)+
  #geom_point(data=a[c(1:3)], aes(x=site_long, y=site_lat, shape = Site), size=3, col = 'black') +
  #scale_shape_manual(values=c(15,17))+
  theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_fixed(ratio=1, xlim = c(-90,-60), ylim = c(25,50))+
  scale_y_continuous(breaks=seq(35,43, 1))
map_distance

#save figure
ggsave(file="~/Documents/Manuscripts/AstrangiaPhysiology/Map_Figure/map_distance.pdf", map_distance, width = 5, height = 5, units = c("in"), useDingbats=FALSE)


