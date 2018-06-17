#####################################################################################
# Author: Hannah E Aichelman
# Last updated: February 28, 2017
# This script reads in all VA and RI Hobo temperature loggers, calculates the % of time spent between a range of 2C increment temperatures, and plots the VA vs. RI comparison (Figure 2B)
#####################################################################################

# set working directory
setwd("~/ODU_MS/ODU_MS/FieldData/TemperatureData/ForPlots/ForTempCalculationPlots")
getwd()

# library needed packages
library(ggplot2)


#Read in first VA logger data
hobo<-read.delim("Tower_reef_1.txt")
hobo$DateTime<-strptime(hobo$DateTime, format="%m/%d/%y %H:%M:%S")
colnames(hobo)[1] <- "RecNum"
hobo$TempCorr<-hobo$Temp-0
head(hobo)
#Read in second VA logger data
hobo2<-read.delim("tower_reef_3_clean.txt")
hobo2$DateTime<-strptime(hobo2$DateTime, format="%m/%d/%y %H:%M")
colnames(hobo2)[1] <- "RecNum"
hobo2$TempCorr<-hobo2$Temp-0.232    # calibration 2 read +0.232 when calibrated in ice
head(hobo2)
#Read in third VA logger data
hobo3<-read.delim("TR_052717_1_clean.txt")
hobo3$DateTime<-strptime(hobo3$DateTime, format="%m/%d/%y %H:%M:%S")
colnames(hobo3)[3] <- "Temp"
hobo3$TempCorr<-hobo3$Temp-0.121       #calibration 1 read +0.121 when calibrated in ice
head(hobo3)

#Read in first logger data from RI, updated to include only days included in VA temp logger data
hobosean<-read.delim("FtWetherill_816_717_seangrace_v3.txt")
hobosean$DateTime<-strptime(hobosean$DateTime, format="%m/%d/%y %H:%M:%S")
colnames(hobosean)[1] <- "DateTime"
hobosean$RecNum <- seq.int(nrow(hobosean))
head(hobosean)
hobosean <- data.frame(hobosean$RecNum,hobosean$DateTime,hobosean$Temp)
colnames(hobosean) <- c("RecNum","DateTime","Temp")
#Read in Sean's second logger data from RI, updated to include only days included in VA temp logger data
hobosean2<-read.delim("Fort_Wetherill_Wall1_619-1107_seangrace_v3.txt")
hobosean2$DateTime<-strptime(hobosean2$DateTime, format="%m/%d/%y %H:%M:%S")
colnames(hobosean2)[1] <- "RecNum"
head(hobosean2)
hobosean2 <- data.frame(hobosean2$RecNum,hobosean2$DateTime,hobosean2$Temp)
colnames(hobosean2) <- c("RecNum","DateTime","Temp")
###########################################

#Combine all VA temp data and figure out total time loggers recorded temp
allVA <- rbind(hobo,hobo2,hobo3)
head(allVA)
#find max, min, mean, sd temperatures in all VA hobo logger data
max(allVA$TempCorr, na.rm=TRUE)
min(allVA$TempCorr, na.rm=TRUE)
mean(allVA$TempCorr, na.rm=TRUE)
sd(allVA$TempCorr, na.rm=TRUE)

# > which.min(allVA$Temp)
# [1] 7242
# > allVA[7242,]
     # RecNum            DateTime  Temp
# 7242   1790 2017-01-25 03:15:00 7.481
# > which.max(allVA$Temp)
# [1] 1058
# > allVA[1058,]
     # RecNum            DateTime   Temp
# 1058   1058 2016-09-03 12:15:00 24.992


#subset summer months to find average summer temp for VA data
VAAug <- subset(allVA, grepl("*-08-[:digit:]*", allVA$DateTime))
VAJun <- subset(allVA, grepl("*-06-[:digit:]*", allVA$DateTime))
VAJul <- subset(allVA, grepl("*-07-[:digit:]*", allVA$DateTime))
allVASummer <- rbind(VAJun, VAJul, VAAug)
mean(allVASummer$TempCorr)
sd(allVASummer$TempCorr)

countVA <- nrow(allVA) #one row = 15 minutes worth of temp data, but assuming that the temp is constant in those 15 mins
totalminsVA <- countVA*15 #loggers recording temp every 15 mins in VA


#plot a histogram of counts of VA temp, specify breaks
breaks <- c(0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32)
#breaks <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32)
VAhist <- hist(allVA$TempCorr, breaks=breaks)
VAhist$totaltime <- (VAhist$counts)*15 #multiply counts of VA histogram by 15, loggers recording every 15 minutes
VAhist$percenttime <- (VAhist$totaltime/totalminsVA)*100
sum(VAhist$percenttime) #check that all add up to 100%
#hist(VAhist$percenttime, breaks=breaks)
#plot VA histogram data using barplot() function
bp <- barplot(VAhist$percenttime, col="black", names.arg=breaks, border=TRUE)

#Had to keep RI temp data separate because they are not recording temperature at the same interval
RI1 <- hobosean
RI2 <- hobosean2
head(RI1)
head(RI2)
countRI1 <- nrow(RI1)
countRI2 <- nrow(RI2)

allRI <- rbind(hobosean,hobosean2)
head(allRI)
mean(allRI$Temp)
sd(allRI$Temp)

#determine max and min values for the RI hobo logger data
max(RI1$Temp, na.rm=TRUE)
max(RI2$Temp, na.rm=TRUE)
min(RI1$Temp, na.rm=TRUE)
min(RI2$Temp, na.rm=TRUE)

# > which.min(RI1$Temp)
# [1] 8748
# > RI1[8748,]
     # RecNum            DateTime Temp
# 8748   8748 2017-02-14 04:40:00 3.13
# > which.max(RI1$Temp)
# [1] 509
# > RI1[509,]
    # RecNum            DateTime  Temp
# 509    509 2016-08-28 15:32:00 23.15

#subset summer months to find average summer temp for VA data
RIAug <- subset(allRI, grepl("*-08-[:digit:]*", allRI$DateTime))
RIJun <- subset(allRI, grepl("*-06-[:digit:]*", allRI$DateTime))
RIJul <- subset(allRI, grepl("*-07-[:digit:]*", allRI$DateTime))
allRISummer <- rbind(RIJun, RIJul, RIAug)
mean(allRISummer$Temp)
sd(allRISummer$Temp)

totalminsRI1 <- countRI1*16 #loggers recording temp values every 16 mins in first RI logger
totalminsRI2 <- countRI2*15 #loggers recording temp values every 15 mins in second RI logger
totalminsRI = totalminsRI1+totalminsRI2

#First RI logger, every 16 mins
RI1hist <- hist(RI1$Temp, breaks=breaks)
RI1hist$totaltime <- (RI1hist$counts)*16 #multiply counts of first RI histogram by 16, loggers recording every 16 minutes
RI1hist$percenttime <- (RI1hist$totaltime/totalminsRI)*100
#sum(RI1hist$percenttime) #check that all add up to 100%
hist(RI1hist$percenttime, breaks=breaks)
#plot VA histogram data using barplot() function
bp <- barplot(RI1hist$percenttime, col="black", names.arg=breaks, border=TRUE)

#Second RI logger, every 15 mins
RI2hist <- hist(RI2$Temp, breaks=breaks)
RI2hist$totaltime <- (RI2hist$counts)*15 #multiply counts of second RI histogram by 15, loggers recording every 16 minutes
RI2hist$percenttime <- (RI2hist$totaltime/totalminsRI)*100
#sum(RI1percenttime) #check that all add up to 100%
hist(RI2hist $percenttime, breaks=breaks)
#plot VA histogram data using barplot() function
bp <- barplot(RI2hist $percenttime, col="black", names.arg=breaks, border=TRUE)

#add together percent time spent between breaks for each RI logger
RIpercenttime <- RI1hist$percenttime + RI2hist$percenttime

#create new data frame with VA and RI percent time spent at each break in temperature
PercentTime_data <- data.frame(breaks,VAhist$percenttime,RIpercenttime) #this doesn't work because breaks is different size than the VA and RI histograms. Ended up writing it manually into a .csv file with the name "PercentTime.csv"
write.csv(PercentTime_data, "~/ODU_MS/ODU_MS/FieldData/TemperatureData/ForPlots/ForTempCalculationPlots/PercentTime.csv")

#had to go in manually and change the data to a format that ggplot can understand. 
PercentTime <- read.csv("~/ODU_MS/ODU_MS/FieldData/TemperatureData/ForPlots/ForTempCalculationPlots/PercentTime.csv")
head(PercentTime)
  # X breaks state   percent
# 1 1      2    VA  0.000000
# 2 2      4    VA  0.000000
# 3 3      6    VA  0.000000
# 4 4      8    VA 12.984378
# 5 5     10    VA 20.592412
# 6 6     12    VA  4.334889

PercentTime$breaks <- factor(PercentTime$breaks, levels = c("0-2","2-4","4-6","6-8","8-10","10-12","12-14","14-16","16-18","18-20","20-22","22-24","24-26","26-28","28-30","30-32"))

Fig.Temp <- ggplot(data=PercentTime, aes(x=breaks, y=percent, fill=state)) +
geom_bar(stat="identity", colour = "black", position=position_dodge())+
xlab("Temperature (°C)") + ylab("Percent of Time Measured")+
scale_fill_manual(values=c("blue","red"))+
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(color = 'black', size=12, angle=90), #Set the text angle
  		axis.text.y = element_text(color = 'black', size=12),
  		axis.title = element_text(color='black', size=14),
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_rect(colour='grey95'),
        legend.background = element_rect(colour='black'),
        legend.position=c(.8,.8)) + #set legend location
  ggtitle("In Situ Temperature Environment") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 14, 
                                  hjust = 0))

Fig.Temp #view plot

ggsave(file="~/ODU_MS/ODU_MS/FieldData/TemperatureData/ForPlots/ForTempCalculationPlots/InSituTempEnviron_2Cbreaks.pdf", Fig.Temp, width = 6, height = 6, units = c("in"), useDingbats=FALSE)


#####################################################
#Old way of calculating time breaks and plotting percent time spent
#Count number of occurrences of temperature values within 5C intervals for VA
VA_0_5 <- sum(allVA$Temp >= 0 & allVA$Temp < 5)
VA_5_10 <- sum(allVA$Temp >= 5 & allVA$Temp < 10)
VA_10_15 <- sum(allVA$Temp >= 10 & allVA$Temp < 15)
VA_15_20 <- sum(allVA$Temp >= 15 & allVA$Temp < 20)
VA_20_25 <- sum(allVA$Temp >= 20 & allVA$Temp < 25)
VA_25_30 <- sum(allVA$Temp >= 25 & allVA$Temp < 30)
VA_30_35 <- sum(allVA$Temp >= 30 & allVA$Temp < 35)

#Count number of occurrences of temperature values within 5C intervals for first RI log
RI1_0_5 <- sum(RI1$Temp >= 0 & RI1$Temp < 5)
RI1_5_10 <- sum(RI1$Temp >= 5 & RI1$Temp < 10)
RI1_10_15 <- sum(RI1$Temp >= 10 & RI1$Temp < 15)
RI1_15_20 <- sum(RI1$Temp >= 15 & RI1$Temp < 20)
RI1_20_25 <- sum(RI1$Temp >= 20 & RI1$Temp < 25)
RI1_25_30 <- sum(RI1$Temp >= 25 & RI1$Temp < 30)
RI1_30_35 <- sum(RI1$Temp >= 30 & RI1$Temp < 35)

#Count number of occurrences of temperature values within 5C intervals for second RI log
RI2_0_5 <- sum(RI2$Temp >= 0 & RI2$Temp < 5)
RI2_5_10 <- sum(RI2$Temp >= 5 & RI2$Temp < 10)
RI2_10_15 <- sum(RI2$Temp >= 10 & RI2$Temp < 15)
RI2_15_20 <- sum(RI2$Temp >= 15 & RI2$Temp < 20)
RI2_20_25 <- sum(RI2$Temp >= 20 & RI2$Temp < 25)
RI2_25_30 <- sum(RI2$Temp >= 25 & RI2$Temp < 30)
RI2_30_35 <- sum(RI2$Temp >= 30 & RI2$Temp < 35)

#Figure out fraction of time spent at each temperature for VA and RI loggers6
VA1 <- (VA0_5*15)/totalminsVA
VA2 <- (VA5_10*15)/totalminsVA
VA3 <- (VA10_15*15)/totalminsVA
VA4 <- (VA15_20*15)/totalminsVA
VA5 <- (VA20_25*15)/totalminsVA
VA6 <- (VA25_30*15)/totalminsVA
VA7 <- (VA30_35*15)/totalminsVA
#VA1+VA2+VA3+VA4+VA5+VA6+VA7 #check they add up to 1, they do

RI1_1 <- (RI1_0_5*16)/totalminsRI
RI1_2 <- (RI1_5_10*16)/totalminsRI
RI1_3 <- (RI1_10_15*16)/totalminsRI
RI1_4 <- (RI1_15_20*16)/totalminsRI
RI1_5 <- (RI1_20_25*16)/totalminsRI
RI1_6 <- (RI1_25_30*16)/totalminsRI
RI1_7 <- (RI1_30_35*16)/totalminsRI

RI2_1 <- (RI2_0_5*15)/totalminsRI
RI2_2 <- (RI2_5_10*15)/totalminsRI
RI2_3 <- (RI2_10_15*15)/totalminsRI
RI2_4 <- (RI2_15_20*15)/totalminsRI
RI2_5 <- (RI2_20_25*15)/totalminsRI
RI2_6 <- (RI2_25_30*15)/totalminsRI
RI2_7 <- (RI2_30_35*15)/totalminsRI

RI1_1+RI1_2+RI1_3+RI1_4+RI1_5+RI1_6+RI1_7+RI2_1+RI2_2+RI2_3+RI2_4+RI2_5+RI2_6+RI2_7 #check they add up to 1, they do

#Combine RI percentages for each temperature range
RI1 <- RI1_1 + RI2_1
RI2 <- RI1_2 + RI2_2
RI3 <- RI1_3 + RI2_3
RI4 <- RI1_4 + RI2_4
RI5 <- RI1_5 + RI2_5
RI6 <- RI1_6 + RI2_6
RI7 <- RI1_7 + RI2_7
 
#Took all of the count and percent time data calculated above and copied into new spreadsheet called "TempCounts.txt"
counts<-read.delim("TempCounts_20180228.txt")
head(counts)

#rearrange temp range values so they plot in order
counts$temp <- factor(counts$temp, levels = c("0-5","5-10","10-15","15-20","20-25","25-30","30-35"))

#Now plot it...

Fig.Temp <- ggplot(data=counts, aes(x=temp, y=percent, fill=state)) +
geom_bar(stat="identity", colour = "black", position=position_dodge())+
xlab("Temperature (°C)") + ylab("Percent of Time Measured")+
scale_fill_manual(values=c("skyblue2","tomato1"))+
#scale_fill_manual(values=c("black","grey"))+
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(color = 'black', size=12), #Set the text angle
  		axis.text.y = element_text(color = 'black', size=12),
  		axis.title = element_text(color='black', size=14),
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_rect(colour='grey95'),
        legend.background = element_rect(colour='black'),
        legend.position=c(.8,.8)) + #set legend location
  ggtitle("In Situ Temperature Environment") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 14, 
                                  hjust = 0))

Fig.Temp #view plot

