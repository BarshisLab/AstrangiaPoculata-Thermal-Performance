#Author: Hannah Aichelman and Dan Barshis
#Last updated: February 2018
#This script plots the hobo logger data from the physiology ramp experiments to compare variation in temperature ramps between experiment days. 

#set working directory to location of hobo files

setwd("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/HoboLoggers_ForPlots")

#library packages needed for the code
library(ggplot2)
library(chron)
library(RColorBrewer)

par(mfrow=c(2,2))
### Cold Holobiont Ramps ###
#read in all cleaned hobo files for heat ramps
setwd("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/HoboLoggers_ForPlots/ColdHolobiont")

hoboCold1<-read.delim("ColdRamp07122017_clean.txt")
hoboCold1$Time<-sapply(strsplit(as.character(hoboCold1$DateTime), " "),"[[", 2)
hoboCold1$Time<-strptime(hoboCold1$Time, format="%R")
hoboCold1$DateTime<-strptime(hoboCold1$DateTime, format="%m/%d/%y %H:%M:%S")
colnames(hoboCold1)[3] <- "Temp"
head(hoboCold1)
hoboCold1$Temp_corr <- hoboCold1$Temp - 1.5 #correct based on glass thermometer calibration

hoboCold2<-read.delim("ColdRamp07132017_clean.txt")
hoboCold2$Time<-sapply(strsplit(as.character(hoboCold2$DateTime), " "),"[[", 2)
hoboCold2$Time<-strptime(hoboCold2$Time, format="%R")
hoboCold2$DateTime<-strptime(hoboCold2$DateTime, format="%m/%d/%y %H:%M:%S")
colnames(hoboCold2)[3] <- "Temp"
head(hoboCold2)
hoboCold2$Temp_corr <- hoboCold2$Temp - 1.5 #correct based on glass thermometer calibration

hoboCold3<-read.delim("ColdRamp07142017_clean.txt")
hoboCold3$Time<-sapply(strsplit(as.character(hoboCold3$DateTime), " "),"[[", 2)
hoboCold3$Time<-strptime(hoboCold3$Time, format="%R")
hoboCold3$DateTime<-strptime(hoboCold3$DateTime, format="%m/%d/%y %H:%M:%S")
colnames(hoboCold3)[3] <- "Temp"
head(hoboCold3)
hoboCold3$Temp_corr <- hoboCold3$Temp - 1.5 #correct based on glass thermometer calibration

hoboCold4<-read.delim("ColdRamp07162017_clean.txt")
hoboCold4$Time<-sapply(strsplit(as.character(hoboCold4$DateTime), " "),"[[", 2)
hoboCold4$Time<-strptime(hoboCold4$Time, format="%R")
hoboCold4$DateTime<-strptime(hoboCold4$DateTime, format="%m/%d/%y %H:%M:%S")
colnames(hoboCold4)[3] <- "Temp"
head(hoboCold4)
hoboCold4$Temp_corr <- hoboCold4$Temp - 1.5 #correct based on glass thermometer calibration


#USE THIS ONE:
pdf("Cold Ramp Experimental Tank Temperature by Day.pdf")
plot(hoboCold1$Time,hoboCold1$Temp_corr, type="l", lwd=4, col="#A6611A", ylab="Tank Temp (°C)", xlab="Time", ylim=c(4,20), yaxt = "n")
axis(2,at=c(0,6,9,12,15,18))
points(hoboCold2$Time, hoboCold2$Temp_corr, type="l", col="#DFC27D", lwd=4)
points(hoboCold3$Time, hoboCold3$Temp_corr, type="l", col="#80CDC1", lwd=4)
points(hoboCold4$Time, hoboCold4$Temp_corr, type="l", col="#018571", lwd=4)
legend("topright", c("Cold Ramp 1", "Cold Ramp 2", "Cold Ramp 3", "Cold Ramp 4"), lty=1, lwd=4, col=c("#A6611A", "#DFC27D", "#80CDC1","#018571"))
box()
dev.off()


### Heat Holobiont Ramps ###
setwd("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/HoboLoggers_ForPlots/HeatHolobiont")

#read in all cleaned hobo files for heat ramps
hoboHeat1<-read.delim("HeatRamp07072017_clean.txt")
hoboHeat1$Time<-sapply(strsplit(as.character(hoboHeat1$DateTime), " "),"[[", 2)
hoboHeat1$Time<-strptime(hoboHeat1$Time, format="%R")
hoboHeat1$DateTime<-strptime(hoboHeat1$DateTime, format="%m/%d/%y %H:%M:%S")
colnames(hoboHeat1)[3] <- "Temp"
head(hoboHeat1)
hoboHeat1$Temp_corr <- hoboHeat1$Temp - 1.377 #correct based on glass thermometer calibration

hoboHeat2<-read.delim("HeatRamp07082017_clean.txt")
hoboHeat2$Time<-sapply(strsplit(as.character(hoboHeat2$DateTime), " "),"[[", 2)
hoboHeat2$Time<-strptime(hoboHeat2$Time, format="%R")
hoboHeat2$DateTime<-strptime(hoboHeat2$DateTime, format="%m/%d/%y %H:%M:%S")
colnames(hoboHeat2)[3] <- "Temp"
head(hoboHeat2)
hoboHeat2$Temp_corr <- hoboHeat2$Temp - 1.377 #correct based on glass thermometer calibration

hoboHeat3<-read.delim("HeatRamp07102017_clean.txt")
hoboHeat3$Time<-sapply(strsplit(as.character(hoboHeat3$DateTime), " "),"[[", 2)
hoboHeat3$Time<-strptime(hoboHeat3$Time, format="%R")
hoboHeat3$DateTime<-strptime(hoboHeat3$DateTime, format="%m/%d/%y %H:%M:%S")
colnames(hoboHeat3)[3] <- "Temp"
head(hoboHeat3)
hoboHeat3$Temp_corr <- hoboHeat3$Temp - 1.377 #correct based on glass thermometer calibration

hoboHeat4<-read.delim("HeatRamp07112017_clean.txt")
hoboHeat4$Time<-sapply(strsplit(as.character(hoboHeat4$DateTime), " "),"[[", 2)
hoboHeat4$Time<-strptime(hoboHeat4$Time, format="%R")
hoboHeat4$DateTime<-strptime(hoboHeat4$DateTime, format="%m/%d/%y %H:%M:%S")
colnames(hoboHeat4)[3] <- "Temp"
head(hoboHeat4)
hoboHeat4$Temp_corr <- hoboHeat4$Temp - 1.377 #correct based on glass thermometer calibration


#USE THIS ONE:
pdf("Heat Ramp Experimental Tank Temperature by Day.pdf")
plot(hoboHeat1$Time, hoboHeat1$Temp_corr, type="l", lwd=4, col="#A6611A", ylab="Tank Temp (°C)", xlab="Time", ylim=c(18,34), yaxt = "n")
axis(2,at=c(18,22,26,29,32))
points(hoboHeat2$Time, hoboHeat2$Temp_corr, type="l", col="#DFC27D", lwd=4)
points(hoboHeat3$Time, hoboHeat3$Temp_corr, type="l", col="#80CDC1", lwd=4)
points(hoboHeat4$Time, hoboHeat4$Temp_corr, type="l", col="#018571", lwd=4)
legend("topleft", c("Heat Ramp 1", "Heat Ramp 2", "Heat Ramp 3", "Heat Ramp 4"), lty=1, lwd=4, col=c("#A6611A", "#DFC27D", "#80CDC1","#018571"))
dev.off()


################################################################################################
## Skeleton Analysis Starts Here ##
### Heat Skeleton Ramps ###
# set working directory to where heat skeleton cleaned files are saved
setwd("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/HoboLoggers_ForPlots/HeatRampSkeleton")

#read in all cleaned hobo files for heat ramps
hoboHeat1<-read.delim("HeatRampSkeleton07172017_clean.txt")
hoboHeat1$Time<-sapply(strsplit(as.character(hoboHeat1$DateTime), " "),"[[", 2)
hoboHeat1$Time<-strptime(hoboHeat1$Time, format="%R")
hoboHeat1$DateTime<-strptime(hoboHeat1$DateTime, format="%m/%d/%y %H:%M:%S")
colnames(hoboHeat1)[3] <- "Temp"
head(hoboHeat1)
hoboHeat1$Temp_corr <- hoboHeat1$Temp - 1.472 #correct based on glass thermometer calibration

hoboHeat2<-read.delim("HeatRampSkeleton07182017_clean.txt")
hoboHeat2$Time<-sapply(strsplit(as.character(hoboHeat2$DateTime), " "),"[[", 2)
hoboHeat2$Time<-strptime(hoboHeat2$Time, format="%R")
hoboHeat2$DateTime<-strptime(hoboHeat2$DateTime, format="%m/%d/%y %H:%M:%S")
colnames(hoboHeat2)[3] <- "Temp"
head(hoboHeat2)
hoboHeat2$Temp_corr <- hoboHeat2$Temp - 1.472 #correct based on glass thermometer calibration

hoboHeat3<-read.delim("HeatRampSkeleton07192017_clean.txt")
hoboHeat3$Time<-sapply(strsplit(as.character(hoboHeat3$DateTime), " "),"[[", 2)
hoboHeat3$Time<-strptime(hoboHeat3$Time, format="%R")
hoboHeat3$DateTime<-strptime(hoboHeat3$DateTime, format="%m/%d/%y %H:%M:%S")
colnames(hoboHeat3)[3] <- "Temp"
head(hoboHeat3)
hoboHeat3$Temp_corr <- hoboHeat3$Temp - 1.472 #correct based on glass thermometer calibration

hoboHeat4<-read.delim("HeatRampSkeleton07212017_clean.txt")
hoboHeat4$Time<-sapply(strsplit(as.character(hoboHeat4$DateTime), " "),"[[", 2)
hoboHeat4$Time<-strptime(hoboHeat4$Time, format="%R")
hoboHeat4$DateTime<-strptime(hoboHeat4$DateTime, format="%m/%d/%y %H:%M:%S")
colnames(hoboHeat4)[3] <- "Temp"
head(hoboHeat4)
hoboHeat4$Temp_corr <- hoboHeat4$Temp - 1.472 #correct based on glass thermometer calibration


#USE THIS ONE:
pdf("Skeleton Heat Ramp Experimental Tank Temperature by Day.pdf")
plot(hoboHeat1$Time, hoboHeat1$Temp_corr, type="l", lwd=4, col="#A6611A", ylab="Tank Temp (°C)", xlab="Time", ylim=c(18,34), yaxt="n")
axis(2,at=c(18,22,26,29,32))
points(hoboHeat2$Time, hoboHeat2$Temp_corr, type="l", col="#DFC27D", lwd=4)
points(hoboHeat3$Time, hoboHeat3$Temp_corr, type="l", col="#80CDC1", lwd=4)
points(hoboHeat4$Time, hoboHeat4$Temp_corr, type="l", col="#018571", lwd=4)
legend("topleft", c("Heat Skeleton Ramp 1", "Heat Skeleton Ramp 2", "Heat Skeleton Ramp 3", "Heat Skeleton Ramp 4"), lty=1, lwd=4, col=c("#A6611A", "#DFC27D", "#80CDC1","#018571"))
box()
dev.off()

### Cold Skeleton Ramps ###
# set working directory to where cold skeleton cleaned files are saved
setwd("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/HoboLoggers_ForPlots/ColdRampSkeleton")

#read in all cleaned hobo files for cold skeleton ramps
hoboCold1<-read.delim("ColdRampSkeleton07222017_clean.txt")
hoboCold1$Time<-sapply(strsplit(as.character(hoboCold1$DateTime), " "),"[[", 2)
hoboCold1$Time<-strptime(hoboCold1$Time, format="%R")
hoboCold1$DateTime<-strptime(hoboCold1$DateTime, format="%m/%d/%y %H:%M:%S")
colnames(hoboCold1)[3] <- "Temp"
head(hoboCold1)
hoboCold1$Temp_corr <- hoboCold1$Temp - 1.42 #correct based on glass thermometer calibration

hoboCold2<-read.delim("ColdRampSkeleton07242017_clean.txt")
hoboCold2$Time<-sapply(strsplit(as.character(hoboCold2$DateTime), " "),"[[", 2)
hoboCold2$Time<-strptime(hoboCold2$Time, format="%R")
hoboCold2$DateTime<-strptime(hoboCold2$DateTime, format="%m/%d/%y %H:%M:%S")
colnames(hoboCold2)[3] <- "Temp"
head(hoboCold2)
hoboCold2$Temp_corr <- hoboCold2$Temp - 1.42 #correct based on glass thermometer calibration

hoboCold3<-read.delim("ColdRampSkeleton07252017_clean.txt")
hoboCold3$Time<-sapply(strsplit(as.character(hoboCold3$DateTime), " "),"[[", 2)
hoboCold3$Time<-strptime(hoboCold3$Time, format="%R")
hoboCold3$DateTime<-strptime(hoboCold3$DateTime, format="%m/%d/%y %H:%M:%S")
colnames(hoboCold3)[3] <- "Temp"
head(hoboCold3)
hoboCold3$Temp_corr <- hoboCold3$Temp - 1.42 #correct based on glass thermometer calibration

hoboCold4<-read.delim("ColdRampSkeleton07262017_clean.txt")
hoboCold4$Time<-sapply(strsplit(as.character(hoboCold4$DateTime), " "),"[[", 2)
hoboCold4$Time<-strptime(hoboCold4$Time, format="%R")
hoboCold4$DateTime<-strptime(hoboCold4$DateTime, format="%m/%d/%y %H:%M:%S")
colnames(hoboCold4)[3] <- "Temp"
head(hoboCold4)
hoboCold4$Temp_corr <- hoboCold4$Temp - 1.42 #correct based on glass thermometer calibration


#USE THIS ONE:
pdf("Skeleton Cold Ramp Experimental Tank Temperature by Day.pdf")
plot(hoboCold1$Time,hoboCold1$Temp_corr, type="l", lwd=4, col="#A6611A", ylab="Tank Temp (°C)", xlab="Time", ylim=c(4,20), yaxt="n")
axis(2,at=c(0,6,9,12,15,18))
points(hoboCold2$Time, hoboCold2$Temp_corr, type="l", col="#DFC27D", lwd=4)
points(hoboCold3$Time, hoboCold3$Temp_corr, type="l", col="#80CDC1", lwd=4)
points(hoboCold4$Time, hoboCold4$Temp_corr, type="l", col="#018571", lwd=4)
legend("topright", c("Skeleton Cold Ramp 1", "Skeleton Cold Ramp 2", "Skeleton Cold Ramp 3", "Skeleton Cold Ramp 4"), lty=1, lwd=4, col=c("#A6611A", "#DFC27D", "#80CDC1","#018571"))
dev.off()


