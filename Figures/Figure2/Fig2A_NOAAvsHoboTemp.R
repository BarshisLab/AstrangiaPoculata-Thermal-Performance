#This script reads in NOAA buoy data and Hobo temperature logger data and plots both over the time considered 
#Written by D. Barshis and H. Aichelman
#Last updated May 2018

setwd("~/ODU_MS/ODU_MS/FieldData/TemperatureData/ForPlots")

getwd()

# read in VA noaa buoy data
boing<-read.delim("2016_Mar-2017_Sep_VA.txt", na.strings='999')
boing<-cbind("datetime"=paste(paste(boing$X.YY..yr, boing$MM.mo, boing$DD.dy, sep="-")," ",paste(boing$hh.hr,boing$mm.mn, sep=":"),sep=""),boing)
boing$datetime<-strptime(boing$datetime, format="%Y-%m-%d %H:%M")
boing<-boing[order(boing$datetime),]

####Note the na.strings='999' argument
# read in RI noaa buoy data
boing2<-read.delim("2016_Mar-2017_Sep_Newport_NWPR1.txt", na.strings='999')
boing2<-cbind("datetime"=paste(paste(boing2$X.YY..yr, boing2$MM.mo, boing2$DD.dy, sep="-")," ",paste(boing2$hh.hr,boing2$mm.mn, sep=":"),sep=""),boing2)
boing2$datetime<-strptime(boing2$datetime, format="%Y-%m-%d %H:%M")
boing2<-boing2[order(boing2$datetime),]
boing2$WTMP.degC[boing2$WTMP.degC==999]<-NA

#read in hobo loggers and format DateTime
hobo<-read.delim("Tower_reef_1.txt")
hobo$DateTime<-strptime(hobo$DateTime, format="%m/%d/%y %H:%M:%S")
hobo2<-read.delim("tower_reef_3_clean.txt") #SN 20014300, Calibration 2
hobo2$DateTime<-strptime(hobo2$DateTime, format="%m/%d/%y %H:%M")
hobo2$TempCorr<-hobo2$Temp-0.232    # calibration 2 read +0.232 when calibrated in ice
hobo3<-read.delim("TR_052717_1_clean.txt") #SN 20014302, Calibration 1
hobo3$DateTime<-strptime(hobo3$DateTime, format="%m/%d/%y %H:%M:%S")
hobo3$TempCorr<-hobo3$TR_052717_1-0.121       #calibration 1 read +0.121 when calibrated in ice
hobosean<-read.delim("FtWetherill_816_717_seangrace_v2.txt")
hobosean$DateTime<-strptime(hobosean$DateTime, format="%m/%d/%y %H:%M:%S")
hobosean2<-read.delim("Fort_Wetherill_Wall1_619-1107_seangrace_v2.txt")
hobosean2$DateTime<-strptime(hobosean2$DateTime, format="%m/%d/%y %H:%M:%S")

#Plot NOAA and Hobo temp data over time considered at both VA and RI sites
#USE THIS ONE (includes CORRECTED VA loggers, updated x-axis labels and plot of Topt lines)
pdf("2016-2017_VA_vs_RI_NOAA-Hobo_corrVA.pdf")
plot(boing$datetime,boing$WTMP.degC, type="l", ylab="Temperature (°C)", main="2016-2017 Water Temp VA vs RI", ylim=c(0,30), xaxt='n', xlab='')
points(boing2$datetime,boing2$WTMP.degC, type="l", col="grey")
points(hobo$DateTime, hobo$Temp, type="l", col="red")
points(hobo2$DateTime, hobo2$TempCorr, type="l", col="red")
points(hobo3$DateTime, hobo3$TempCorr, type="l", col="red")
points(hobosean$DateTime, hobosean$Temp, type="l", col="blue")
points(hobosean2$DateTime, hobosean2$Temp, type="l", col="blue")
#abline(h=23.04, col="blue", lwd=3) #RI-brown Topt average
#abline(h=21.35, col="blue", lwd=3, lty=2) #RI-white Topt average
#abline(h=26.81, col="red", lwd=3) #VA-brown Topt average
#abline(h=28.23, col="red", lwd=3, lty=2) #VA-white Topt average
ticks<-seq(from=boing$datetime[1], to=as.POSIXct("2017-10-01 00:28:00 EST"), by='1 month')
lbl<-strftime(ticks, format="%b-%y")
axis(side=1, at=ticks, labels=F)
text(ticks, par("usr")[3] - 1, srt=60, adj=1 , labels=lbl, xpd = T, cex=1)
legend("bottomright", c("NOAA-VA", "hobo-VA", "NOAA-RI", "hobo-RI"), lty=1, lwd=3, col=c("black",  "red", "grey", "blue"), )
dev.off()

#older version of the plot, does not have the angled month-year date labels
pdf("2016-2017_VA_vs_RI_NOAA-Hobo.pdf")
plot(boing$datetime,boing$WTMP.degC, type="l", ylab="Water Temp C", xlab="Date:Time", main="2016-2017 Water Temp VA vs RI", ylim=c(0,30))
#abline(h=13.3)
points(boing2$datetime,boing2$WTMP.degC, type="l", col="grey")
points(hobo$DateTime, hobo$Temp, type="l", col="red")
points(hobo2$DateTime, hobo2$Temp, type="l", col="red")
points(hobo3$DateTime, hobo3$TR_052717_1, type="l", col="red")
points(hobosean$DateTime, hobosean$Temp, type="l", col="blue")
legend("topleft", c("NOAA-VA", "hobo-VA", "NOAA-RI", "hobo-RI"), lty=1, lwd=3, col=c("black",  "red", "grey", "blue"), )
dev.off()

# to plot only noaa buoy data:
pdf("2016-2017_VA_vs_RI_NOAA.pdf")
plot(boing$datetime,boing$WTMP.degC, type="l", ylab="Temperature (°C)", main="2016-2017 Water Temp VA vs RI", ylim=c(0,30), xaxt='n', xlab='')
points(boing2$datetime,boing2$WTMP.degC, type="l", col="grey")
ticks<-seq(from=boing$datetime[1], to=as.POSIXct("2017-10-01 00:28:00 EST"), by='1 month')
lbl<-strftime(ticks, format="%b-%y")
axis(side=1, at=ticks, labels=F)
text(ticks, par("usr")[3] - 1, srt=60, adj=1 , labels=lbl, xpd = T, cex=1)
legend("topleft", c("NOAA-VA", "NOAA-RI"), lty=1, lwd=3, col=c("black", "grey"), )
dev.off()



# some additional temp calculations

head(boing)
head(boing2)
mean(boing$WTMP.degC, na.rm=TRUE)
sd(boing$WTMP.degC, na.rm=TRUE)

mean(boing2$WTMP.degC, na.rm=TRUE)
sd(boing2$WTMP.degC, na.rm=TRUE)

t.test(boing$WTMP.degC,boing2$WTMP.degC, paired=TRUE, conf.level=0.95) #can't do this because the data frames are different sizes
ks.test(boing$WTMP.degC,boing2$WTMP.degC)

	# Two-sample Kolmogorov-Smirnov test

# data:  boing$WTMP.degC and boing2$WTMP.degC
# D = 0.29942, p-value < 2.2e-16
# alternative hypothesis: two-sided