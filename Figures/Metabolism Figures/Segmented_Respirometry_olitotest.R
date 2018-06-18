#Title: Photosynthesis and Respiration Calculations
#Author: HM Putnam, updated by Hannah Aichelman to analyze the data for the Astrangia respirometry experiments
# original script found at the following page: https://github.com/hputnam/Mcap_PGA_TGA/blob/master/RAnalysis/Scripts/Respirometry.R
#Date Last Modified: October 2017
#See Readme file for details (on H Putnam github page)

rm(list=ls()) #clears workspace 

#Read in required libraries
##### Include Versions of libraries

#library(devtools)
#install_github('colin-olito/LoLinR')

library("ggplot2")
library("segmented")
library("plotrix")
library("gridExtra")
library("LoLinR")
library("lubridate")
library("chron")
library("Hmisc")

#Required Data files

# Set Working Directory:
setwd("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/20170726_coldrampskeleton4/RAnalysis/") #set working

##### Ramp Experiment Data #####
path.p <- ("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/20170726_coldrampskeleton4/RAnalysis/") #the location of all your respirometry files for particular ramp experiment
file.names<-list.files(path = paste0(path.p,"/Photo"), pattern = "txt$") #list all csv/txt file names in the folder. These files have already been corrected using python script
Photo.R <- data.frame(matrix(NA, nrow=length(file.names)*2, ncol=3)) #generate a 3 column dataframe with specific column names
colnames(Photo.R) <- c("Fragment.ID","Intercept", "umol.L.sec")

for(i in 1:length(file.names)) { # for every file in list start at the first and run this following function
  Photo.Data1 <-read.table(file.path(paste0(path.p,"/Photo"),file.names[i]), header=T, sep="\t", na.string="NA", as.is=TRUE, fileEncoding="latin1") #reads in the data files
  Photo.Data1  <- Photo.Data1[,c(2,17,7)] #subset columns of interest (time, O2, temp)
  Photo.Data1$Time.hh.mm.ss <- as.POSIXct(Photo.Data1$Time.hh.mm.ss,format="%H:%M:%S", tz = "") #convert time from character to time
#  brk <- which(diff(Photo.Data1$Time) > 30) #look for breaks in time of 30 seconds or more
  Photo <- Photo.Data1
  respname<- paste0(sub("P.*", "", file.names[i]),"R",sub(".*P", "", file.names[i]))
  Resp.Data1 <-read.table(file.path(paste0(path.p,"/Resp"),respname), header=T, sep="\t", na.string="NA", as.is=TRUE, fileEncoding="latin1") #reads in the data files
  Resp.Data1  <- Resp.Data1[,c(2,17,7)] #subset columns of interest (time, O2, temp)
  Resp.Data1$Time.hh.mm.ss <- as.POSIXct(Resp.Data1$Time.hh.mm.ss,format="%H:%M:%S", tz = "")   
  Resp <- Resp.Data1 #subset by break in time stamp keeping everything before break
  lt.levs <- list(Photo, Resp) #list levels of segmentation

for(j in 1:length(lt.levs)){    #go through Photo then Resp for each pair
  if(j==1){
  	z="P"
  }
  else{
  	z="R"
  }
  Photo.Data <- as.data.frame(lt.levs[j])
  n<-dim(Photo.Data)[1] #identify length of data
  Photo.Data <-Photo.Data [8:n,] #start at data point ~2 minute in to avoid excess noise from start of run
  n<-dim(Photo.Data )[1] #list length of trimmed data
  Photo.Data $sec <- seq(from=1,to=n*15, by=15) #set seconds by one from start to finish of run
  
  #Save plot prior to and after data thinning to make sure thinning is not too extreme
  rename <- paste0(sub("P.*", "", file.names[i]),strsplit(sub(".*P-", "", file.names[i]),"_",)[[1]][1])
  pdf(paste0("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/20170726_coldrampskeleton4/RAnalysis/Output/",rename,"_",z,"_thinning.pdf"))
  par(omi=rep(0.3, 4)) #set size of the outer margins in inches
  par(mfrow=c(1,2)) #set number of rows and columns in multi plot graphic
  plot(SalCorcO2..µmol.L. ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE) #plot data as a function of time
  usr  <-  par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
  whiteGrid()
  box()
  points(Photo.Data $SalCorcO2..µmol.L. ~ Photo.Data $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
  axis(1)
  axis(2, las=1)

###thins the data before regressing  
  Photo.Data   <-  thinData(Photo.Data , by=2)$newData1 #thin data by every 2 points/30sec
  Photo.Data $sec <- seq(from=1,to=nrow(Photo.Data)*30, by=30) #maintain numeric values for time
  plot(SalCorcO2..µmol.L. ~ sec, data=Photo.Data , xlab='Time (seconds)', ylab=substitute(' O'[2]~' (µmol/L)'), type='n', axes=FALSE) #plot thinned data
  usr  <-  par('usr')
  rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
  whiteGrid()
  box()
  points(Photo.Data $SalCorcO2..µmol.L. ~ Photo.Data $sec, pch=16, col=transparentColor('dodgerblue2', 0.6), cex=1.1)
  axis(1)
  axis(2, las=1)
  dev.off()
  
  Regs  <-  rankLocReg(xall=Photo.Data $sec, yall=Photo.Data $SalCorcO2..µmol.L., alpha=0.3, 
                       method="pc", verbose=TRUE) 
                       #alpha of 0.3 means that with ~80 data points per run a minimum of 6 minutes of data is included in the linear regression
  pdf(paste0("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/20170726_coldrampskeleton4/RAnalysis/Output/",rename,"_",z,"_regression_thin.pdf"))
  plot(Regs)
  dev.off()

  pdf(paste0("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/20170726_coldrampskeleton4/RAnalysis/Output/",rename,"_",z,"_regressionranks_thin.pdf"))
  outputRankLocRegPlot(Regs)
  dev.off()  

  s <- seq(0,nrow(Photo.R),length(lt.levs)) #to order the file output sequence in correct order in data frame
  Photo.R[j+s[i],2:3] <- Regs$allRegs[1,c(4,5)] #inserts slope and intercept in the dataframe
  Photo.R[j+s[i],1] <- paste0(rename,"_",z) #stores the file name in the Date column
}
}

Photo.R       #Photo.R is a data frame with 3 columns (Fragment.ID, Intercept and metabolic rate in umol.L.sec) and 90 rows (P and R in each chamber at each temperature)
# Format of Fragment.ID is: 20170707_18C_ch10_P or _R


#####Calculated raw values######

#####Correcting for blank/drift and surface area####
PHO <- Photo.R[seq(1, nrow(Photo.R), 2), ] #sets every other line in Photo.R as P value 
RES <- Photo.R[seq(2, nrow(Photo.R), 2), ] #sets every other line in Photo.R as R

PHO$Fragment.ID <- sub(".*P-", "", PHO$Fragment.ID)
RES$Fragment.ID <- sub(".*R-", "", RES$Fragment.ID)

#Load Fragment Sample Info
Sample.Info <- read.csv(file="Nubbin_Sample_Info_20170726.csv", header=T) #read in volume and sample.info data, including symbiotic state, origin, and temp of measurement
Sample.Info$Vol.L <- Sample.Info$Chamber.Vol.mL/1000 #calculate volume

#Merge P and R data with sample info data by the Fragment.ID column
Resp <- merge(RES,Sample.Info, by="Fragment.ID")
Photo <- merge(PHO,Sample.Info, by="Fragment.ID")

#Account for chamber volume  to convert R and P rates to units of umol s-1
Resp$umol.sec <- Resp$umol.L.sec*Resp$Vol.L
Photo$umol.sec <- Photo$umol.L.sec*Photo$Vol.L

#Account for photo blank rate (chamber with no coral), separate by temperature
photo.blnk <- aggregate(umol.sec ~ Temp+Genotype, data=Photo, mean) #changed from umol.L.sec to umol.sec...is this right??
#photo.blnk is data frame with 3 columns, Temp Genotype and umol.sec
photo.Blanks <- subset(photo.blnk, Genotype == "Blank") #subsets blank rate from each temp
photoBlank6 <-photo.Blanks[1,3]
photoBlank9 <-photo.Blanks[2,3]
photoBlank12 <-photo.Blanks[3,3]
photoBlank15 <-photo.Blanks[4,3]
photoBlank18 <-photo.Blanks[5,3]

#Account for respiration blank rate, separate by temperature
resp.blnk <- aggregate(umol.sec ~ Temp+Genotype, data=Resp, mean)
resp.Blanks <- subset(resp.blnk, Genotype == "Blank")
respBlank6 <-resp.Blanks[1,3]
respBlank9 <-resp.Blanks[2,3]
respBlank12 <-resp.Blanks[3,3]
respBlank15 <-resp.Blanks[4,3]
respBlank18 <-resp.Blanks[5,3]

#Subtract Blank, separately for each temperature and then re-combine into one data frame (Photo and Resp)
#Photosynthesis
Photo.6 <- subset(Photo, Temp == "6")
Photo.9 <- subset(Photo, Temp == "9")
Photo.12 <- subset(Photo, Temp == "12")
Photo.15 <- subset(Photo, Temp == "15")
Photo.18 <- subset(Photo, Temp == "18")

Photo.6$umol.sec.corr <- Photo.6$umol.sec-photoBlank6
Photo.9$umol.sec.corr <- Photo.9$umol.sec-photoBlank9
Photo.12$umol.sec.corr <- Photo.12$umol.sec-photoBlank12
Photo.15$umol.sec.corr <- Photo.15$umol.sec-photoBlank15
Photo.18$umol.sec.corr <- Photo.18$umol.sec-photoBlank18

Photo <- rbind(Photo.6,Photo.9,Photo.12,Photo.15,Photo.18)

#Respiration
Resp.6 <- subset(Resp, Temp == "6")
Resp.9 <- subset(Resp, Temp == "9")
Resp.12 <- subset(Resp, Temp == "12")
Resp.15 <- subset(Resp, Temp == "15")
Resp.18 <- subset(Resp, Temp == "18")

Resp.6$umol.sec.corr <- Resp.6$umol.sec-respBlank6
Resp.9$umol.sec.corr <- Resp.9$umol.sec-respBlank9
Resp.12$umol.sec.corr <- Resp.12$umol.sec-respBlank12
Resp.15$umol.sec.corr <- Resp.15$umol.sec-respBlank15
Resp.18$umol.sec.corr <- Resp.18$umol.sec-respBlank18

Resp <- rbind(Resp.6,Resp.9,Resp.12,Resp.15,Resp.18)

#normalize to surface area and h-1
Resp$umol.cm2.hr <- (Resp$umol.sec.corr*3600)/Resp$Surf.Area.cm2
Photo$umol.cm2.hr <- (Photo$umol.sec.corr*3600)/Photo$Surf.Area.cm2

####### Now you have Corrected Photosynthesis and Respiration Results #######

#remove blanks
Photo <- subset(Photo, Genotype!= "Blank")
Resp <- subset(Resp, Genotype!= "Blank")

#Write out photosynthesis and respiration data to separate .csv files
write.csv(Photo, row.names=FALSE, file="~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/20170726_coldrampskeleton4/RAnalysis/Output/Photo_Resp_Output/Photosynthesis.rates.csv")
write.csv(Resp, row.names=FALSE, file="~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/20170726_coldrampskeleton4/RAnalysis/Output/Photo_Resp_Output/Respiration.rates.csv")

#######

#Analysis of all holobiont and skeleton data by type of ramp (hot/cold)
#Only use this section of code if you have already combined all individual ramp days into one master data sheet (copy and paste all Photo and Resp .csv files written out above into one master)
setwd("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData")

Photo.holobiont.heat <- read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Photosynthesis_master_holobiont_heat.csv")
Resp.holobiont.heat <- read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respiration_master_holobiont_heat.csv")

Photo.holobiont.cold <- read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Photosynthesis_master_holobiont_cold.csv")
Resp.holobiont.cold <- read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respiration_master_holobiont_cold.csv")

Photo.skeleton.heat <- read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Photosynthesis_master_skeleton_heat.csv")
Resp.skeleton.heat <- read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respiration_master_skeleton_heat.csv")

Photo.skeleton.cold <- read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Photosynthesis_master_skeleton_cold.csv")
Resp.skeleton.cold <- read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respiration_master_skeleton_cold.csv")

Resp.holobiont.all <- read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respiration_master_holobiont.csv")

Photo.holobiont.all <- read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Photosynthesis_master_holobiont.csv")
###########

#bind together corrected data and important frag into, then calculate gross photosynthesis 
#Pnet -- Rdark
#Keeping heat and cold ramp data separate, and skeleton and holobiont data separate within temp condition

##########Combined Heat and Cold Ramps
head(Resp.holobiont.all)
head(Photo.holobiont.all)

resp.data.holobiont <- cbind(Photo.holobiont.all[,c(4,5,6,10,11,12,16)],Resp.holobiont.all[,c(16)]) #double check column numbers here!!!
colnames(resp.data.holobiont)[7] <- "Pnet_umol.cm2.hr_holo"
colnames(resp.data.holobiont) [8] <- "Rdark_umol.cm2.hr_holo"

resp.data.holobiont <-cbind(resp.data.holobiont[0:5],"colorstate"=paste(resp.data.holobiont$Origin,resp.data.holobiont$Color, sep="_"),resp.data.holobiont[6:ncol(resp.data.holobiont)])

resp.data.holobiont <-cbind(resp.data.holobiont[0:3],"GxTemp"=paste(resp.data.holobiont$Genotype,resp.data.holobiont$Temp, sep="_"),resp.data.holobiont[4:ncol(resp.data.holobiont)])

##########Heat Ramps
### Holobiont Heat Data
head(Photo.holobiont.heat)
head(Resp.holobiont.heat)

resp.data.heat.holobiont <- cbind(Photo.holobiont.heat[,c(4,5,6,10,11,12,16)],Resp.holobiont.heat[,c(16)]) #double check column numbers here!!!
colnames(resp.data.heat.holobiont)[7] <- "Pnet_umol.cm2.hr_holo"
colnames(resp.data.heat.holobiont) [8] <- "Rdark_umol.cm2.hr_holo"

resp.data.heat.holobiont <-cbind(resp.data.heat.holobiont[0:5],"colorstate"=paste(resp.data.heat.holobiont$Origin_holo,resp.data.heat.holobiont$Color_holo, sep="_"),resp.data.heat.holobiont[6:ncol(resp.data.heat.holobiont)])

resp.data.heat.holobiont <-cbind(resp.data.heat.holobiont[0:3],"GxTemp"=paste(resp.data.heat.holobiont$Genotype_holo,resp.data.heat.holobiont$Temp_holo, sep="_"),resp.data.heat.holobiont[4:ncol(resp.data.heat.holobiont)])

#calculate Pgross
resp.data.heat.holobiont$Pgross_umol.cm2.hr_holo <- resp.data.heat.holobiont$Pnet_umol.cm2.hr_holo-resp.data.heat.holobiont$Rdark_umol.cm2.hr_holo #resp.data.heat.holobiont is 160 lines long here

### Skeleton Heat Data

head(Photo.skeleton.heat)
head(Resp.skeleton.heat)

resp.data.heat.skeleton <- cbind(Photo.skeleton.heat[,c(4,5,6,10,11,12,16)],Resp.skeleton.heat[,c(16)]) #double check column numbers here!!!
colnames(resp.data.heat.skeleton)[7] <- "Pnet_umol.cm2.hr_skel"
colnames(resp.data.heat.skeleton) [8] <- "Rdark_umol.cm2.hr_skel"

resp.data.heat.skeleton <-cbind(resp.data.heat.skeleton[0:5],"colorstate"=paste(resp.data.heat.skeleton$Origin_skel,resp.data.heat.skeleton$Color_skel, sep="_"),resp.data.heat.skeleton[6:ncol(resp.data.heat.skeleton)])

resp.data.heat.skeleton <-cbind(resp.data.heat.skeleton[0:3],"GxTemp"=paste(resp.data.heat.skeleton$Genotype_skel, resp.data.heat.skeleton$Temp_skel, sep="_"), resp.data.heat.skeleton[4:ncol(resp.data.heat.skeleton)])

#calculate Pgross
resp.data.heat.skeleton$Pgross_umol.cm2.hr_skel <- resp.data.heat.skeleton$Pnet_umol.cm2.hr_skel-resp.data.heat.skeleton $Rdark_umol.cm2.hr_skel #resp.data.heat.skeleton is 160 lines long here

Heat.Data <- merge(resp.data.heat.skeleton,resp.data.heat.holobiont, by="GxTemp") #155 lines long, lose the samples that don't have matching genotypes

Heat.Data$Pnet_umol.cm2.hr_corr <- Heat.Data$Pnet_umol.cm2.hr_holo-Heat.Data$Pnet_umol.cm2.hr_skel
Heat.Data$Rdark_umol.cm2.hr_corr <- Heat.Data$Rdark_umol.cm2.hr_holo-Heat.Data$Rdark_umol.cm2.hr_skel
Heat.Data$Pgross_umol.cm2.hr_corr <- Heat.Data$Pgross_umol.cm2.hr_holo-Heat.Data$Pgross_umol.cm2.hr_skel

write.csv(Heat.Data, file="~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Heat_Data.csv")

##########Cold Ramps
### Holobiont Cold Data
head(Photo.holobiont.cold)
head(Resp.holobiont.cold)

resp.data.cold.holobiont <- cbind(Photo.holobiont.cold[,c(4,5,6,10,11,12,16)],Resp.holobiont.cold[,c(16)]) #double check column numbers here!!!
colnames(resp.data.cold.holobiont)[7] <- "Pnet_umol.cm2.hr_holo"
colnames(resp.data.cold.holobiont) [8] <- "Rdark_umol.cm2.hr_holo"

resp.data.cold.holobiont <-cbind(resp.data.cold.holobiont[0:5],"colorstate"=paste(resp.data.cold.holobiont$Origin_holo, resp.data.cold.holobiont$Color_holo, sep="_"), resp.data.cold.holobiont[6:ncol(resp.data.cold.holobiont)])

resp.data.cold.holobiont <-cbind(resp.data.cold.holobiont[0:3],"GxTemp"=paste(resp.data.cold.holobiont$Genotype_holo, resp.data.cold.holobiont $Temp_holo, sep="_"), resp.data.cold.holobiont[4:ncol(resp.data.cold.holobiont)])

#calculate Pgross
resp.data.cold.holobiont$Pgross_umol.cm2.hr_holo <- resp.data.cold.holobiont$Pnet_umol.cm2.hr_holo-resp.data.cold.holobiont$Rdark_umol.cm2.hr_holo

### Skeleton Cold Data
head(Photo.skeleton.cold)
head(Resp.skeleton.cold)

resp.data.cold.skeleton <- cbind(Photo.skeleton.cold[,c(4,5,6,10,11,12,16)],Resp.skeleton.cold[,c(16)]) #double check column numbers here!!!
colnames(resp.data.cold.skeleton)[7] <- "Pnet_umol.cm2.hr_skel"
colnames(resp.data.cold.skeleton) [8] <- "Rdark_umol.cm2.hr_skel"

resp.data.cold.skeleton <-cbind(resp.data.cold.skeleton[0:5],"colorstate"=paste(resp.data.cold.skeleton$Origin_skel, resp.data.cold.skeleton $Color_skel, sep="_"), resp.data.cold.skeleton[6:ncol(resp.data.cold.skeleton)])

resp.data.cold.skeleton <-cbind(resp.data.cold.skeleton[0:3],"GxTemp"=paste(resp.data.cold.skeleton$Genotype_skel, resp.data.cold.skeleton $Temp_skel, sep="_"), resp.data.cold.skeleton[4:ncol(resp.data.cold.skeleton)])

#calculate Pgross
resp.data.cold.skeleton$Pgross_umol.cm2.hr_skel <- resp.data.cold.skeleton$Pnet_umol.cm2.hr_skel-resp.data.cold.skeleton $Rdark_umol.cm2.hr_skel


Cold.Data <- merge(resp.data.cold.skeleton,resp.data.cold.holobiont, by="GxTemp")

#calculate holobiont - skeleton by genotype calculations (Pn, Pg, Rdark), gives coral animal rates only
Cold.Data$Pnet_umol.cm2.hr_corr <- Cold.Data$Pnet_umol.cm2.hr_holo-Cold.Data$Pnet_umol.cm2.hr_skel 
Cold.Data$Rdark_umol.cm2.hr_corr <- Cold.Data$Rdark_umol.cm2.hr_holo-Cold.Data$Rdark_umol.cm2.hr_skel
Cold.Data$Pgross_umol.cm2.hr_corr <- Cold.Data$Pgross_umol.cm2.hr_holo-Cold.Data$Pgross_umol.cm2.hr_skel

write.csv(Cold.Data, file="~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Cold_Data.csv")

########################################################################################

########################################################################################

#Calculate mean, se, sd, n and plot by conditions

###########Combined Heat and Cold data
Rdark.mean_holo <- aggregate(Rdark_umol.cm2.hr_holo ~ colorstate*Temp, data=resp.data.holobiont, mean)
Rdark.se_holo <- aggregate(Rdark_umol.cm2.hr_holo ~ colorstate*Temp, data= resp.data.holobiont, std.error)
Rdark.sd_holo <- aggregate(Rdark_umol.cm2.hr_holo ~ colorstate*Temp, data= resp.data.holobiont, sd)
Rdark.n_holo <- aggregate(Rdark_umol.cm2.hr_holo ~ colorstate*Temp, data= resp.data.holobiont, length)

Rdark.means_holo <- cbind(Rdark.mean_holo, Rdark.se_holo$Rdark_umol.cm2.hr_holo, Rdark.sd_holo$Rdark_umol.cm2.hr_holo, Rdark.n_holo$Rdark_umol.cm2.hr_holo) #combine mean and standard error results
colnames(Rdark.means_holo) <- c("ColorState", "Temp", "mean", "se", "sd", "n")  #rename columns to describe contents
Rdark.means_holo$Metric <- "Holobiont Dark Respiration"
Rdark.means_holo$abs_mean <- abs(Rdark.means_holo$mean)

########### Heat Data
###Net Photosynthesis
#Holobiont
Pnet.mean_holo <- aggregate(Pnet_umol.cm2.hr_holo ~ colorstate.y*Temp_holo, data=Heat.Data, mean)
Pnet.se_holo <- aggregate(Pnet_umol.cm2.hr_holo ~ colorstate.y*Temp_holo, data= Heat.Data, std.error)
Pnet.sd_holo <- aggregate(Pnet_umol.cm2.hr_holo ~ colorstate.y*Temp_holo, data= Heat.Data, sd)
Pnet.n_holo <- aggregate(Pnet_umol.cm2.hr_holo ~ colorstate.y*Temp_holo, data= Heat.Data, length)

Pnet.means_holo <- cbind(Pnet.mean_holo, Pnet.se_holo$Pnet_umol.cm2.hr_holo, Pnet.sd_holo$Pnet_umol.cm2.hr_holo, Pnet.n_holo$Pnet_umol.cm2.hr_holo) #combine mean and standard error results
colnames(Pnet.means_holo) <- c("ColorState", "Temp", "mean", "se", "sd", "n")  #rename columns to describe contents
Pnet.means_holo$Metric <- "Holobiont Net Photosynthesis"

#Skeleton
Pnet.mean_skel <- aggregate(Pnet_umol.cm2.hr_skel ~ colorstate.x*Temp_skel, data=Heat.Data, mean)
Pnet.se_skel <- aggregate(Pnet_umol.cm2.hr_skel ~ colorstate.x*Temp_skel, data= Heat.Data, std.error)
Pnet.sd_skel <- aggregate(Pnet_umol.cm2.hr_skel ~ colorstate.x*Temp_skel, data= Heat.Data, sd)
Pnet.n_skel <- aggregate(Pnet_umol.cm2.hr_skel ~ colorstate.x*Temp_skel, data= Heat.Data, length)

Pnet.means_skel <- cbind(Pnet.mean_skel, Pnet.se_skel$Pnet_umol.cm2.hr_skel, Pnet.sd_skel$Pnet_umol.cm2.hr_skel, Pnet.n_skel$Pnet_umol.cm2.hr_skel) #combine mean and standard error results
colnames(Pnet.means_skel) <- c("ColorState", "Temp", "mean", "se", "sd", "n")  #rename columns to describe contents
Pnet.means_skel$Metric <- "Skeleton Net Photosynthesis"

#Corrected (animal only)
Pnet.mean_corr <- aggregate(Pnet_umol.cm2.hr_corr ~ colorstate.x*Temp_skel, data=Heat.Data, mean)
Pnet.se_corr <- aggregate(Pnet_umol.cm2.hr_corr ~ colorstate.x*Temp_skel, data= Heat.Data, std.error)
Pnet.sd_corr <- aggregate(Pnet_umol.cm2.hr_corr ~ colorstate.x*Temp_skel, data= Heat.Data, sd)
Pnet.n_corr <- aggregate(Pnet_umol.cm2.hr_corr ~ colorstate.x*Temp_skel, data= Heat.Data, length)

Pnet.means_corr <- cbind(Pnet.mean_corr, Pnet.se_corr$Pnet_umol.cm2.hr_corr, Pnet.sd_corr$Pnet_umol.cm2.hr_corr, Pnet.n_corr$Pnet_umol.cm2.hr_corr) #combine mean and standard error results
colnames(Pnet.means_corr) <- c("ColorState", "Temp", "mean", "se", "sd", "n")  #rename columns to describe contents
Pnet.means_corr$Metric <- "Corrected Net Photosynthesis"

###Dark Respiration
#Holobiont
Rdark.mean_holo <- aggregate(Rdark_umol.cm2.hr_holo ~ colorstate.y*Temp_holo, data=Heat.Data, mean)
Rdark.se_holo <- aggregate(Rdark_umol.cm2.hr_holo ~ colorstate.y*Temp_holo, data= Heat.Data, std.error)
Rdark.sd_holo <- aggregate(Rdark_umol.cm2.hr_holo ~ colorstate.y*Temp_holo, data= Heat.Data, sd)
Rdark.n_holo <- aggregate(Rdark_umol.cm2.hr_holo ~ colorstate.y*Temp_holo, data= Heat.Data, length)

Rdark.means_holo <- cbind(Rdark.mean_holo, Rdark.se_holo$Rdark_umol.cm2.hr_holo, Rdark.sd_holo$Rdark_umol.cm2.hr_holo, Rdark.n_holo$Rdark_umol.cm2.hr_holo) #combine mean and standard error results
colnames(Rdark.means_holo) <- c("ColorState", "Temp", "mean", "se", "sd", "n")  #rename columns to describe contents
Rdark.means_holo$Metric <- "Holobiont Dark Respiration"
Rdark.means_holo$abs_mean <- abs(Rdark.means_holo$mean)

#Skeleton
Rdark.mean_skel <- aggregate(Rdark_umol.cm2.hr_skel ~ colorstate.x*Temp_skel, data=Heat.Data, mean)
Rdark.se_skel <- aggregate(Rdark_umol.cm2.hr_skel ~ colorstate.x*Temp_skel, data= Heat.Data, std.error)
Rdark.sd_skel <- aggregate(Rdark_umol.cm2.hr_skel ~ colorstate.x*Temp_skel, data= Heat.Data, sd)
Rdark.n_skel <- aggregate(Rdark_umol.cm2.hr_skel ~ colorstate.x*Temp_skel, data= Heat.Data, length)

Rdark.means_skel <- cbind(Rdark.mean_skel, Rdark.se_skel$Rdark_umol.cm2.hr_skel, Rdark.sd_skel$Rdark_umol.cm2.hr_skel, Rdark.n_skel$Rdark_umol.cm2.hr_skel) #combine mean and standard error results
colnames(Rdark.means_skel) <- c("ColorState", "Temp", "mean", "se", "sd", "n")  #rename columns to describe contents
Rdark.means_skel$Metric <- "Skeleton Dark Respiration"
Rdark.means_skel$abs_mean <- abs(Rdark.means_skel$mean)

#Corrected (animal only)
Rdark.mean_corr <- aggregate(Rdark_umol.cm2.hr_corr ~ colorstate.x*Temp_skel, data=Heat.Data, mean)
Rdark.se_corr <- aggregate(Rdark_umol.cm2.hr_corr ~ colorstate.x*Temp_skel, data= Heat.Data, std.error)
Rdark.sd_corr <- aggregate(Rdark_umol.cm2.hr_corr ~ colorstate.x*Temp_skel, data= Heat.Data, sd)
Rdark.n_corr <- aggregate(Rdark_umol.cm2.hr_corr ~ colorstate.x*Temp_skel, data= Heat.Data, length)

Rdark.means_corr <- cbind(Rdark.mean_corr, Rdark.se_corr$Rdark_umol.cm2.hr_corr, Rdark.sd_corr$Rdark_umol.cm2.hr_corr, Rdark.n_corr$Rdark_umol.cm2.hr_corr) #combine mean and standard error results
colnames(Rdark.means_corr) <- c("ColorState", "Temp", "mean", "se", "sd", "n")  #rename columns to describe contents
Rdark.means_corr$Metric <- "Corrected Dark Respiration"
Rdark.means_corr$abs_mean <- abs(Rdark.means_corr$mean)

###Gross Photosynthesis
#Holobiont
Pgross.mean_holo <- aggregate(Pgross_umol.cm2.hr_holo ~ colorstate.y*Temp_holo, data=Heat.Data, mean)
Pgross.se_holo <- aggregate(Pgross_umol.cm2.hr_holo ~ colorstate.y*Temp_holo, data= Heat.Data, std.error)
Pgross.sd_holo <- aggregate(Pgross_umol.cm2.hr_holo ~ colorstate.y*Temp_holo, data= Heat.Data, sd)
Pgross.n_holo <- aggregate(Pgross_umol.cm2.hr_holo ~ colorstate.y*Temp_holo, data= Heat.Data, length)

Pgross.means_holo <- cbind(Pgross.mean_holo, Pgross.se_holo$Pgross_umol.cm2.hr_holo, Pgross.sd_holo$Pgross_umol.cm2.hr_holo, Pgross.n_holo$Pgross_umol.cm2.hr_holo) #combine mean and standard error results
colnames(Pgross.means_holo) <- c("ColorState", "Temp", "mean", "se", "sd", "n")  #rename columns to describe contents
Pgross.means_holo$Metric <- "Holobiont Gross Photosynthesis"

#Skeleton
Pgross.mean_skel <- aggregate(Pgross_umol.cm2.hr_skel ~ colorstate.x*Temp_skel, data=Heat.Data, mean)
Pgross.se_skel <- aggregate(Pgross_umol.cm2.hr_skel ~ colorstate.x*Temp_skel, data= Heat.Data, std.error)
Pgross.sd_skel <- aggregate(Pgross_umol.cm2.hr_skel ~ colorstate.x*Temp_skel, data= Heat.Data, sd)
Pgross.n_skel <- aggregate(Pgross_umol.cm2.hr_skel ~ colorstate.x*Temp_skel, data= Heat.Data, length)

Pgross.means_skel <- cbind(Pgross.mean_skel, Pgross.se_skel$Pgross_umol.cm2.hr_skel, Pgross.sd_skel$Pgross_umol.cm2.hr_skel, Pgross.n_skel$Pgross_umol.cm2.hr_skel) #combine mean and standard error results
colnames(Pgross.means_skel) <- c("ColorState", "Temp", "mean", "se", "sd", "n")  #rename columns to describe contents
Pgross.means_skel$Metric <- "Skeleton Gross Photosynthesis"

#Corrected (animal only)
Pgross.mean_corr <- aggregate(Pgross_umol.cm2.hr_corr ~ colorstate.x*Temp_skel, data=Heat.Data, mean)
Pgross.se_corr <- aggregate(Pgross_umol.cm2.hr_corr ~ colorstate.x*Temp_skel, data= Heat.Data, std.error)
Pgross.sd_corr <- aggregate(Pgross_umol.cm2.hr_corr ~ colorstate.x*Temp_skel, data= Heat.Data, sd)
Pgross.n_corr <- aggregate(Pgross_umol.cm2.hr_corr ~ colorstate.x*Temp_skel, data= Heat.Data, length)

Pgross.means_corr <- cbind(Pgross.mean_corr, Pgross.se_corr$Pgross_umol.cm2.hr_corr, Pgross.sd_corr$Pgross_umol.cm2.hr_corr, Pgross.n_corr$Pgross_umol.cm2.hr_corr) #combine mean and standard error results
colnames(Pgross.means_corr) <- c("ColorState", "Temp", "mean", "se", "sd", "n")  #rename columns to describe contents
Pgross.means_corr$Metric <- "Corrected Gross Photosynthesis"

#####################Cold Data
###Net Photosynthesis
#Holobiont
Pnet.mean_holo <- aggregate(Pnet_umol.cm2.hr_holo ~ colorstate.y*Temp_holo, data=Cold.Data, mean)
Pnet.se_holo <- aggregate(Pnet_umol.cm2.hr_holo ~ colorstate.y*Temp_holo, data= Cold.Data, std.error)
Pnet.sd_holo <- aggregate(Pnet_umol.cm2.hr_holo ~ colorstate.y*Temp_holo, data= Cold.Data, sd)
Pnet.n_holo <- aggregate(Pnet_umol.cm2.hr_holo ~ colorstate.y*Temp_holo, data= Cold.Data, length)

Pnet.means_holo <- cbind(Pnet.mean_holo, Pnet.se_holo$Pnet_umol.cm2.hr_holo, Pnet.sd_holo$Pnet_umol.cm2.hr_holo, Pnet.n_holo$Pnet_umol.cm2.hr_holo) #combine mean and standard error results
colnames(Pnet.means_holo) <- c("ColorState", "Temp", "mean", "se", "sd", "n")  #rename columns to describe contents
Pnet.means_holo$Metric <- "Holobiont Net Photosynthesis"

#Skeleton
Pnet.mean_skel <- aggregate(Pnet_umol.cm2.hr_skel ~ colorstate.x*Temp_skel, data= Cold.Data, mean)
Pnet.se_skel <- aggregate(Pnet_umol.cm2.hr_skel ~ colorstate.x*Temp_skel, data= Cold.Data, std.error)
Pnet.sd_skel <- aggregate(Pnet_umol.cm2.hr_skel ~ colorstate.x*Temp_skel, data= Cold.Data, sd)
Pnet.n_skel <- aggregate(Pnet_umol.cm2.hr_skel ~ colorstate.x*Temp_skel, data= Cold.Data, length)

Pnet.means_skel <- cbind(Pnet.mean_skel, Pnet.se_skel$Pnet_umol.cm2.hr_skel, Pnet.sd_skel$Pnet_umol.cm2.hr_skel, Pnet.n_skel$Pnet_umol.cm2.hr_skel) #combine mean and standard error results
colnames(Pnet.means_skel) <- c("ColorState", "Temp", "mean", "se", "sd", "n")  #rename columns to describe contents
Pnet.means_skel$Metric <- "Skeleton Net Photosynthesis"

#Corrected (animal only)
Pnet.mean_corr <- aggregate(Pnet_umol.cm2.hr_corr ~ colorstate.x*Temp_skel, data= Cold.Data, mean)
Pnet.se_corr <- aggregate(Pnet_umol.cm2.hr_corr ~ colorstate.x*Temp_skel, data= Cold.Data, std.error)
Pnet.sd_corr <- aggregate(Pnet_umol.cm2.hr_corr ~ colorstate.x*Temp_skel, data= Cold.Data, sd)
Pnet.n_corr <- aggregate(Pnet_umol.cm2.hr_corr ~ colorstate.x*Temp_skel, data= Cold.Data, length)

Pnet.means_corr <- cbind(Pnet.mean_corr, Pnet.se_corr$Pnet_umol.cm2.hr_corr, Pnet.sd_corr$Pnet_umol.cm2.hr_corr, Pnet.n_corr$Pnet_umol.cm2.hr_corr) #combine mean and standard error results
colnames(Pnet.means_corr) <- c("ColorState", "Temp", "mean", "se", "sd", "n")  #rename columns to describe contents
Pnet.means_corr$Metric <- "Corrected Net Photosynthesis"

###Dark Respiration
#Holobiont
Rdark.mean_holo <- aggregate(Rdark_umol.cm2.hr_holo ~ colorstate.y*Temp_holo, data= Cold.Data, mean)
Rdark.se_holo <- aggregate(Rdark_umol.cm2.hr_holo ~ colorstate.y*Temp_holo, data= Cold.Data, std.error)
Rdark.sd_holo <- aggregate(Rdark_umol.cm2.hr_holo ~ colorstate.y*Temp_holo, data= Cold.Data, sd)
Rdark.n_holo <- aggregate(Rdark_umol.cm2.hr_holo ~ colorstate.y*Temp_holo, data= Cold.Data, length)

Rdark.means_holo <- cbind(Rdark.mean_holo, Rdark.se_holo$Rdark_umol.cm2.hr_holo, Rdark.sd_holo$Rdark_umol.cm2.hr_holo, Rdark.n_holo$Rdark_umol.cm2.hr_holo) #combine mean and standard error results
colnames(Rdark.means_holo) <- c("ColorState", "Temp", "mean", "se", "sd", "n")  #rename columns to describe contents
Rdark.means_holo$Metric <- "Holobiont Dark Respiration"
Rdark.means_holo$abs_mean <- abs(Rdark.means_holo$mean)

#Skeleton
Rdark.mean_skel <- aggregate(Rdark_umol.cm2.hr_skel ~ colorstate.x*Temp_skel, data= Cold.Data, mean)
Rdark.se_skel <- aggregate(Rdark_umol.cm2.hr_skel ~ colorstate.x*Temp_skel, data= Cold.Data, std.error)
Rdark.sd_skel <- aggregate(Rdark_umol.cm2.hr_skel ~ colorstate.x*Temp_skel, data= Cold.Data, sd)
Rdark.n_skel <- aggregate(Rdark_umol.cm2.hr_skel ~ colorstate.x*Temp_skel, data= Cold.Data, length)

Rdark.means_skel <- cbind(Rdark.mean_skel, Rdark.se_skel$Rdark_umol.cm2.hr_skel, Rdark.sd_skel$Rdark_umol.cm2.hr_skel, Rdark.n_skel$Rdark_umol.cm2.hr_skel) #combine mean and standard error results
colnames(Rdark.means_skel) <- c("ColorState", "Temp", "mean", "se", "sd", "n")  #rename columns to describe contents
Rdark.means_skel$Metric <- "Skeleton Dark Respiration"
Rdark.means_skel$abs_mean <- abs(Rdark.means_skel$mean)

#Corrected (animal only)
Rdark.mean_corr <- aggregate(Rdark_umol.cm2.hr_corr ~ colorstate.x*Temp_skel, data= Cold.Data, mean)
Rdark.se_corr <- aggregate(Rdark_umol.cm2.hr_corr ~ colorstate.x*Temp_skel, data= Cold.Data, std.error)
Rdark.sd_corr <- aggregate(Rdark_umol.cm2.hr_corr ~ colorstate.x*Temp_skel, data= Cold.Data, sd)
Rdark.n_corr <- aggregate(Rdark_umol.cm2.hr_corr ~ colorstate.x*Temp_skel, data= Cold.Data, length)

Rdark.means_corr <- cbind(Rdark.mean_corr, Rdark.se_corr$Rdark_umol.cm2.hr_corr, Rdark.sd_corr$Rdark_umol.cm2.hr_corr, Rdark.n_corr$Rdark_umol.cm2.hr_corr) #combine mean and standard error results
colnames(Rdark.means_corr) <- c("ColorState", "Temp", "mean", "se", "sd", "n")  #rename columns to describe contents
Rdark.means_corr$Metric <- "Corrected Dark Respiration"
Rdark.means_corr$abs_mean <- abs(Rdark.means_corr$mean)

###Gross Photosynthesis
#Holobiont
Pgross.mean_holo <- aggregate(Pgross_umol.cm2.hr_holo ~ colorstate.y*Temp_holo, data= Cold.Data, mean)
Pgross.se_holo <- aggregate(Pgross_umol.cm2.hr_holo ~ colorstate.y*Temp_holo, data= Cold.Data, std.error)
Pgross.sd_holo <- aggregate(Pgross_umol.cm2.hr_holo ~ colorstate.y*Temp_holo, data= Cold.Data, sd)
Pgross.n_holo <- aggregate(Pgross_umol.cm2.hr_holo ~ colorstate.y*Temp_holo, data= Cold.Data, length)

Pgross.means_holo <- cbind(Pgross.mean_holo, Pgross.se_holo$Pgross_umol.cm2.hr_holo, Pgross.sd_holo$Pgross_umol.cm2.hr_holo, Pgross.n_holo$Pgross_umol.cm2.hr_holo) #combine mean and standard error results
colnames(Pgross.means_holo) <- c("ColorState", "Temp", "mean", "se", "sd", "n")  #rename columns to describe contents
Pgross.means_holo$Metric <- "Holobiont Gross Photosynthesis"

#Skeleton
Pgross.mean_skel <- aggregate(Pgross_umol.cm2.hr_skel ~ colorstate.x*Temp_skel, data= Cold.Data, mean)
Pgross.se_skel <- aggregate(Pgross_umol.cm2.hr_skel ~ colorstate.x*Temp_skel, data= Cold.Data, std.error)
Pgross.sd_skel <- aggregate(Pgross_umol.cm2.hr_skel ~ colorstate.x*Temp_skel, data= Cold.Data, sd)
Pgross.n_skel <- aggregate(Pgross_umol.cm2.hr_skel ~ colorstate.x*Temp_skel, data= Cold.Data, length)

Pgross.means_skel <- cbind(Pgross.mean_skel, Pgross.se_skel$Pgross_umol.cm2.hr_skel, Pgross.sd_skel$Pgross_umol.cm2.hr_skel, Pgross.n_skel$Pgross_umol.cm2.hr_skel) #combine mean and standard error results
colnames(Pgross.means_skel) <- c("ColorState", "Temp", "mean", "se", "sd", "n")  #rename columns to describe contents
Pgross.means_skel$Metric <- "Skeleton Gross Photosynthesis"

#Corrected
Pgross.mean_corr <- aggregate(Pgross_umol.cm2.hr_corr ~ colorstate.x*Temp_skel, data= Cold.Data, mean)
Pgross.se_corr <- aggregate(Pgross_umol.cm2.hr_corr ~ colorstate.x*Temp_skel, data= Cold.Data, std.error)
Pgross.sd_corr <- aggregate(Pgross_umol.cm2.hr_corr ~ colorstate.x*Temp_skel, data= Cold.Data, sd)
Pgross.n_corr <- aggregate(Pgross_umol.cm2.hr_corr ~ colorstate.x*Temp_skel, data= Cold.Data, length)

Pgross.means_corr <- cbind(Pgross.mean_corr, Pgross.se_corr$Pgross_umol.cm2.hr_corr, Pgross.sd_corr$Pgross_umol.cm2.hr_corr, Pgross.n_corr$Pgross_umol.cm2.hr_corr) #combine mean and standard error results
colnames(Pgross.means_corr) <- c("ColorState", "Temp", "mean", "se", "sd", "n")  #rename columns to describe contents
Pgross.means_corr$Metric <- "Corrected Gross Photosynthesis"

#####Figures#####
##Same code used to plot heat and cold data, but have to run scripts above for heat and cold data depending on what you want to plot
#Plots are done separately for holobiont, skeleton, and corrected data. All plotted together on one plot at the end of each section

#Net Photosynthesis (holobiont, skeleton, corrected)
Fig.Pn.holo <-  ggplot(Pnet.means_holo, aes(x=Temp, y=mean,  group=ColorState)) + #set up plot information
  geom_errorbar(aes(x=Temp, ymax=mean+sd, ymin=mean-sd), colour="black", width=.1, position = position_dodge(width = 0.2)) + #add standard error bars about the mean
  geom_point(aes(shape=ColorState), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(17,2,16,1)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Temperature (°C)") + #label x axis
    scale_x_continuous(breaks=c(6,9,12,15,18,22,26,29,32)) +
  ylab(bquote('Net Photosynthesis ('*mu~ 'mol' ~O[2] ~cm^-2~h^-1*')')) + #label y axis
  ylim(-1.5, 1)+
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(color = 'black', size=10), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_rect(colour='grey95'),
        legend.background = element_rect(colour='black'),
        legend.position=c(.2,.2)) + #set legend location
  ggtitle("Net Photosythesis - Holobiont") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 14, 
                                  hjust = 0)) #set title attributes
Fig.Pn.holo #view plot

Fig.Pn.skel <-  ggplot(Pnet.means_skel, aes(x=Temp, y=mean,  group=ColorState)) + #set up plot information
  geom_errorbar(aes(x=Temp, ymax=mean+sd, ymin=mean-sd), colour="black", width=.1, position = position_dodge(width = 0.2)) + #add standard error bars about the mean
  geom_point(aes(shape=ColorState), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(17,2,16,1)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Temperature (°C)") + #label x axis
    scale_x_continuous(breaks=c(6,9,12,15,18,22,26,29,32)) +
  ylab(bquote('Net Photosynthesis ('*mu~ 'mol' ~O[2] ~cm^-2~h^-1*')')) + #label y axis
  ylim(-1.5, 1)+
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(color = 'black', size=10), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_rect(colour='grey95'),
        legend.background = element_rect(colour='black'),
        legend.position="none") + #set legend location
  ggtitle("Net Photosythesis - Skeleton") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 14, 
                                  hjust = 0)) #set title attributes
Fig.Pn.skel #view plot

Fig.Pn.corr <-  ggplot(Pnet.means_corr, aes(x=Temp, y=mean,  group=ColorState)) + #set up plot information
  geom_errorbar(aes(x=Temp, ymax=mean+sd, ymin=mean-sd), colour="black", width=.1, position = position_dodge(width = 0.2)) + #add standard error bars about the mean
  geom_point(aes(shape=ColorState), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(17,2,16,1)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Temperature (°C)") + #label x axis
    scale_x_continuous(breaks=c(6,9,12,15,18,22,26,29,32)) +
  ylab(bquote('Net Photosynthesis ('*mu~ 'mol' ~O[2] ~cm^-2~h^-1*')')) + #label y axis
  ylim(-1.5, 1)+
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(color = 'black', size=10), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_rect(colour='grey95'),
        legend.background = element_rect(colour='black'),
        legend.position="none") + #set legend location
  ggtitle("Net Photosythesis - Corrected") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 14, 
                                  hjust = 0)) #set title attributes
Fig.Pn.corr #view plot

Figs.Pn.cold <- arrangeGrob(Fig.Pn.holo,Fig.Pn.skel,Fig.Pn.corr, ncol=3)
ggsave(file="~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Figs.Pn.cold.pdf", Figs.Pn.cold, width = 11, height = 6, units = c("in"), useDingbats=FALSE)

#Dark Respiration (holobiont, skeleton, corrected)
Fig.Rd.holo <-  ggplot(Rdark.means_holo, aes(x=Temp, y=abs_mean,  group=ColorState)) + #set up plot information
  geom_errorbar(aes(x=Temp, ymax=abs_mean+sd, ymin=abs_mean-sd), colour="black", width=.1, position = position_dodge(width = 0.2)) + #add standard error bars about the mean
  geom_point(aes(shape=ColorState), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(17,2,16,1)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Temperature (°C)") + #label x axis
    scale_x_continuous(breaks=c(6,9,12,15,18,22,26,29,32)) +
  ylab(bquote('Respiration ('*mu~ 'mol' ~O[2] ~cm^-2~h^-1*')')) + #label y axis
  ylim(-.5, 2.5)+
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(color = 'black', size=10), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_rect(colour='grey95'),
        legend.background = element_rect(colour='black'),
        legend.position=c(.2,.8)) + #set legend location
  ggtitle("Dark Respiration - Holobiont") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 14, 
                                  hjust = 0)) #set title attributes
Fig.Rd.holo #view plot
ggsave(file="~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Figs.Rd.all.pdf", Fig.Rd.holo, width = 11, height = 6, units = c("in"), useDingbats=FALSE)

Fig.Rd.skel <-  ggplot(Rdark.means_skel, aes(x=Temp, y=abs_mean,  group=ColorState)) + #set up plot information
  geom_errorbar(aes(x=Temp, ymax=abs_mean+sd, ymin=abs_mean-sd), colour="black", width=.1, position = position_dodge(width = 0.2)) + #add standard error bars about the mean
  geom_point(aes(shape=ColorState), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(17,2,16,1)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Temperature (°C)") + #label x axis
    scale_x_continuous(breaks=c(6,9,12,15,18,22,26,29,32)) +
  ylab(bquote('Respiration ('*mu~ 'mol' ~O[2] ~cm^-2~h^-1*')')) + #label y axis
  ylim(-.5, 1.5)+
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(color = 'black', size=10), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_rect(colour='grey95'),
        legend.background = element_rect(colour='black'),
        legend.position="none") + #set legend location
  ggtitle("Dark Respiration - Skeleton") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 14, 
                                  hjust = 0)) #set title attributes
Fig.Rd.skel #view plot

Fig.Rd.corr <-  ggplot(Rdark.means_corr, aes(x=Temp, y=abs_mean,  group=ColorState)) + #set up plot information
  geom_errorbar(aes(x=Temp, ymax=abs_mean+sd, ymin=abs_mean-sd), colour="black", width=.1, position = position_dodge(width = 0.2)) + #add standard error bars about the mean
  geom_point(aes(shape=ColorState), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(17,2,16,1)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Temperature (°C)") + #label x axis
    scale_x_continuous(breaks=c(6,9,12,15,18,22,26,29,32)) +
  ylab(bquote('Respiration ('*mu~ 'mol' ~O[2] ~cm^-2~h^-1*')')) + #label y axis
  ylim(-.5, 1.5)+
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(color = 'black', size=10), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_rect(colour='grey95'),
        legend.background = element_rect(colour='black'),
        legend.position="none") + #set legend location
  ggtitle("Dark Respiration - Corrected") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 14, 
                                  hjust = 0)) #set title attributes
Fig.Rd.corr #view plot

Figs.Rd.cold <- arrangeGrob(Fig.Rd.holo,Fig.Rd.skel,Fig.Rd.corr, ncol=3)
ggsave(file="~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Figs.Rd.cold.pdf", Figs.Rd.cold, width = 11, height = 6, units = c("in"), useDingbats=FALSE)

#Gross Photosynthesis (holobiont, skeleton, corrected)
Fig.Pg.holo <-  ggplot(Pgross.means_holo, aes(x=Temp, y=mean,  group=ColorState)) + #set up plot information
  geom_errorbar(aes(x=Temp, ymax=mean+sd, ymin=mean-sd), colour="black", width=.1, position = position_dodge(width = 0.2)) + #add standard error bars about the mean
  geom_point(aes(shape=ColorState), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(17,2,16,1)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Temperature (°C)") + #label x axis
    scale_x_continuous(breaks=c(6,9,12,15,18,22,26,29,32)) +
  ylab(bquote('Gross Photosynthesis ('*mu~ 'mol' ~O[2] ~cm^-2~h^-1*')')) + #label y axis
  ylim(-1, 1.5)+
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(color = 'black', size=10), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_rect(colour='grey95'),
        legend.background = element_rect(colour='black'),
        legend.position=c(.2,.8)) + #set legend location
  ggtitle("Gross Photosythesis - Holobiont") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 14, 
                                  hjust = 0)) #set title attributes
Fig.Pg.holo #view plot

Fig.Pg.skel <-  ggplot(Pgross.means_skel, aes(x=Temp, y=mean,  group=ColorState)) + #set up plot information
  geom_errorbar(aes(x=Temp, ymax=mean+sd, ymin=mean-sd), colour="black", width=.1, position = position_dodge(width = 0.2)) + #add standard error bars about the mean
  geom_point(aes(shape=ColorState), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(17,2,16,1)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Temperature (°C)") + #label x axis
    scale_x_continuous(breaks=c(6,9,12,15,18,22,26,29,32)) +
  ylab(bquote('Gross Photosynthesis ('*mu~ 'mol' ~O[2] ~cm^-2~h^-1*')')) + #label y axis
  ylim(-1, 1.5)+
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(color = 'black', size=10), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_rect(colour='grey95'),
        legend.background = element_rect(colour='black'),
        legend.position="none") + #set legend location
  ggtitle("Gross Photosythesis - Skeleton") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 14, 
                                  hjust = 0)) #set title attributes
Fig.Pg.skel #view plot

Fig.Pg.corr <-  ggplot(Pnet.means_corr, aes(x=Temp, y=mean,  group=ColorState)) + #set up plot information
  geom_errorbar(aes(x=Temp, ymax=mean+sd, ymin=mean-sd), colour="black", width=.1, position = position_dodge(width = 0.2)) + #add standard error bars about the mean
  geom_point(aes(shape=ColorState), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(17,2,16,1)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Temperature (°C)") + #label x axis
    scale_x_continuous(breaks=c(6,9,12,15,18,22,26,29,32)) +
  ylab(bquote('Gross Photosynthesis ('*mu~ 'mol' ~O[2] ~cm^-2~h^-1*')')) + #label y axis
  ylim(-1, 1.5)+
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(color = 'black', size=10), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_rect(colour='grey95'),
        legend.background = element_rect(colour='black'),
        legend.position="none") + #set legend location
  ggtitle("Gross Photosythesis - Corrected") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 14, 
                                  hjust = 0)) #set title attributes
Fig.Pg.corr #view plot

Figs.Pg.cold <- arrangeGrob(Fig.Pg.holo,Fig.Pg.skel,Fig.Pg.corr, ncol=3)
ggsave(file="~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Figs.Pg.cold.pdf", Figs.Pg.cold, width = 11, height = 6, units = c("in"), useDingbats=FALSE)

