#Written by Hannah Aichelman
#Last updated March 2, 2018 
#Script to analyze and plot Fv/Fm data from respirometry experiments, including holobiont Fv/Fm (Fig 5) and skeleton Fv/Fm (Fig S4A)

library("ggplot2")
library("segmented")
library("plotrix")
library("gridExtra")
library("lubridate")
library("chron")
library("Hmisc")
library("Rmisc")
library("pgirmess")
library("lsmeans")
library("MASS")

#for all experiments, FvFm values listed as - in the spreadsheet were converted to 0 before reading data in. 
#load in data from heat experiments

setwd('~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/PAM/Cleaned_DataFiles')
heat1<-read.csv('20170707_PAM_cleaned.csv',header=TRUE,sep=",")
heat2<-read.csv('20170708_PAM_cleaned.csv',header=TRUE,sep=",")
heat3<-read.csv('20170710_PAM_cleaned.csv',header=TRUE,sep=",")
heat4<-read.csv('20170711_PAM_cleaned.csv',header=TRUE,sep=",")
#load in data from cold experiments
cold1<-read.csv('20170712_PAM_cleaned.csv',header=TRUE,sep=",")
cold2<-read.csv('20170713_PAM_cleaned.csv',header=TRUE,sep=",")
cold3<-read.csv('20170714_PAM_cleaned.csv',header=TRUE,sep=",")
cold4<-read.csv('20170716_PAM_cleaned.csv',header=TRUE,sep=",")
#load in data from heat skeleton experiments
heatskel1<-read.csv('20170717_PAM_cleaned.csv',header=TRUE,sep=",")
heatskel2<-read.csv('20170718_PAM_cleaned.csv',header=TRUE,sep=",")
heatskel3<-read.csv('20170719_PAM_cleaned.csv',header=TRUE,sep=",")
heatskel4<-read.csv('20170721_PAM_cleaned.csv',header=TRUE,sep=",")
#load in data from cold skeleton experiments
coldskel1<-read.csv('20170722_PAM_cleaned.csv',header=TRUE,sep=",")
coldskel2<-read.csv('20170724_PAM_cleaned.csv',header=TRUE,sep=",")
coldskel3<-read.csv('20170725_PAM_cleaned.csv',header=TRUE,sep=",")
coldskel4<-read.csv('20170726_PAM_cleaned.csv',header=TRUE,sep=",")


#set working directory back to where plots and .csv files will be saved
setwd('~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/PAM')

###################################################################
#### analyze all holobiont data combined, updated March 16, 2018
# updated March 16th so you don't have to write out or read in files, can do all data manipulation inside this script (mostly random sampling of 18C measurements)

holo.data <- rbind(heat1,heat2,heat3,heat4,cold1,cold2,cold3,cold4)
holo.data <-cbind(holo.data[0:11],"colorstate"=paste(holo.data $origin, holo.data $color, sep="_"), holo.data[12]) 

#remove all 18C data so we can later randomly subset
holo.18.data <- subset(holo.data, temp==18) 

#randomly subset 18C data
spoing<-rbind(aggregate(FvFm ~ ID + temp + Date + color + origin + ramptype + genotype, data= holo.18.data, mean)) ###Average replicate measurements for each individual by temp and day

newspoing<-spoing[order(spoing$genotype),]
randspoing<-newspoing[seq(0,62,by=2)+sample(1:2,32, replace=T),] #randspoing is now a dataset with half of the 18C PAM measurements (aggregated by ID) randomly sampled, n=32 measurements in 18C now

table(randspoing$genotype)
table(randspoing$ID)

#Subset all data except for 18C measurements
data.no18 <- subset(holo.data, temp != 18)

data.no18<-aggregate(FvFm ~ ID + temp + Date + color + origin + ramptype + genotype, data=data.no18, mean) ###Average replicate measurements for each individual by temp and day

boingPAM <- rbind(randspoing,data.no18) #combine 18C data with all other data
boingPAM <-cbind(boingPAM[0:7],"colorstate"=paste(boingPAM $origin, boingPAM $color, sep="_"), boingPAM[8])  #add colorstate factor

write.csv(boingPAM, "~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/PAM/FvFm_RandSample_holo.csv") 

#once the data has been randomly sampled, read in for all future looks at the data
boingPAM<-read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/PAM/FvFm_RandSample_holo.csv")
head(boingPAM)

boingPAM$Date<-factor(boingPAM$Date)

#check replication of genotype x temp data for holobiont
replications(FvFm ~ genotype*temp, boingPAM)
replications(FvFm ~ genotype*Date, boingPAM)


aggregate(FvFm ~ color + origin + temp, data=boingPAM, FUN=function(x) shapiro.test(x)$p.value) #shapiro test assesses normality of data. 

boxplot(FvFm ~ colorstate + temp, data=boingPAM, las=2)

#Bartlett test for homogeneity of variance
for(i in names(table(boingPAM$temp))){
	print(paste0("Temp",i))
	print(bartlett.test(FvFm ~ colorstate, data=boingPAM[boingPAM$temp==i,]))
}

#ANOVA
holo.aov<-aov(FvFm ~ origin*color*temp + Error(genotype/temp), data=boingPAM) 
summary(holo.aov)

boingPAM$temp<-factor(boingPAM$temp)
holo.aov<-aov(FvFm ~ origin*color*temp + Error(genotype/temp), data=boingPAM) 
print(lsmeans(holo.aov, list(pairwise ~ temp)), adjust = c("tukey"))
print(lsmeans(holo.aov, list(pairwise ~ color)), adjust = c("tukey"))
print(lsmeans(holo.aov, list(pairwise ~ color | temp)), adjust = c("tukey"))

#Summarize data to plot
all.means_holo <- summarySE(boingPAM, measurevar="FvFm", groupvars=c("colorstate","temp"))
names(all.means_holo)[4] <- "mean"

#plot all colorstates data using ggplot
all.means_holo$newtemp<-as.numeric(levels(all.means_holo$temp)[all.means_holo$temp])

Fig.holo.FvFm <-  ggplot(all.means_holo, aes(x=newtemp, y=mean,  group=colorstate)) + #set up plot information
  geom_errorbar(aes(x=newtemp, ymax=mean+ci, ymin=mean-ci), colour="black", width=.1, position = position_dodge(width = 0.6)) + #add standard error bars about the mean
  geom_point(aes(shape=colorstate), position = position_dodge(width = 0.6), size=5) + #plot points
  scale_shape_manual(values=c(17,2,16,1)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Temperature (°C)") + #label x axis
    scale_x_continuous(breaks=c(6,9,12,15,18,22,26,29,32)) +
  ylab("Fv/Fm") + #label y axis
  ylim(-0.1, 0.7)+
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(color = 'black', size=10), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.position=c(.08,.88)) + #set legend location
  ggtitle("Fv/Fm - Holobiont") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 14, 
                                  hjust = 0)) #set title attributes
Fig.holo.FvFm #view plot

ggsave(file="~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/PAM/FvFm_Holobiont.CI.pdf", Fig.holo.FvFm, width = 8, height = 6, units = c("in"), useDingbats=FALSE)

# or create the plot using base R plot
pdf("DansFvFmHoloPlot_nolines.pdf",9,7, useDingbats=FALSE)
par(mar=c(4.5,4.5,1,1))
x.pos<-as.numeric(levels(all.means_holo$temp)[all.means_holo$temp])+c(-0.5)
sym=17 #RI brown
indv=1:9
errbar(x.pos[indv], all.means_holo$mean[indv], all.means_holo$mean[indv]+ all.means_holo$ci[indv], all.means_holo$mean[indv]-all.means_holo$ci[indv],  ylim=c(-0.1, 0.7), xlim=c(5, 33), pch=sym, type="p", cex=2.5, ylab="Fv/Fm", xlab="Temperature (°C)" , xaxt='n') #type = "p" for no lines, type = "b" for lines connecting points

sym=2 #RI white
indv=10:18
errbar(x.pos[indv], all.means_holo$mean[indv], all.means_holo$mean[indv]+ all.means_holo$ci[indv], all.means_holo$mean[indv]-all.means_holo$ci[indv],  pch=sym, type="p", cex=2.5,add=TRUE, xaxt='n')

x.pos<-as.numeric(levels(all.means_holo$temp)[all.means_holo$temp])+c(0.5)
sym=16 #VA brown
indv=19:27
errbar(x.pos[indv], all.means_holo$mean[indv], all.means_holo$mean[indv]+ all.means_holo$ci[indv], all.means_holo$mean[indv]-all.means_holo$ci[indv],  pch=sym, type="p", cex=2.5,add=TRUE, xaxt='n')

sym=1 #VA white
indv=28:36
errbar(x.pos[indv], all.means_holo$mean[indv], all.means_holo$mean[indv]+ all.means_holo$ci[indv], all.means_holo$mean[indv]-all.means_holo$ci[indv],  pch=sym, type="p", cex=2.5,add=TRUE, xaxt='n')

axis(1,at=c(6,9,12,15,18,22,26,29,32), labels=c(6,9,12,15,18,22,26,29,32))
legend("topleft",c("RI_brown","RI_white","VA_brown","VA_white"), pch=c(17,2,16,1), pt.cex=1.5, bty="n")
text(10,seq(0.68,0.61,by=-0.033),c("*** Color, brown > white", "** Temp", "** Color*Temp"),adj=c(0,0))
dev.off()

#### analyze all skeleton data combined (Fig S4A)

skel.data <- rbind(heatskel1,heatskel2,heatskel3,heatskel4,coldskel1,coldskel2,coldskel3,coldskel4)
skel.data <-cbind(skel.data[0:11],"colorstate"=paste(skel.data $origin, skel.data $color, sep="_"), skel.data[12]) 

#remove all 18C data so we can later randomly subset
skel.18.data <- subset(skel.data, temp==18) 

#randomly subset 18C data
spoing<-rbind(aggregate(FvFm ~ ID + temp + Date + color + origin + ramptype + genotype, data= skel.18.data, mean)) ###Average replicate measurements for each individual by temp and day

newspoing<-spoing[order(spoing$genotype),]
randspoing<-newspoing[seq(0,62,by=2)+sample(1:2,32, replace=T),] #randspoing is now a dataset with half of the 18C PAM measurements (aggregated by ID) randomly sampled, n=32 measurements in 18C now

table(randspoing$genotype)
table(randspoing$ID)

#read in all other holobiont data, excluding 18C data
data.no18 <- subset(skel.data, temp != 18)

data.no18<-aggregate(FvFm ~ ID + temp + Date + color + origin + ramptype + genotype, data=data.no18, mean) ###Average replicate measurements for each individual by temp and day

boingPAM <- rbind(randspoing,data.no18) #combine 18C data with all other data
boingPAM <-cbind(boingPAM[0:7],"colorstate"=paste(boingPAM $origin, boingPAM $color, sep="_"), boingPAM[8])  #add colorstate factor

write.csv(boingPAM, "~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/PAM/FvFm_RandSample_skel.csv") 
#once the data has been randomly sampled, read in for all future looks at the data
boingPAM<-read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/PAM/FvFm_RandSample_skel.csv")

boingPAM$Date<-factor(boingPAM$Date)

#check replications of skeleton genotype and temperature data. 
replications(FvFm ~ genotype*temp, boingPAM)

replications(FvFm ~ genotype*Date, boingPAM)

aggregate(FvFm ~ color + origin + temp, data=boingPAM, FUN=function(x) shapiro.test(x)$p.value) #shapiro test assesses normality of data. 

boxplot(FvFm ~ colorstate + temp, data=boingPAM, las=2)

#Bartlett test for homogeneity of variance
for(i in names(table(boingPAM$temp))){
	print(paste0("Temp",i))
	print(bartlett.test(FvFm ~ colorstate, data=boingPAM[boingPAM$temp==i,]))
}

#ANOVA
skel.aov<-aov(FvFm ~ origin*color*temp + Error(genotype/temp), data=boingPAM) 
summary(skel.aov)

boingPAM$temp<-factor(boingPAM$temp)
skel.aov<-aov(FvFm ~ origin*color*temp + Error(genotype/temp), data=boingPAM) 
print(lsmeans(skel.aov, list(pairwise ~ origin)), adjust = c("tukey"))
print(lsmeans(skel.aov, list(pairwise ~ color)), adjust = c("tukey"))
print(lsmeans(skel.aov, list(pairwise ~ temp)), adjust = c("tukey"))
print(lsmeans(skel.aov, list(pairwise ~ origin | color)), adjust = c("tukey"))
print(lsmeans(skel.aov, list(pairwise ~  origin | temp | color)), adjust = c("tukey"))

#Summarize data to plot
all.means_skel <- summarySE(boingPAM, measurevar="FvFm", groupvars=c("colorstate","temp"))
names(all.means_skel)[4] <- "mean"

#plot using ggplot
all.means_skel$newtemp<-as.numeric(levels(all.means_skel$temp)[all.means_skel$temp])

Fig.skel.FvFm <-  ggplot(all.means_skel, aes(x=newtemp, y=mean,  group=colorstate)) + #set up plot information
  geom_errorbar(aes(x=newtemp, ymax=mean+ci, ymin=mean-ci), colour="black", width=.1, position = position_dodge(width = 0.6)) + #add standard error bars about the mean
  geom_point(aes(shape=colorstate), position = position_dodge(width = 0.6), size=5) + #plot points
  scale_shape_manual(values=c(17,2,16,1)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Temperature (°C)") + #label x axis
    scale_x_continuous(breaks=c(6,9,12,15,18,22,26,29,32)) +
  ylab("Fv/Fm") + #label y axis
  ylim(-.2, 0.7)+
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(color = 'black', size=10), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.position=c(.08,.88)) + #set legend location
  ggtitle("Fv/Fm - Skeleton") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 14, 
                                  hjust = 0)) #set title attributes
Fig.skel.FvFm #view plot

ggsave(file="~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/PAM/FvFm_Skeleton.CI.pdf", Fig.skel.FvFm, width = 8, height = 6, units = c("in"), useDingbats=FALSE)

#or plot using base R script
pdf("DansFvFmSkelPlot.pdf",9,7, useDingbats=FALSE)
par(mar=c(4.5,4.5,1,1))
x.pos<-as.numeric(levels(all.means_skel$temp)[all.means_skel$temp])+c(-0.5)
sym=17 #RI brown
indv=1:9
errbar(x.pos[indv], all.means_skel$mean[indv], all.means_skel$mean[indv]+ all.means_skel$ci[indv], all.means_skel$mean[indv]-all.means_skel$ci[indv],  ylim=c(-0.1, 0.7), xlim=c(5, 33), pch=sym, type="b", cex=2.5, ylab="Fv/Fm", xlab="Temperature (°C)" , xaxt='n')

sym=2 #RI white
indv=10:18
errbar(x.pos[indv], all.means_skel$mean[indv], all.means_skel$mean[indv]+ all.means_skel$ci[indv], all.means_skel$mean[indv]-all.means_skel$ci[indv],  pch=sym, type="b", cex=2.5,add=TRUE, xaxt='n')

x.pos<-as.numeric(levels(all.means_skel$temp)[all.means_skel$temp])+c(0.5)
sym=16 #VA brown
indv=19:27
errbar(x.pos[indv], all.means_skel$mean[indv], all.means_skel$mean[indv]+ all.means_skel$ci[indv], all.means_skel$mean[indv]-all.means_skel$ci[indv],  pch=sym, type="b", cex=2.5,add=TRUE, xaxt='n')

sym=1 #VA white
indv=28:36
errbar(x.pos[indv], all.means_skel$mean[indv], all.means_skel$mean[indv]+ all.means_skel$ci[indv], all.means_skel$mean[indv]-all.means_skel$ci[indv],  pch=sym, type="b", cex=2.5,add=TRUE, xaxt='n')

axis(1,at=c(6,9,12,15,18,22,26,29,32), labels=c(6,9,12,15,18,22,26,29,32))
legend("topleft",c("RI_brown","RI_white","VA_brown","VA_white"), pch=c(17,2,16,1), pt.cex=1.5, bty="n")
text(10,seq(0.7,0.6,by=-0.033),c("** Origin, RI > VA","* Color, brown > white", "* Origin*Color", "* Origin*Color*Temp"),adj=c(0,0))
dev.off()



