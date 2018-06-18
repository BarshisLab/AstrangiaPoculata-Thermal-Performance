## Author: Hannah Aichelman
## Last Updated: March 23, 2018

# This script plots and does statistical analysis of respirometry experiments using a repeated measures ANOVA and post-hoc Tukey's HSD tests. Analysis includes holobiont, skeleton, and corrected data for dark respiration and gross photosynthesis. 

# load in all required packages
library("ggplot2")
library("segmented")
library("plotrix")
library("gridExtra")
library("lubridate")
library("chron")
library("Hmisc")
library("pgirmess")
library("lsmeans")
library("MASS")
library("Rmisc")

#These spreadsheets were created from the master Heat.csv and Cold.csv files that are the output from the Olito script from Hollie Putnam. These are the same spreadsheets being used in Padfield analysis.

holoresp_data <- read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/TPCFitting_Padfield_Analysis/AllHoloRespData_ForPadfield.csv") 

skelresp_data <- read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/TPCFitting_Padfield_Analysis/AllSkelRespData_ForPadfield.csv") #This spreadsheet has been updated from Heat.csv and Cold.csv to include all skeleton respiration rates. Additionally, changed the genotype info so it reflects the true genotype included, not what was originally changed (by me) to make replication even.

corrresp_data <- read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/TPCFitting_Padfield_Analysis/AllCorrRespData_ForPadfield.csv") #This spreadsheet has been updated from Heat.csv and Cold.csv to include the corrected respiration rate. The other data associated with these rates (temp, origin, colorstate, genotype) are all the same as the holobiont data. This disregards the fact that some of the last skeleton ramp samples (ramp4) were not the same genotype for heat and cold ramps. 

photo_data <- read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/TPCFitting_Padfield_Analysis/AllPData_ForPadfield.csv")

skelphoto_data <- read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/TPCFitting_Padfield_Analysis/AllSkelPData_ForPadfield.csv")

#now set working directory to where all of the ANOVA results will be saved
setwd("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/Respirometry_ANOVA/")

###################################################################
## Holobiont respiration data combined, updated March 14, 2018 ####

###Because we are analyzing the heat and cold ramps together, we have twice the amount of 18C measurements as all other temperatures. The following block of code randomly samples the 18C data in half
holo.18.data <- subset(holoresp_data, TempC_holo==18) 

newspoing<-holo.18.data[order(holo.18.data$Genotype_holo),]
randspoing<-newspoing[seq(0,62,by=2)+sample(1:2,32, replace=T),] #randspoing is now a dataset with half of the 18C PAM measurements (aggregated by ID) randomly sampled, n=32 measurements in 18C now

#Create data frame excluding 18C data from all holobiont respiration data
data.no18 <- subset(holoresp_data, TempC_holo != 18)

boingR <- rbind(randspoing,data.no18) #combine 18C data with all other data
write.csv(boingR, "~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/Respirometry_ANOVA/HoloR_RandSample.csv") #write out the randomly sampled file to use again in all future looks at data, this file is included in the 'Data For Figures' folder on github
boingR <- read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/Respirometry_ANOVA/HoloR_RandSample.csv")

boingR$Date_holo<-factor(boingR$Date_holo)

replications(abs_R ~ Genotype_holo, boingR)

replications(abs_R ~ Genotype_holo*Date_holo, boingR)

boxplot(abs_R ~ colorstate + TempC_holo, data=boingR, las=2)

#Bartlett test for homogeneity of variance
for(i in names(table(boingR$TempC_holo))){
	print(paste0("Temp",i))
	print(bartlett.test(abs_R ~ colorstate, data=boingR[boingR$TempC_holo ==i,]))
}

holoR.aov<-aov(abs_R ~ Origin_holo*Color_holo*TempC_holo + Error(Genotype_holo/TempC_holo), data=boingR) 
summary(holoR.aov)

boingR$TempC_holo<-factor(boingR$TempC_holo)
holoR.aov<-aov(abs_R ~ Origin_holo*Color_holo*TempC_holo + Error(Genotype_holo/TempC_holo), data=boingR) 
print(lsmeans(holoR.aov, list(pairwise ~ Origin_holo)), adjust = c("tukey"))
print(lsmeans(holoR.aov, list(pairwise ~ Color_holo)), adjust = c("tukey"))
print(lsmeans(holoR.aov, list(pairwise ~ TempC_holo)), adjust = c("tukey"))
print(lsmeans(holoR.aov, list(pairwise ~ Origin_holo | TempC_holo)), adjust = c("tukey"))
print(lsmeans(holoR.aov, list(pairwise ~  Color_holo | TempC_holo)), adjust = c("tukey"))


# Summarize data using summarySE
all.means_holo <- summarySE(boingR, measurevar="abs_R", groupvars=c("colorstate","TempC_holo"))
names(all.means_holo)[4] <- "mean"

#plot all colorstates data using ggplot
all.means_holo$newtemp<-as.numeric(levels(all.means_holo$TempC_holo)[all.means_holo$TempC_holo])

Fig.holo.R <-  ggplot(all.means_holo, aes(x=newtemp, y=mean,  group=colorstate)) + #set up plot information
  geom_errorbar(aes(x=newtemp, ymax=mean+ci, ymin=mean-ci), colour="black", width=.1, position = position_dodge(width = 0.6)) + #add standard error bars about the mean
  geom_point(aes(shape=colorstate), position = position_dodge(width = 0.6), size=5) + #plot points
  scale_shape_manual(values=c(17,2,16,1)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Temperature (°C)") + #label x axis
    scale_x_continuous(breaks=c(6,9,12,15,18,22,26,29,32)) +
  ylab(bquote('Respiration ('*mu~ 'mol' ~O[2] ~cm^-2~h^-1*')')) + #label y axis
  ylim(-.25, 2.5)+
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
        legend.position=c(.08,.85)) + #set legend location
  ggtitle("Respiration - Holobiont") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 14, 
                                  hjust = 0)) #set title attributes
Fig.holo.R #view plot

ggsave(file="~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/Respirometry_ANOVA/Fig.holo.R.CI.pdf", Fig.holo.R, width = 8, height = 6, units = c("in"), useDingbats=FALSE)

#plot data using base R plot
pdf("DansRHoloPlot.pdf",9,7, useDingbats=FALSE)
par(mar=c(4.5,4.5,1,1))
x.pos<-as.numeric(levels(all.means_holo$TempC_holo)[all.means_holo$TempC_holo])+c(-0.5)
sym=17 #RI brown
indv=1:9
errbar(x.pos[indv], all.means_holo$mean[indv], all.means_holo$mean[indv]+ all.means_holo$ci[indv], all.means_holo$mean[indv]-all.means_holo$ci[indv],  ylim=c(-0.025,2.5), xlim=c(5, 33), pch=sym, type="b", cex=2.5, ylab=bquote('Respiration ('*mu~ 'mol' ~O[2]~cm^-2~h^-1*')'), xlab="Temperature (°C)" , xaxt='n')

sym=2 #RI white
indv=10:18
errbar(x.pos[indv], all.means_holo$mean[indv], all.means_holo$mean[indv]+ all.means_holo$ci[indv], all.means_holo$mean[indv]-all.means_holo$ci[indv],  pch=sym, type="b", cex=2.5,add=TRUE, xaxt='n')

x.pos<-as.numeric(levels(all.means_holo$TempC_holo)[all.means_holo$TempC_holo])+c(0.5)
sym=16 #VA brown
indv=19:27
errbar(x.pos[indv], all.means_holo$mean[indv], all.means_holo$mean[indv]+ all.means_holo$ci[indv], all.means_holo$mean[indv]-all.means_holo$ci[indv],  pch=sym, type="b", cex=2.5,add=TRUE, xaxt='n')

sym=1 #VA white
indv=28:36
errbar(x.pos[indv], all.means_holo$mean[indv], all.means_holo$mean[indv]+ all.means_holo$ci[indv], all.means_holo$mean[indv]-all.means_holo$ci[indv],  pch=sym, type="b", cex=2.5,add=TRUE, xaxt='n')

axis(1,at=c(6,9,12,15,18,22,26,29,32), labels=c(6,9,12,15,18,22,26,29,32))
legend("topleft",c("RI_brown","RI_white","VA_brown","VA_white"), pch=c(17,2,16,1), pt.cex=1.5, bty="n")
text(5,seq(1.9,1.6,by=-0.3/2),c("** Origin, RI > VA","** Color, white > brown", "*** Temp"),adj=c(0,0))
text(c(6,15,18,22,26,29),-.05,'*',cex=2)
dev.off()

##########################################################
## Skeleton Respiration Data ##

###Because we are analyzing the heat and cold ramps together, we have twice the amount of 18C measurements as all other temperatures. The following block of code randomly samples the 18C data in half
skel.18.data <- subset(skelresp_data, Temp_skel==18) 
#write.csv(skel.18.data, "18C_holobiont_data.csv")

newspoing<-skel.18.data[order(skel.18.data$Genotype_skel),]
randspoing<-newspoing[seq(0,62,by=2)+sample(1:2,32, replace=T),] #randspoing is now a dataset with half of the 18C PAM measurements (aggregated by ID) randomly sampled, n=32 measurements in 18C now

#Create data frame excluding 18C data from all skeleton respiration data
data.no18 <- subset(skelresp_data, Temp_skel != 18)
boingR <- rbind(randspoing,data.no18) #combine 18C data with all other data

write.csv(boingR, "~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/Respirometry_ANOVA/SkelR_RandSample.csv") #write out the randomly sampled file to use again in all future looks at data, this file is included in the 'Data For Figures' folder on github
boingR <- read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/Respirometry_ANOVA/SkelR_RandSample.csv")

boingR$Date_skel<-factor(boingR$Date_skel)

replications(abs_R_skel ~ Genotype_skel*Temp_skel, boingR)
replications(abs_R_skel ~ Genotype_skel*Date_skel, boingR)

boxplot(abs_R_skel ~ colorstate + Temp_skel, data=boingR, las=2)

#Bartlett test for homogeneity of variance
for(i in names(table(boingR$Temp_skel))){
	print(paste0("Temp",i))
	print(bartlett.test(abs_R_skel ~ colorstate, data=boingR[boingR$Temp_skel ==i,]))
}

skelR.aov<-aov(abs_R_skel ~ Origin_skel*Color_skel*Temp_skel + Error(Genotype_skel/Temp_skel), data=boingR) 
summary(skelR.aov)

boingR$Temp_skel<-factor(boingR$Temp_skel)
skelR.aov<-aov(abs_R_skel ~ Origin_skel*Color_skel*Temp_skel + Error(Genotype_skel/Temp_skel), data=boingR) 
print(lsmeans(skelR.aov, list(pairwise ~ Temp_skel)), adjust = c("tukey"))


#or use summarySE, this way you have confidence interval to plot later
all.means_skel <- summarySE(boingR, measurevar="abs_R_skel", groupvars=c("colorstate","Temp_skel"))
names(all.means_skel)[4] <- "mean"

#plot all colorstates data using ggplot
all.means_skel$newtemp<-as.numeric(levels(all.means_skel$Temp_skel)[all.means_skel$Temp_skel])

Fig.skel.R <-  ggplot(all.means_skel, aes(x=newtemp, y=mean,  group=colorstate)) + #set up plot information
  geom_errorbar(aes(x=newtemp, ymax=mean+ci, ymin=mean-ci), colour="black", width=.1, position = position_dodge(width = 0.8)) + #add standard error bars about the mean
  geom_point(aes(shape=colorstate), position = position_dodge(width = 0.8), size=5) + #plot points
  scale_shape_manual(values=c(17,2,16,1)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Temperature (°C)") + #label x axis
    scale_x_continuous(breaks=c(6,9,12,15,18,22,26,29,32)) +
  ylab(bquote('Respiration ('*mu~ 'mol' ~O[2] ~cm^-2~h^-1*')')) + #label y axis
  ylim(-.25, 2.5)+
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(color = 'black', size=10), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.title=element_blank(),
        legend.position=c(.08,.85)) + #set legend location
  ggtitle("Respiration - Skeleton") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 14, 
                                  hjust = 0)) #set title attributes
Fig.skel.R #view plot

ggsave(file="~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/Respirometry_ANOVA/Fig.skel.R.CI.pdf", Fig.skel.R, width = 8, height = 6, units = c("in"), useDingbats=FALSE) 

#plot using base R script
pdf("DansRSkelPlot.pdf",9,7, useDingbats=FALSE)
par(mar=c(4.5,4.5,1,1))
x.pos<-as.numeric(levels(all.means_skel$Temp_skel)[all.means_skel$Temp_skel])+c(-0.5)
sym=17 #RI brown
indv=1:9
errbar(x.pos[indv], all.means_skel$mean[indv], all.means_skel$mean[indv]+ all.means_skel$ci[indv], all.means_skel$mean[indv]-all.means_skel$ci[indv],  ylim=c(-0.025,2.5), xlim=c(5, 33), pch=sym, type="b", cex=2.5, ylab=bquote('Respiration ('*mu~ 'mol' ~O[2]~cm^-2~h^-1*')'), xlab="Temperature (°C)" , xaxt='n')

sym=2 #RI white
indv=10:18
errbar(x.pos[indv], all.means_skel$mean[indv], all.means_skel$mean[indv]+ all.means_skel$ci[indv], all.means_skel$mean[indv]-all.means_skel$ci[indv],  pch=sym, type="b", cex=2.5,add=TRUE, xaxt='n')

x.pos<-as.numeric(levels(all.means_skel$Temp_skel)[all.means_skel$Temp_skel])+c(0.5)
sym=16 #VA brown
indv=19:27
errbar(x.pos[indv], all.means_skel$mean[indv], all.means_skel$mean[indv]+ all.means_skel$ci[indv], all.means_skel$mean[indv]-all.means_skel$ci[indv],  pch=sym, type="b", cex=2.5,add=TRUE, xaxt='n')

sym=1 #VA white
indv=28:36
errbar(x.pos[indv], all.means_skel$mean[indv], all.means_skel$mean[indv]+ all.means_skel$ci[indv], all.means_skel$mean[indv]-all.means_skel$ci[indv],  pch=sym, type="b", cex=2.5,add=TRUE, xaxt='n')

axis(1,at=c(6,9,12,15,18,22,26,29,32), labels=c(6,9,12,15,18,22,26,29,32))
legend("topleft",c("RI_brown","RI_white","VA_brown","VA_white"), pch=c(17,2,16,1), pt.cex=1.5, bty="n")
text(5,1.9,c("*** Temp"),adj=c(0,0))
dev.off()

##########################################################
## Corrected Respiration Data (holobiont - skeleton R) ##

corr.18.data <- subset(corrresp_data, TempC_holo==18) 
#write.csv(holo.18.data, "18C_holobiont_data.csv")

newspoing<-corr.18.data[order(corr.18.data$Genotype_holo),]
randspoing<-newspoing[seq(0,62,by=2)+sample(1:2,32, replace=T),] #randspoing is now a dataset with half of the 18C PAM measurements (aggregated by ID) randomly sampled, n=32 measurements in 18C now

#Create data frame excluding 18C data from all holobiont respiration data
data.no18 <- subset(corrresp_data, TempC_holo != 18)

boingR <- rbind(randspoing,data.no18) #combine 18C data with all other data

write.csv(boingR, "~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/Respirometry_ANOVA/CorrR_RandSample.csv") #write out the randomly sampled file to use again in all future looks at data, this file is included in the 'Data For Figures' folder on github
boingR <- read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/Respirometry_ANOVA/CorrR_RandSample.csv")

boingR$Date_holo<-factor(boingR$Date_holo)

replications(abs_Rcorr ~ Genotype_holo*TempC_holo, boingR)

replications(abs_Rcorr ~ Genotype_holo*Date_holo, boingR)

boxplot(abs_Rcorr ~ colorstate + TempC_holo, data=boingR, las=2)

#Bartlett test for homogeneity of variance
for(i in names(table(boingR$TempC_holo))){
	print(paste0("Temp",i))
	print(bartlett.test(abs_Rcorr ~ colorstate, data=boingR[boingR$TempC_holo ==i,]))
}
#Fligner test for homogeneity of variance
for(i in names(table(boingR$TempC_holo))){
	print(paste0("Temp",i))
	print(fligner.test(abs_Rcorr ~ colorstate, data=boingR[boingR$TempC_holo ==i,]))
}

corrR.aov<-aov(abs_Rcorr ~ Origin_holo*Color_holo*TempC_holo + Error(Genotype_holo/TempC_holo), data=boingR) 
summary(corrR.aov)

boingR$TempC_holo<-factor(boingR$TempC_holo)
corrR.aov<-aov(abs_Rcorr ~ Origin_holo*Color_holo*TempC_holo + Error(Genotype_holo/TempC_holo), data=boingR) 
print(lsmeans(corrR.aov, list(pairwise ~ Origin_holo)), adjust = c("tukey"))
print(lsmeans(corrR.aov, list(pairwise ~ Color_holo)), adjust = c("tukey"))
print(lsmeans(corrR.aov, list(pairwise ~ TempC_holo)), adjust = c("tukey"))
print(lsmeans(corrR.aov, list(pairwise ~ Origin_holo | TempC_holo)), adjust = c("tukey"))
#(6C 1.34, 9C 1.46, 12C 1.42, 15C 1.49, 18C 1.50, 22C 1.63, 1.40 26C, 29C 1.17, 32C 1.15) ratio of RI:VA Rcorr at each temperature

# Use summarySE to summarize data
all.means_corr <- summarySE(boingR, measurevar="abs_Rcorr", groupvars=c("colorstate","TempC_holo"))
names(all.means_corr)[4] <- "mean"

#load in all Topt predictions to plot for corrected respiration (Figure 3)
all.means_corr$newtemp<-as.numeric(levels(all.means_corr$TempC_holo)[all.means_corr$TempC_holo])

predsVAbrown <- read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/TPCFitting_Padfield_Analysis/Preds_VAbrown_CorrR.csv")
predsVAbrown$unlogfitted <- exp(1)^predsVAbrown$.fitted
predsVAbrown$unloglwrCI <- exp(1)^predsVAbrown$lwr_CI
predsVAbrown$unloguprCI <- exp(1)^predsVAbrown$upr_CI
predsVAbrown$TempC <- predsVAbrown$TempK_holo-273.15
#predsVAbrown$newtemp<-as.numeric(levels(predsVAbrown$TempC)[predsVAbrown$TempC])

predsVAwhite <- read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/TPCFitting_Padfield_Analysis/Preds_VAwhite_CorrR.csv")
predsVAwhite$unlog <- exp(1)^predsVAwhite$.fitted
predsVAwhite$unlogfitted <- exp(1)^predsVAwhite$.fitted
predsVAwhite$unloglwrCI <- exp(1)^predsVAwhite$lwr_CI
predsVAwhite$unloguprCI <- exp(1)^predsVAwhite$upr_CI
predsVAwhite$TempC <- predsVAwhite$TempK_holo-273.15

predsRIbrown <- read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/TPCFitting_Padfield_Analysis/Preds_RIbrown_CorrR.csv")
predsRIbrown$unlog <- exp(1)^predsRIbrown$.fitted
predsRIbrown$unlogfitted <- exp(1)^predsRIbrown$.fitted
predsRIbrown$unloglwrCI <- exp(1)^predsRIbrown$lwr_CI
predsRIbrown$unloguprCI <- exp(1)^predsRIbrown$upr_CI
predsRIbrown$TempC <- predsRIbrown$TempK_holo-273.15

predsRIwhite <- read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/TPCFitting_Padfield_Analysis/Preds_RIwhite_CorrR.csv")
predsRIwhite$unlog <- exp(1)^predsRIwhite$.fitted
predsRIwhite$unlogfitted <- exp(1)^predsRIwhite$.fitted
predsRIwhite$unloglwrCI <- exp(1)^predsRIwhite$lwr_CI
predsRIwhite$unloguprCI <- exp(1)^predsRIwhite$upr_CI
predsRIwhite$TempC <- predsRIwhite$TempK_holo-273.15

ggplot(predsVAbrown, aes(x=TempC, y=unlogfitted))+
geom_line(aes(x=TempC, y=unlogfitted), colour="black")+
geom_line(aes(x=TempC, y=unloguprCI), colour="red")+
geom_line(aes(x=TempC, y=unloglwrCI), colour="blue")

#plot using ggplot
Fig.corr.R <-  ggplot(all.means_corr, aes(x=newtemp, y=mean,  group=colorstate)) + #set up plot information
  geom_errorbar(aes(x=newtemp, ymax=mean+ci, ymin=mean-ci), colour="black", width=.1, position = position_dodge(width = 0.6)) + #add standard error bars about the mean
  geom_point(aes(shape=colorstate), position = position_dodge(width = 0.6), size=5) + #plot points
  scale_shape_manual(values=c(17,2,16,1)) + #sets point shape manually
  geom_line(predsVAbrown, aes(x=TempC, y=unloguprCI)) +
  #geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines connecting points
  #geom_vline(xintercept=29.201, colour="red", linetype="dotted",size=1)+ #VA white Topt estimate
  #geom_vline(xintercept=29.749, colour="red", linetype="solid",size=1)+ #VA brown Topt estimate
  #geom_vline(xintercept=26.333, colour="blue", linetype="dotted",size=1)+ #RI white Topt estimate
  #geom_vline(xintercept=27.975, colour="blue", linetype="solid",size=1)+ #RI brown Topt estimate
  xlab("Temperature (°C)") + #label x axis
    scale_x_continuous(breaks=c(6,9,12,15,18,22,26,29,32)) +
  ylab(bquote('Respiration ('*mu~ 'mol' ~O[2] ~cm^-2~h^-1*')')) + #label y axis
  ylim(-0.25, 2.5)+
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(color = 'black', size=10), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        axis.text.y = element_text(color = 'black', size=10), 
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.position=c(.08,.85)) + #set legend location
  ggtitle("Respiration - Corrected") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 14, 
                                  hjust = 0)) #set title attributes
Fig.corr.R #view plot

ggsave(file="~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/Respirometry_ANOVA/Fig.corr.R.CI.pdf", Fig.corr.R, width = 8, height = 6, units = c("in"), useDingbats=FALSE)


# Plot using base R script
pdf("DansRCorrPlot.pdf",9,7, useDingbats=FALSE)
par(mar=c(4.5,4.5,1,1))
x.pos<-as.numeric(levels(all.means_corr$TempC_holo)[all.means_corr$TempC_holo])+c(-0.5)
sym=17
indv=1:9
errbar(x.pos[indv], all.means_corr$mean[indv], all.means_corr$mean[indv]+all.means_corr$ci[indv], all.means_corr$mean[indv]-all.means_corr$ci[indv],  ylim=c(-0.025,2.5), xlim=c(5, 33), pch=sym, type="p", cex=2.5, ylab=bquote('Respiration ('*mu~ 'mol' ~O[2]~cm^-2~h^-1*')'), xlab="Temperature (°C)" , xaxt='n')
points(predsRIbrown$TempC, predsRIbrown$unlogfitted, type="l", col="blue", lwd=2)
#points(predsRIbrown$TempC, predsRIbrown$unloglwrCI, type="l", lty=3, col="blue", lwd=2)
#points(predsRIbrown$TempC, predsRIbrown$unloguprCI, type="l", lty=3, col="blue", lwd=2)

sym=2
indv=10:18
errbar(x.pos[indv], all.means_corr$mean[indv], all.means_corr$mean[indv]+all.means_corr$ci[indv], all.means_corr$mean[indv]-all.means_corr$ci[indv],  pch=sym, type="p", cex=2.5,add=TRUE, xaxt='n')
points(predsRIwhite$TempC, predsRIwhite$unlogfitted, type="l", lty=3, col="blue", lwd=2)
#points(predsRIwhite$TempC, predsRIwhite$unloglwrCI, type="l", lty=3, col="blue", lwd=2)
#points(predsRIwhite$TempC, predsRIwhite$unloguprCI, type="l", lty=3, col="blue", lwd=2)

x.pos<-as.numeric(levels(all.means_corr$TempC_holo)[all.means_corr$TempC_holo])+c(0.5)
sym=16
indv=19:27
errbar(x.pos[indv], all.means_corr$mean[indv], all.means_corr$mean[indv]+all.means_corr$ci[indv], all.means_corr$mean[indv]-all.means_corr$ci[indv],  pch=sym, type="p", cex=2.5,add=TRUE, xaxt='n')
points(predsVAbrown$TempC, predsVAbrown$unlogfitted, type="l", col="red", lwd=2)
#points(predsVAbrown$TempC, predsVAbrown$unloglwrCI, type="l", lty=3, col="blue", lwd=2)
#points(predsVAbrown$TempC, predsVAbrown$unloguprCI, type="l", lty=3, col="blue", lwd=2)

sym=1
indv=28:36
errbar(x.pos[indv], all.means_corr$mean[indv], all.means_corr$mean[indv]+all.means_corr$ci[indv], all.means_corr$mean[indv]-all.means_corr$ci[indv],  pch=sym, type="p", cex=2.5,add=TRUE, xaxt='n')
points(predsVAwhite$TempC, predsVAwhite$unlogfitted, type="l", lty=3, col="red", lwd=2)
#points(predsVAwhite$TempC, predsVAwhite$unloglwrCI, type="l", lty=3, col="blue", lwd=2)
#points(predsVAwhite$TempC, predsVAwhite$unloguprCI, type="l", lty=3, col="blue", lwd=2)

axis(1,at=c(6,9,12,15,18,22,26,29,32), labels=c(6,9,12,15,18,22,26,29,32))
legend("topleft",c("RI_brown","RI_white","VA_brown","VA_white"), pch=c(17,2,16,1), pt.cex=1.5, bty="n")
legend(3.5,2.2,c("RI_brown fit","RI_white fit","VA_brown fit","VA_white fit"), lty=c(1,3,1,3), col=c("blue","blue","red","red"),lwd=2, bty="n")
text(5,seq(1,0.7,by=-0.3/2),c("*Origin, RI > VA","* Color, white > brown", "*** Temp"),adj=c(0,0))
text(c(6,15,18,22,26),-.05,'*',cex=2)
dev.off()

###################################################################
## Holobiont photosynthesis data combined, updated March 21, 2018 ####

###Because we are analyzing the heat and cold ramps together, we have twice the amount of 18C measurements as all other temperatures. The following block of code randomly samples the 18C data in half
holo.18.data <- subset(photo_data, TempC_holo==18) 

newspoing<-holo.18.data[order(holo.18.data$Genotype_holo),]
randspoing<-newspoing[seq(0,62,by=2)+sample(1:2,32, replace=T),] #randspoing is now a dataset with half of the 18C PAM measurements (aggregated by ID) randomly sampled, n=32 measurements in 18C now

#Create data frame excluding 18C data from all holobiont respiration data
data.no18 <- subset(photo_data, TempC_holo != 18)

boingR <- rbind(randspoing,data.no18) #combine 18C data with all other data

write.csv(boingR, "~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/Respirometry_ANOVA/HoloP_RandSample.csv") #write out the randomly sampled file to use again in all future looks at data, this file is included in the 'Data For Figures' folder on github
boingR <- read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/Respirometry_ANOVA/HoloP_RandSample.csv")

boingR$Date_holo<-factor(boingR$Date_holo)

replications(Pgross_holo ~ Genotype_holo, boingR)

replications(Pgross_holo ~ Genotype_holo*Date_holo, boingR)

boxplot(Pgross_holo ~ colorstate + TempC_holo, data=boingR, las=2)

#Bartlett test for homogeneity of variance
for(i in names(table(boingR$TempC_holo))){
	print(paste0("Temp",i))
	print(bartlett.test(Pgross_holo ~ colorstate, data=boingR[boingR$TempC_holo ==i,]))
}

holoP.aov<-aov(Pgross_holo ~ Origin_holo*Color_holo*TempC_holo + Error(Genotype_holo/TempC_holo), data=boingR) 
summary(holoP.aov)

boingR$TempC_holo<-factor(boingR$TempC_holo)
holoP.aov<-aov(Pgross_holo ~ Origin_holo*Color_holo*TempC_holo + Error(Genotype_holo/TempC_holo), data=boingR) 
print(lsmeans(holoP.aov, list(pairwise ~ Color_holo)), adjust = c("tukey"))
print(lsmeans(holoP.aov, list(pairwise ~ TempC_holo)), adjust = c("tukey"))
print(lsmeans(holoP.aov, list(pairwise ~ Color_holo | TempC_holo)), adjust = c("tukey"))


#or use summarySE, this way you have confidence interval to plot later
all.means_holo <- summarySE(boingR, measurevar="Pgross_holo", groupvars=c("colorstate","TempC_holo"))
names(all.means_holo)[4] <- "mean"

#plot all colorstates data using ggplot
all.means_holo$newtemp<-as.numeric(levels(all.means_holo$TempC_holo)[all.means_holo$TempC_holo])

Fig.holo.P <-  ggplot(all.means_holo, aes(x=newtemp, y=mean,  group=colorstate)) + #set up plot information
  geom_errorbar(aes(x=newtemp, ymax=mean+ci, ymin=mean-ci), colour="black", width=.1, position = position_dodge(width = 0.6)) + #add standard error bars about the mean
  geom_point(aes(shape=colorstate), position = position_dodge(width = 0.6), size=5) + #plot points
  scale_shape_manual(values=c(17,2,16,1)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Temperature (°C)") + #label x axis
    scale_x_continuous(breaks=c(6,9,12,15,18,22,26,29,32)) +
  ylab(bquote('Pgross ('*mu~ 'mol' ~O[2] ~cm^-2~h^-1*')')) + #label y axis
  ylim(-.5, 1.75)+
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
        legend.position=c(.08,.85)) + #set legend location
  ggtitle("Gross Photosynthesis - Holobiont") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 14, 
                                  hjust = 0)) #set title attributes
Fig.holo.P #view plot

ggsave(file="~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/Respirometry_ANOVA/Fig.holo.Pg.CI.pdf", Fig.holo.P, width = 8, height = 6, units = c("in"), useDingbats=FALSE)

# Plot using base R script
pdf("DansPgHoloPlot_nolines.pdf",9,7)
par(mar=c(4.5,4.5,1,1))
x.pos<-as.numeric(levels(all.means_holo$TempC_holo)[all.means_holo$TempC_holo])+c(-0.5)
sym=17
indv=1:9 #RI brown
errbar(x.pos[indv], all.means_holo$mean[indv], all.means_holo$mean[indv]+ all.means_holo$ci[indv], all.means_holo$mean[indv]-all.means_holo$ci[indv],  ylim=c(-.5, 1.75), xlim=c(5, 33), pch=sym, type="p", cex=2.5, ylab=bquote('Gross Photosynthesis ('*mu~ 'mol' ~O[2]~cm^-2~h^-1*')'), xlab="Temperature (°C)" , xaxt='n') #type = b to connect dots

sym=2
indv=10:18 #RI white
errbar(x.pos[indv], all.means_holo$mean[indv], all.means_holo$mean[indv]+ all.means_holo$ci[indv], all.means_holo$mean[indv]-all.means_holo$ci[indv],  pch=sym, type="p", cex=2.5,add=TRUE, xaxt='n')

x.pos<-as.numeric(levels(all.means_holo$TempC_holo)[all.means_holo$TempC_holo])+c(0.5)
sym=16
indv=19:27 #VA brown
errbar(x.pos[indv], all.means_holo$mean[indv], all.means_holo$mean[indv]+ all.means_holo$ci[indv], all.means_holo$mean[indv]-all.means_holo$ci[indv],  pch=sym, type="p", cex=2.5,add=TRUE, xaxt='n')

sym=1
indv=28:36 #VA white
errbar(x.pos[indv], all.means_holo$mean[indv], all.means_holo$mean[indv]+ all.means_holo$ci[indv], all.means_holo$mean[indv]-all.means_holo$ci[indv],  pch=sym, type="p", cex=2.5,add=TRUE, xaxt='n')

axis(1,at=c(6,9,12,15,18,22,26,29,32), labels=c(6,9,12,15,18,22,26,29,32))
legend("topleft",c("RI_brown","RI_white","VA_brown","VA_white"), pch=c(17,2,16,1), pt.cex=1.5, bty="n")
text(5,seq(1.2,0.8,by=-0.3/2),c("*** Color, white > brown", "*** Temp", "** Color*Temp"),adj=c(0,0))
text(c(6,12,15,18,22,26,29,32),-.52,'+',cex=2)
dev.off()


###################################################################
## Skeleton photosynthesis data combined, updated March 21, 2018 ####

###Because we are analyzing the heat and cold ramps together, we have twice the amount of 18C measurements as all other temperatures. The following block of code randomly samples the 18C data in half
skel.18.data <- subset(skelphoto_data, Temp_skel==18) 

newspoing<-skel.18.data[order(skel.18.data$Genotype_skel),]
randspoing<-newspoing[seq(0,62,by=2)+sample(1:2,32, replace=T),] #randspoing is now a dataset with half of the 18C PAM measurements (aggregated by ID) randomly sampled, n=32 measurements in 18C now

#Create data frame excluding 18C data from all holobiont respiration data
data.no18 <- subset(skelphoto_data, Temp_skel != 18)

boingR <- rbind(randspoing,data.no18) #combine 18C data with all other data

write.csv(boingR, "~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/Respirometry_ANOVA/SkelP_RandSample.csv") ##write out the randomly sampled file to use again in all future looks at data, this file is included in the 'Data For Figures' folder on github
boingR <- read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/Respirometry_ANOVA/SkelP_RandSample.csv")

boingR$Date_skel<-factor(boingR$Date_skel)

replications(Pgross_skel ~ Genotype_skel*Temp_skel, boingR)

replications(Pgross_skel ~ Genotype_skel*Date_skel, boingR)

boxplot(Pgross_skel ~ colorstate + Temp_skel, data=boingR, las=2)

#Bartlett test for homogeneity of variance
for(i in names(table(boingR$Temp_skel))){
	print(paste0("Temp",i))
	print(bartlett.test(Pgross_skel ~ colorstate, data=boingR[boingR$Temp_skel ==i,]))
}

skelP.aov<-aov(Pgross_skel ~ Origin_skel*Color_skel*Temp_skel + Error(Genotype_skel/Temp_skel), data=boingR) 
summary(skelP.aov)

boingR$Temp_skel<-factor(boingR$Temp_skel)
skelP.aov<-aov(Pgross_skel ~ Origin_skel*Color_skel*Temp_skel + Error(Genotype_skel/Temp_skel), data=boingR) 
print(lsmeans(skelP.aov, list(pairwise ~ Temp_skel)), adjust = c("tukey"))

#or use summarySE, this way you have confidence interval to plot later
all.means_skel <- summarySE(boingR, measurevar="Pgross_skel", groupvars=c("colorstate","Temp_skel"))
names(all.means_skel)[4] <- "mean"

#plot all colorstates data using ggplot
all.means_skel$newtemp<-as.numeric(levels(all.means_skel$Temp_skel)[all.means_skel$Temp_skel])

Fig.skel.P <-  ggplot(all.means_skel, aes(x=newtemp, y=mean,  group=colorstate)) + #set up plot information
  geom_errorbar(aes(x=newtemp, ymax=mean+ci, ymin=mean-ci), colour="black", width=.1, position = position_dodge(width = 0.8)) + #add standard error bars about the mean
  geom_point(aes(shape=colorstate), position = position_dodge(width = 0.8), size=5) + #plot points
  scale_shape_manual(values=c(17,2,16,1)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Temperature (°C)") + #label x axis
    scale_x_continuous(breaks=c(6,9,12,15,18,22,26,29,32)) +
  ylab(bquote('Pgross ('*mu~ 'mol' ~O[2] ~cm^-2~h^-1*')')) + #label y axis
  ylim(-.58, 1.75)+
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
        legend.position=c(.08,.85)) + #set legend location
  ggtitle("Gross Photosynthesis - Skeleton") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 14, 
                                  hjust = 0)) #set title attributes
Fig.skel.P #view plot

ggsave(file="~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/Respirometry_ANOVA/Fig.skel.Pg.CI.pdf", Fig.skel.P, width = 8, height = 6, units = c("in"), useDingbats=FALSE)


#Plot using base R script
pdf("DansPgSkelPlot.pdf",9,7, useDingbats=FALSE)
par(mar=c(4.5,4.5,1,1))
x.pos<-as.numeric(levels(all.means_skel$Temp_skel)[all.means_skel$Temp_skel])+c(-0.5)
sym=17 #RI brown
indv=1:9
errbar(x.pos[indv], all.means_skel$mean[indv], all.means_skel$mean[indv]+ all.means_skel$ci[indv], all.means_skel$mean[indv]-all.means_skel$ci[indv],  ylim=c(-.6, 1.75), xlim=c(5, 33), pch=sym, type="b", cex=2.5, ylab=bquote('Gross Photosynthesis ('*mu~ 'mol' ~O[2]~cm^-2~h^-1*')'), xlab="Temperature (°C)" , xaxt='n')

sym=2 #RI white
indv=10:18
errbar(x.pos[indv], all.means_skel$mean[indv], all.means_skel$mean[indv]+ all.means_skel$ci[indv], all.means_skel$mean[indv]-all.means_skel$ci[indv],  pch=sym, type="b", cex=2.5,add=TRUE, xaxt='n')

x.pos<-as.numeric(levels(all.means_skel$Temp_skel)[all.means_skel$Temp_skel])+c(0.5)
sym=16 #VA brown
indv=19:27
errbar(x.pos[indv], all.means_skel$mean[indv], all.means_skel$mean[indv]+ all.means_skel$ci[indv], all.means_skel$mean[indv]-all.means_skel$ci[indv],  pch=sym, type="b", cex=2.5,add=TRUE, xaxt='n')

sym=1 #VA white
indv=28:36
errbar(x.pos[indv], all.means_skel$mean[indv], all.means_skel$mean[indv]+ all.means_skel$ci[indv], all.means_skel$mean[indv]-all.means_skel$ci[indv],  pch=sym, type="b", cex=2.5,add=TRUE, xaxt='n')

axis(1,at=c(6,9,12,15,18,22,26,29,32), labels=c(6,9,12,15,18,22,26,29,32))
legend("topleft",c("RI_brown","RI_white","VA_brown","VA_white"), pch=c(17,2,16,1), pt.cex=1.5, bty="n")
text(5,1.0,c("* Temp"),adj=c(0,0))
dev.off()


###################################################################
## Corrected photosynthesis data combined, updated March 21, 2018 ####

###Because we are analyzing the heat and cold ramps together, we have twice the amount of 18C measurements as all other temperatures. The following block of code randomly samples the 18C data in half
corr.18.data <- subset(photo_data, TempC_holo==18) 

newspoing<-corr.18.data[order(corr.18.data$Genotype_holo),]
randspoing<-newspoing[seq(0,62,by=2)+sample(1:2,32, replace=T),] #randspoing is now a dataset with half of the 18C PAM measurements (aggregated by ID) randomly sampled, n=32 measurements in 18C now

#Create data frame excluding 18C data from all holobiont respiration data
data.no18 <- subset(photo_data, TempC_holo != 18)

boingR <- rbind(randspoing,data.no18) #combine 18C data with all other data

write.csv(boingR, "~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/Respirometry_ANOVA/CorrP_RandSample.csv") #write out the randomly sampled file to use again in all future looks at data, this file is included in the 'Data For Figures' folder on github
boingR <- read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/Respirometry_ANOVA/CorrP_RandSample.csv")

boingR$Date_holo<-factor(boingR$Date_holo)

replications(Pgross_corr ~ Genotype_holo*TempC_holo, boingR)

replications(Pgross_corr ~ Genotype_holo*Date_holo, boingR)

boxplot(Pgross_corr ~ colorstate + TempC_holo, data=boingR, las=2)

#Bartlett test for homogeneity of variance
for(i in names(table(boingR$TempC_holo))){
	print(paste0("Temp",i))
	print(bartlett.test(Pgross_corr ~ colorstate, data=boingR[boingR$TempC_holo ==i,]))
}

corrP.aov<-aov(Pgross_corr ~ Origin_holo*Color_holo*TempC_holo + Error(Genotype_holo/TempC_holo), data=boingR) 
summary(corrP.aov)

boingR$TempC_holo<-factor(boingR$TempC_holo)
corrP.aov<-aov(Pgross_corr ~ Origin_holo*Color_holo*TempC_holo + Error(Genotype_holo/TempC_holo), data=boingR) 
print(lsmeans(corrP.aov, list(pairwise ~ Color_holo)), adjust = c("tukey"))
print(lsmeans(corrP.aov, list(pairwise ~ TempC_holo)), adjust = c("tukey"))
print(lsmeans(corrP.aov, list(pairwise ~ Color_holo | TempC_holo)), adjust = c("tukey"))

#or use summarySE, this way you have confidence interval to plot later
all.means_corr <- summarySE(boingR, measurevar="Pgross_corr", groupvars=c("colorstate","TempC_holo"))
names(all.means_corr)[4] <- "mean"

#plot all colorstates data using ggplot
all.means_corr$newtemp<-as.numeric(levels(all.means_corr$TempC_holo)[all.means_corr$TempC_holo])

Fig.corr.P <-  ggplot(all.means_corr, aes(x=newtemp, y=mean,  group=colorstate)) + #set up plot information
  geom_errorbar(aes(x=newtemp, ymax=mean+ci, ymin=mean-ci), colour="black", width=.1, position = position_dodge(width = 0.6)) + #add standard error bars about the mean
  geom_point(aes(shape=colorstate), position = position_dodge(width = 0.6), size=5) + #plot points
  scale_shape_manual(values=c(17,2,16,1)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Temperature (°C)") + #label x axis
    scale_x_continuous(breaks=c(6,9,12,15,18,22,26,29,32)) +
  ylab(bquote('Pgross ('*mu~ 'mol' ~O[2] ~cm^-2~h^-1*')')) + #label y axis
  	scale_y_continuous(breaks=c(-0.5,0.0,0.5,1.0,1.5))+
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
        legend.position=c(.08,.85)) + #set legend location
  ggtitle("Gross Photosynthesis - Corrected") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 14, 
                                  hjust = 0)) #set title attributes

Fig.corr.P #view plot

ggsave(file="~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/Respirometry_ANOVA/Fig.corr.Pg.CI.pdf", Fig.corr.P, width = 8, height = 6, units = c("in"), useDingbats=FALSE)

# Plot using base R script
pdf("DansPgCorrPlot.pdf",9,7, useDingbats=FALSE)
par(mar=c(4.5,4.5,1,1))
x.pos<-as.numeric(levels(all.means_corr$TempC_holo)[all.means_corr$TempC_holo])+c(-0.5)
sym=17
indv=1:9
errbar(x.pos[indv], all.means_corr$mean[indv], all.means_corr$mean[indv]+all.means_corr$ci[indv], all.means_corr$mean[indv]-all.means_corr$ci[indv],  ylim=c(-.6, 1.75), xlim=c(5, 33), pch=sym, type="b", cex=2.5, ylab=bquote('Gross Photosynthesis ('*mu~ 'mol' ~O[2]~cm^-2~h^-1*')'), xlab="Temperature (°C)" , xaxt='n')

sym=2
indv=10:18
errbar(x.pos[indv], all.means_corr$mean[indv], all.means_corr$mean[indv]+all.means_corr$ci[indv], all.means_corr$mean[indv]-all.means_corr$ci[indv],  pch=sym, type="b", cex=2.5,add=TRUE, xaxt='n')

x.pos<-as.numeric(levels(all.means_corr$TempC_holo)[all.means_corr$TempC_holo])+c(0.5)
sym=16
indv=19:27
errbar(x.pos[indv], all.means_corr$mean[indv], all.means_corr$mean[indv]+all.means_corr$ci[indv], all.means_corr$mean[indv]-all.means_corr$ci[indv],  pch=sym, type="b", cex=2.5,add=TRUE, xaxt='n')

sym=1
indv=28:36
errbar(x.pos[indv], all.means_corr$mean[indv], all.means_corr$mean[indv]+all.means_corr$ci[indv], all.means_corr$mean[indv]-all.means_corr$ci[indv],  pch=sym, type="b", cex=2.5,add=TRUE, xaxt='n')

axis(1,at=c(6,9,12,15,18,22,26,29,32), labels=c(6,9,12,15,18,22,26,29,32))
legend("topleft",c("RI_brown","RI_white","VA_brown","VA_white"), pch=c(17,2,16,1), pt.cex=1.5, bty="n")
text(5,seq(1.2,0.8,by=-0.3/2),c("*** Color, brown > white", "** Temp", "* Color*Temp"),adj=c(0,0))
text(c(6,22,26,29),-.64,'+',cex=2)
dev.off()

###################################################################
