#Author: Hannah Aichelman
#Date Last Updated: February 23, 2018
#This script analyzes and plots red channel intensity data calculated from photos of corals from the Winters et al. 2009 Matlab script

getwd() 
setwd("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/RGBAnalysis") #set working directory

#library all packages needed for the code
library(Hmisc)
library(plotrix)
library(lsmeans)
library(Rmisc)

data=read.csv('RedChannelIntensity_data.csv') #read the data into R

head(data) #take a look at the first few lines of data
#create a new column of data called "colorstate", which is combined origin and color. we use this later to group the data for plotting.
newdata<-cbind(data[0:4],"colorstate"=paste(data$Origin,data$Color, sep="_"),data[5:ncol(data)])  
head(newdata)


#subset data into data from holobiont and skeleton photos to analyze separately
holobiontdata <- subset(newdata, PhotoType=='holobiont')
holobiontdata

skeletondata <- subset(newdata, PhotoType=='skeleton')
skeletondata


# ANOVAS

# difference in holobiont color based on origin and color of coral
holobiontdata <- subset(holobiontdata,PhotoTime=="BeforeRamp") #including only holobiont data before ramp

#Bartlett test for homogeneity of variance
bartlett.test(AvgRedChannel_reanalysis ~ colorstate, data=holobiontdata)

color.aov <- aov(AvgRedChannel_reanalysis ~ Origin*Color, data=holobiontdata)
summary(color.aov)
print(lsmeans(color.aov, list(pairwise ~ Origin)), adjust = c("tukey"))
print(lsmeans(color.aov, list(pairwise ~ Color)), adjust = c("tukey"))

# difference in skeleton color based on origin and color, including only skeleton data before the ramp
#Bartlett test for homogeneity of variance
bartlett.test(AvgRedChannel_reanalysis ~ colorstate, data=skeletondata)

skeleton.aov <- aov(AvgRedChannel_reanalysis ~ Origin*Color, data = skeletondata)
summary(skeleton.aov)
print(lsmeans(skeleton.aov, list(pairwise ~ Origin)), adjust = c("tukey"))
print(lsmeans(skeleton.aov, list(pairwise ~ Origin | Color)), adjust = c("tukey"))


# difference based on time of photo (pre vs. post experiment)
holo.skel.aov <- aov(AvgRedChannel_reanalysis ~ PhotoTime, data=newdata)
summary(holo.skel.aov)
print(lsmeans(holo.skel.aov, list(pairwise ~ PhotoTime)), adjust = c("tukey"))

# difference between holobiont and skeleton before ramp
beforeramp <- subset(newdata,PhotoTime=="BeforeRamp" | PhotoTime=="AfterAirbrush")
holo.skel.aov <- aov(AvgRedChannel_reanalysis ~ PhotoType, data=beforeramp)
summary(holo.skel.aov)
print(lsmeans(holo.skel.aov, list(pairwise ~ PhotoType)), adjust = c("tukey"))


### Analyze and plot all data (holobiont and skeleton) combined ###

Color.mean <- aggregate(AvgRedChannel_reanalysis ~ colorstate*Origin*PhotoType, data=newdata, mean)
Color.se <- aggregate(AvgRedChannel_reanalysis ~ colorstate*Origin*PhotoType, data=newdata, std.error)
Color.n <- aggregate(AvgRedChannel_reanalysis ~ colorstate*Origin*PhotoType, data= newdata, length)
Color.means <- cbind(Color.mean, Color.se$AvgRedChannel_reanalysis, Color.n$AvgRedChannel_reanalysis) #combine mean and standard error results
colnames(Color.means) <- c("ColorState", "Origin", "PhotoType","mean", "se", "n")   #rename columns to describe contents
Color.means$Metric <- "Avg Red Channel Intensity - reanalysis"

#or use summarySE, this way you have confidence interval to plot later
Color.means <- summarySE(newdata, measurevar="AvgRedChannel_reanalysis", groupvars=c("colorstate","Origin","PhotoType"))
names(Color.means)[5] <- "mean"


# now plot the average red channel intensity #
Fig.all <-  ggplot(Color.means, aes(x=PhotoType, y=mean,  group=colorstate)) + #set up plot information
  geom_errorbar(aes(x=PhotoType, ymax= mean+ci, ymin= mean-ci), colour="black", width=.1, position = position_dodge(width = 0.2)) + #add standard error bars about the mean
  geom_point(aes(shape=colorstate), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(17,2,16,1)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab('Origin') + #label x axis
  ylab(expression(atop('Average Intensity of the Red Channel',paste('(0=darkest, 255=brightest)')))) + #label y axis
  ylim(110,210)+
  theme_bw() + #Set the background color
  theme(axis.text = element_text(color = 'black', size=10), #Set the text size and color
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_rect(colour='grey95'),
        legend.background = element_rect(colour='black'),
        legend.position=c(0.1,0.2)) + #set legend location
  ggtitle("Chl Proxy") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 14, 
                                  hjust = 0)) #set title attributes
Fig.all #view plot


### Analyze and plot holobiont data only ###

Color.mean <- aggregate(AvgRedChannel_reanalysis ~ colorstate*Origin, data=holobiontdata, mean)
Color.se <- aggregate(AvgRedChannel_reanalysis ~ colorstate*Origin, data=holobiontdata, std.error)
Color.sd <- aggregate(AvgRedChannel_reanalysis ~ colorstate*Origin, data= skeletondata, sd)
Color.n <- aggregate(AvgRedChannel_reanalysis ~ colorstate*Origin, data=holobiontdata, length)
Color.means <- cbind(Color.mean, Color.se$AvgRedChannel_reanalysis, Color.sd$AvgRedChannel_reanalysis, Color.n$AvgRedChannel_reanalysis) #combine mean and standard error results
colnames(Color.means) <- c("ColorState", "Origin","mean", "se", "sd", "n")   #rename columns to describe contents
Color.means$Metric <- "Avg Red Channel Intensity - reanalysis"

#or use summarySE, this way you have confidence interval to plot later
Color.means <- summarySE(holobiontdata, measurevar="AvgRedChannel_reanalysis", groupvars=c("colorstate","Origin"))
names(Color.means)[4] <- "mean"

Fig.holobiont <-  ggplot(Color.means, aes(x=Origin, y=mean,  group=colorstate)) + #set up plot information
  geom_errorbar(aes(x=Origin, ymax= mean+ci, ymin= mean-ci), colour="black", width=.1, position = position_dodge(width = 0.2)) + #add standard error bars about the mean
  geom_point(aes(shape=colorstate), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(17,2,16,1)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab('Origin') + #label x axis
  ylab(expression(atop('Average Intensity of the Red Channel',paste('(0=darkest, 255=brightest)')))) + #label y axis
  ylim(110,230)+
  theme_bw() + #Set the background color
  theme(axis.text = element_text(color = 'black', size=10), #Set the text size and color
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_rect(colour='grey95'),
        legend.background = element_rect(colour='black'),
        legend.position=c(0.12,0.85)) + #set legend location
  ggtitle("Holobiont Chl Proxy") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 14, 
                                  hjust = 0)) #set title attributes
Fig.holobiont #view plot
ggsave(file="~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/RGBAnalysis/RedChannelIntensity_holobiont.CI.pdf", Fig.holobiont, width = 6, height = 6, units = c("in"), useDingbats=FALSE)


### Analyze and plot the skeleton data only ###

Color.mean <- aggregate(AvgRedChannel_reanalysis ~ colorstate*Origin, data=skeletondata, mean)
Color.se <- aggregate(AvgRedChannel_reanalysis ~ colorstate*Origin, data= skeletondata, std.error)
Color.sd <- aggregate(AvgRedChannel_reanalysis ~ colorstate*Origin, data= skeletondata, sd)
Color.n <- aggregate(AvgRedChannel_reanalysis ~ colorstate*Origin, data= skeletondata, length)
Color.means <- cbind(Color.mean, Color.se$AvgRedChannel, Color.sd$AvgRedChannel, Color.n$AvgRedChannel) #combine mean and standard error results
colnames(Color.means) <- c("ColorState", "Origin","mean", "se", "sd", "n")   #rename columns to describe contents
Color.means$Metric <- "Avg Red Channel Intensity"

#or use summarySE, this way you have confidence interval to plot later
Color.means <- summarySE(skeletondata, measurevar="AvgRedChannel_reanalysis", groupvars=c("colorstate","Origin"))
names(Color.means)[4] <- "mean"

Fig.skeleton <-  ggplot(Color.means, aes(x=Origin, y=mean,  group=colorstate)) + #set up plot information
  geom_errorbar(aes(x=Origin, ymax= mean+ci, ymin= mean-ci), colour="black", width=.1, position = position_dodge(width = 0.2)) + #add standard error bars about the mean
  geom_point(aes(shape=colorstate), position = position_dodge(width = 0.2), size=4, lwd=10) + #plot points
  scale_shape_manual(values=c(17,2,16,1)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab('Origin') + #label x axis
  ylab(expression(atop('Average Intensity of the Red Channel',paste('(0=darkest, 255=brightest)')))) + #label y axis
  ylim(110,230)+
  theme_bw() + #Set the background color
  theme(axis.text = element_text(color = 'black', size=10), #Set the text size and color
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_rect(colour='grey95'),
        legend.background = element_rect(colour='black'),
        legend.position=c(.12, .85)) + #set legend location
  ggtitle("Skeleton Chl Proxy") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 14, 
                                  hjust = 0)) #set title attributes
Fig.skeleton #view plot

ggsave(file="~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/RGBAnalysis/RedChannelIntensity_skeleton.CI.pdf", Fig.skeleton, width = 6, height = 6, units = c("in"), useDingbats=FALSE)


Figs.RGB <- arrangeGrob(Fig.all, Fig.holobiont,Fig.skeleton, ncol=3)
ggsave(file="~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/RGBAnalysis/Figs.RGB.pdf", Figs.RGB, width = 11, height = 6, units = c("in"))


### Analyze and plot the before/after ramp comparison ###

Color.mean <- aggregate(AvgRedChannel_reanalysis ~ colorstate*Origin*PhotoTime, data=holobiontdata, mean)
Color.se <- aggregate(AvgRedChannel_reanalysis ~ colorstate*Origin*PhotoTime, data=holobiontdata, std.error)
Color.n <- aggregate(AvgRedChannel_reanalysis ~ colorstate*Origin*PhotoTime, data=holobiontdata, length)

Color.means <- cbind(Color.mean, Color.se$AvgRedChannel_reanalysis, Color.n$AvgRedChannel_reanalysis) #combine mean and standard error results
colnames(Color.means) <- c("ColorState", "Origin","PhotoTime","mean", "se", "n")   #rename columns to describe contents
Color.means$Metric <- "Avg Red Channel Intensity - reanalysis"

Color.means$PhotoTime <- factor(Color.means$PhotoTime, levels = c("BeforeRamp","AfterRamp"))

Fig.beforeafter <-  ggplot(Color.means, aes(x=PhotoTime, y=mean,  group=ColorState)) + #set up plot information
  geom_errorbar(aes(x=PhotoTime, ymax= mean+se, ymin= mean-se), colour="black", width=.1, position = position_dodge(width = 0.2)) + #add standard error bars about the mean
  geom_point(aes(shape=ColorState), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(17,2,16,1)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab('Photo Time') + #label x axis
  ylab(expression(atop('Average Intensity of the Red Channel',paste('(0=darkest, 255=brightest)')))) + #label y axis
  ylim(110,210)+
  theme_bw() + #Set the background color
  theme(axis.text = element_text(color = 'black', size=10), #Set the text size and color
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_rect(colour='grey95'),
        legend.background = element_rect(colour='black'),
        legend.position=c(.2, .2)) + #set legend location
  ggtitle("Before/After Experiment Comparison") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 14, 
                                  hjust = 0)) #set title attributes
Fig.beforeafter #view plot
