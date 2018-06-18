#Author: Hannah Aichelman
#Last updated: April 2, 2018
# Analysis of Topt estimates obtained using updated Padfield approach

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

# read in data
#corrected R
setwd("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/TPCFitting_Padfield_Analysis/DataOutput/Corr_R")
data <- read.csv("ParameterEstimateCompiled_CorrR.csv")

# Topt from population averages, both R and P:
setwd("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/TPCFitting_Padfield_Analysis/DataOutput/")
data <- read.csv("PopAv_Topt_Estimates_Compiled.csv")

head(data)


# change some variables to fators
data$term<-factor(data$term)
data$genotype<-factor(data$genotype)
data$color<-factor(data$color)
data$colorstate<-factor(data$colorstate)


### Analysis used when doing Topt estimates for every individual:
data$Topt_C <- data$Topt - 273.15
data <-cbind(data[0:8],"colorstate"=paste(data$origin, data$color, sep="_"), data[9]) 

Topt <- subset(data, term == 'mu')

#look at estimates for colorstaates separately
Toptvab <- subset(Topt, colorstate == 'va_brown')
Toptvaw <- subset(Topt, colorstate == 'va_white')
Toptrib <- subset(Topt, colorstate == 'ri_brown')
Toptriw <- subset(Topt, colorstate == 'ri_white')

#plot all estimates to see how well the bootstrapping approach is estimating the Topt
ggplot(Topt, aes(x=colorstate, y=Topt_C, color=colorstate)) +
	geom_boxplot(outlier.color = "black", outlier.shape = 16, outlier.size = 2)+
	theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank())
	labs(title = "Topt Estimates")

#scatterplot	
ggplot(data, aes(x=colorstate, y=Topt_C, color=colorstate)) +
	geom_point(outlier.color = "black", outlier.shape = 16, outlier.size = 2)+
	labs(title = "Topt Estimates")


#get rid of negative estimates, since they are unreasonable estimates from the bootstrapping approach
Topt <- subset(Topt, Topt_C > 0)

# boxplot of Topt based on colorstate
# this is the most recent version of the plot in an attempt to match it with the respiration temperature axis
Fig.boxplot <- ggplot(Topt, aes(x=colorstate, y=Topt_C)) + 
	geom_boxplot()+
	ylab("Temperature (°C)")+
	expand_limits(y=c(6,32))+
	scale_y_continuous(breaks=c(6,9,12,15,18,22,26,29,32)) +
	xlab("Population")+	
	theme_bw() + #Set the background color
  	theme(axis.text.x = element_text(color = 'black', size=10), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_rect(colour='grey95'),
        legend.background = element_rect(colour='black'),
        legend.position=c(.09,.85)) #set legend location	  
Fig.boxplot
# save the plot so the dimensions of the pdf files of respiration and this boxplot are the same
ggsave(file="~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/TPCFitting_Padfield_Analysis/DataOutput/Fig.boxplot.pdf", Fig.boxplot, width = 6, height = 8, units = c("in"), useDingbats=FALSE)

#and another way to plot it
boxplot(Topt_C ~ colorstate, data=Topt, las=1, xlab = "Temp", ylab = "population", ylim=c(6,32), yaxt='n')
axis(2,at=c(6,9,12,15,18,22,26,29,32)) #designate specific axes so we can line up exactly with the temperature axis of metabolic rate figures
#save plot
quartz.save("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/TPCFitting_Padfield_Analysis/DataOutput/Topt_R_Uncorr_Trim_forPPT.pdf", type="pdf", width = 6, height = 8, device=dev.cur())

#Bartlett test for homogeneity of variance
bartlett.test(Topt_C ~ colorstate, data=Topt)

#Anova
Topt.aov<-aov(Topt ~ origin*color, data=Topt)
summary(Topt.aov)
print(lsmeans(Topt.aov, list(pairwise ~ origin)), adjust = c("tukey"))


##############
Topt.mean <- aggregate(Topt_C ~ colorstate, data= Topt, mean)
Topt.se <- aggregate(Topt_C ~ colorstate, data= Topt, std.error)
Topt.sd <- aggregate(Topt_C ~ colorstate, data= Topt, sd)
Topt.n <- aggregate(Topt_C ~ colorstate, data= Topt, length)
Topt.means <- cbind(Topt.mean, Topt.se$Topt_C, Topt.sd$Topt_C, Topt.n$Topt_C)
colnames(Topt.means) <- c("colorstate", "mean", "se","sd","n")
Topt.means $Metric <- "Topt"

#plot all colorstates data
Fig.Topt <-  ggplot(Topt.means, aes(x=colorstate, y=mean,  group=colorstate)) + #set up plot information
  geom_errorbar(aes(x=colorstate, ymax=mean+sd, ymin=mean-sd), colour="black", width=.1, position = position_dodge(width = 0.2)) + #add standard error bars about the mean
  geom_point(aes(shape=colorstate), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(17,2,16,1)) + #sets point shape manually
  geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Colorstate") + #label x axis
    #scale_x_discrete(breaks=c(6,9,12,15,18,22,26,29,32)) +
  ylab("Topt (°C)") + #label y axis
  ylim(15, 35)+
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(color = 'black', size=10), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_rect(colour='grey95'),
        legend.background = element_rect(colour='black'),
        legend.position=c(.8,.2)) + #set legend location
  ggtitle("Optimum Temperature") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 14, 
                                  hjust = 0)) #set title attributes
Fig.Topt #view plot

### Analysis used when doing Topt estimates by population
# Convert K back to C
data$Topt_C <- data$Topt - 273.15
data$Topt_lwrCI_C <- data$Topt_lwrCI - 273.15
data$Topt_upr_CI_C <- data$Topt_uprCI - 273.15

head(data)

#scatterplot	
ggplot(data, aes(x=colorstate, y=Topt_C, color=colorstate)) +
	geom_point(outlier.color = "black", outlier.shape = 16, outlier.size = 2)+
	labs(title = "Topt Estimates")

ggplot(data, aes(x=colorstate, y=Topt_C, color=rate)) +
	#geom_errorbar(aes(x=colorstate, ymax=upr_CI_C, ymin=lwr_CI_C))+
	geom_point(size=3)+
	ylim(24,30)+
	theme_bw()+
	theme(axis.text.x = element_text(color = 'black', size=10),
		axis.text.y = element_text(color = 'black', size=10),
		axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.position=c(.08,.85)) + #set legend location
	labs(title = "Topt Estimates")


