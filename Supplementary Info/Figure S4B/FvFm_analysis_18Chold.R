#Last updated May 18, 2018 by Hannah Aichelman
# This script analyzes Fv/Fm data from the 18C hold experiment conducted only on VA Astrangia on December 20, 2016. Using this data to try and understand if Fv/Fm is affected by the amount of time the corals are kept in the chambers. This experiment was run very similarly to the later official TPC experiments, but temp was maintained at 18C throughout

#library packages needed
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
library("Rmisc") # has the summarySE function

#load in data from experiment
setwd('~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/20161220_18holdramp/20161220_pam')
data <- read.csv('20161220_PAM_cleaned.csv', header=TRUE, sep=",")

head(data)
###################################################################

data <-cbind(data[0:5],"colorstate"=paste(data$origin, data$color, sep="_"), data[6]) 

data_agg<-aggregate(FvFm ~ time + temp + color + origin + genotype + colorstate, data=data, mean) ###Average replicate measurements for each individual by temp and day

data_agg$time<-factor(data_agg$time)

#check replication of genotype x temp data for holobiont
replications(FvFm ~ genotype*time, data_agg)

aggregate(FvFm ~ color + origin + time, data=data_agg, FUN=function(x) shapiro.test(x)$p.value) #shapiro test assesses normality of data. 
boxplot(FvFm ~ colorstate + time, data=data_agg, las=2)

#Bartlett test for homogeneity of variance
for(i in names(table(data_agg$time))){
	print(paste0("Time",i))
	print(bartlett.test(FvFm ~ colorstate, data=data_agg[data_agg$temp==i,]))
}

brown_data <- subset(data_agg, data_agg$color=="brown")

#ANOVA
aov<-aov(FvFm ~ time*color + Error(genotype/time), data=data_agg) 
summary(aov)

# > summary(aov)

# Error: genotype
          # Df Sum Sq Mean Sq F value  Pr(>F)   
# color      1 0.5957  0.5957   61.71 0.00142 **
# Residuals  4 0.0386  0.0097                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Error: genotype:time
           # Df  Sum Sq  Mean Sq F value  Pr(>F)   
# time        5 0.01665 0.003330    3.15 0.02943 * 
# time:color  5 0.02294 0.004587    4.34 0.00775 **
# Residuals  20 0.02114 0.001057                   
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

print(lsmeans(aov, list(pairwise ~ color)), adjust = c("tukey"))
# > print(lsmeans(aov, list(pairwise ~ color)), adjust = c("tukey"))
# NOTE: Results may be misleading due to involvement in interactions
# $`lsmeans of color`
 # color     lsmean         SE    df    lower.CL  upper.CL
 # brown 0.34354630 0.02370125 14.16  0.28430164 0.4027910
 # white 0.05110185 0.04042678  5.89 -0.04995069 0.1521544

# Results are averaged over the levels of: time 
# Confidence level used: 0.95 
# Conf-level adjustment: sidak method for 2 estimates 

# $`pairwise differences of contrast`
 # contrast       estimate         SE   df t.ratio p.value
 # brown - white 0.2924444 0.04074047 9.04   7.178  0.0001

# Results are averaged over the levels of: time 

print(lsmeans(aov, list(pairwise ~ time)), adjust = c("tukey"))

# > print(lsmeans(aov, list(pairwise ~ time)), adjust = c("tukey"))
# NOTE: Results may be misleading due to involvement in interactions
# $`lsmeans of time`
 # time    lsmean         SE    df   lower.CL  upper.CL
 # 1    0.2007593 0.02315789  4.00 0.08909870 0.3124198
 # 2    0.2195370 0.02980963 10.11 0.07580374 0.3632703
 # 3    0.1675926 0.02980963 10.11 0.02385930 0.3113259
 # 10   0.2104259 0.02980963 10.11 0.06669263 0.3541592
 # 11   0.1685926 0.02980963 10.11 0.02485930 0.3123259
 # 12   0.2170370 0.02980963 10.11 0.07330374 0.3607703

# Results are averaged over the levels of: color 
# Confidence level used: 0.95 
# Conf-level adjustment: sidak method for 6 estimates 

# $`pairwise differences of contrast`
 # contrast     estimate         SE df t.ratio p.value
 # 1 - 2    -0.018777778 0.01877035 20  -1.000  0.9124
 # 1 - 3     0.033166667 0.01877035 20   1.767  0.5071
 # 1 - 10   -0.009666667 0.01877035 20  -0.515  0.9950
 # 1 - 11    0.032166667 0.01877035 20   1.714  0.5389
 # 1 - 12   -0.016277778 0.01877035 20  -0.867  0.9501
 # 2 - 3     0.051944444 0.01877035 20   2.767  0.1050
 # 2 - 10    0.009111111 0.01877035 20   0.485  0.9962
 # 2 - 11    0.050944444 0.01877035 20   2.714  0.1161
 # 2 - 12    0.002500000 0.01877035 20   0.133  1.0000
 # 3 - 10   -0.042833333 0.01877035 20  -2.282  0.2467
 # 3 - 11   -0.001000000 0.01877035 20  -0.053  1.0000
 # 3 - 12   -0.049444444 0.01877035 20  -2.634  0.1345
 # 10 - 11   0.041833333 0.01877035 20   2.229  0.2685
 # 10 - 12  -0.006611111 0.01877035 20  -0.352  0.9992
 # 11 - 12  -0.048444444 0.01877035 20  -2.581  0.1482


print(lsmeans(aov, list(pairwise ~ color | time)), adjust = c("tukey"))

# use summarySE to summarize data, this way you have confidence interval to plot later
all.means <- summarySE(data_agg, measurevar="FvFm", groupvars=c("colorstate","time"))
names(all.means)[4] <- "mean"


#plot all colorstates data using ggplot
all.means$time <- factor(all.means$time, levels = c("10", "11", "12", "1", "2", "3"))

Fig.FvFm <-  ggplot(all.means, aes(x=time, y=mean,  group=colorstate)) + #set up plot information
  geom_errorbar(aes(x=time, ymax=mean+ci, ymin=mean-ci), colour="black", width=.1, position = position_dodge(width = 0.6)) + #add standard error bars about the mean
  geom_point(aes(shape=colorstate), position = position_dodge(width = 0.6), size=5) + #plot points
  scale_shape_manual(values=c(16,1)) + #sets point shape manually
  #geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Time") + #label x axis
    #scale_x_discrete(limits=c(10,11,12,1,2,3)) +
  ylab("Fv/Fm") + #label y axis
  ylim(-0.1, 0.7)+
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(color = 'black', size=10), #Set the text angle
        axis.text.y = element_text(color = 'black', size=10),
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.position=c(.08,.1)) + #set legend location
  #ggtitle("Fv/Fm - 18C hold") + #add a main title
  #theme(plot.title = element_text(face = 'bold', 
                                  size = 14, 
                                  hjust = 0)) #set title attributes
Fig.FvFm #view plot

ggsave(file="~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/20161220_18holdramp/20161220_pam/FvFm_CI.pdf", Fig.FvFm, width = 6, height = 6, units = c("in"), useDingbats=FALSE)


# Or make the figure using base R plotting script:
pdf("DansRCorrPlot.pdf",9,7)
par(mar=c(4.5,4.5,1,1))
x.pos<-as.factor(levels(all.means$time)[all.means$time])+c(-0.5)
x.pos<-factor(all.means$time, levels = c("10", "11", "12", "1", "2", "3"))
sym=16
indv=1:6 #VA brown
errbar(x.pos[indv], all.means$mean[indv], all.means$mean[indv]+ all.means$ci[indv], all.means$mean[indv]-all.means$ci[indv],  ylim=c(-.5, 1), pch=sym, type="p", cex=2.5, ylab=bquote('Fv/Fm'), xlab="Time" , xaxt='n') #type = b to connect dots

sym=1
indv=7:12 #VA white
errbar(x.pos[indv], all.means$mean[indv], all.means$mean[indv]+ all.means$ci[indv], all.means$mean[indv]-all.means$ci[indv],  pch=sym, type="p", cex=2.5,add=TRUE, xaxt='n')
#points(predsVAwhite$TempC, predsVAwhite$unlogfitted, type="l", lty=3, col="red", lwd=2)
#points(predsVAwhite$TempC, predsVAwhite$unloglwrCI, type="l", lty=3, col="blue", lwd=2)
#points(predsVAwhite$TempC, predsVAwhite$unloguprCI, type="l", lty=3, col="blue", lwd=2)

axis(1,at=c(6,9,12,15,18,22,26,29,32), labels=c(6,9,12,15,18,22,26,29,32))
legend("topleft",c("RI_brown","RI_white","VA_brown","VA_white"), pch=c(17,2,16,1), pt.cex=1.5, bty="n")
#legend(9,2.6,c("RI_brown fit","RI_white fit","VA_brown fit","VA_white fit"), lty=c(1,3,1,3), col=c("blue","blue","red","red"),lwd=2, bty="n")
text(5,seq(1.2,0.8,by=-0.3/2),c("*** Color, white > brown", "*** Temp", "** Color*Temp"),adj=c(0,0))
text(c(6,12,15,18,22,26,29,32),-.52,'+',cex=2)
dev.off()

