# This script takes physiology data (R, P, coral color and Fv/Fm) from VA and RI Astrangia poculata and constructs genetic variance-covariance matrices (G matrices) separately for each population. Then, the Random Skewers method is used to test differences between the two G matrices. Finally, since the VA and RI G matrices are not significantly different from one another, all physiology data is combined to boost confidence in the genetic architecture of the physiology traits considered here. 
# Code originally written by Rachel Wright, updated by Hannah Aichelman 
# Last updated January 2019 


setwd("~/Documents/ODU_MS/ExperimentalTanks/ExperimentData/MCMCglmmAnalysis")

library(tidyverse) # for plotting (ggplot2) and data wrangling (tidyr + dplyr)
library(MCMCglmm) # for MCMCglmm
library(pheatmap) # for pretty heatmaps
library(evolvability) # for evolvability metrics
library(reshape2) # for melt
library(Rmisc) # for summarySE
library(evolqg) # evolutionary quantitative genetics, for comparing G matrices

######## Read in physiology data and compile into single .csv file ("AllData.csv")
#read in final versions of all data (coral color, heat respirometry, cold respirometry, and Fv/Fm)
# colordata <- read.csv("~/Documents/ODU_MS/ExperimentalTanks/ExperimentData/RGBAnalysis/RedChannelIntensity_data.csv")
# fvfmdata <- read.csv("~/Documents/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/PAM/FvFm_RandSample_holo.csv")
# corrRdata <- read.csv("~/Documents/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/Respirometry_ANOVA/CorrR_RandSample.csv")
# holoPdata <- read.csv("~/Documents/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/Respirometry_ANOVA/HoloP_RandSample.csv")

#use tidyverse and dplyr to filter out the data we are interested in here for each of these traits
#write out .csv files for each cleaned dataset to compile in excel 
# colordatafilter <- colordata %>%
	# dplyr::select(FragmentId, Origin, Color, AvgRedChannel_reanalysis, PhotoTime) %>%
	# filter(PhotoTime=="BeforeRamp")
# summary(colordatafilter)
# write.csv(colordatafilter, file = "colordatafilter.csv")

# corrRdatafilter <- corrRdata %>%
	# dplyr::select(Genotype_holo, TempC_holo, Origin_holo, Color_holo, abs_Rcorr) %>%
	# filter(TempC_holo==18 | TempC_holo==22)
# summary(corrRdatafilter)
# write.csv(corrRdatafilter, file = "corrRdatafilter.csv")

# holoPdatafilter <- holoPdata %>%
	# dplyr::select(Genotype_holo, TempC_holo, Origin_holo, Color_holo, Pgross_holo) %>%
	# filter(TempC_holo==18 | TempC_holo==22)
# summary(holoPdatafilter)
# write.csv(holoPdatafilter, file = "holoPdatafilter.csv")

# fvfmdatafilter <- fvfmdata %>%
	# dplyr::select(ID, genotype, color, origin, FvFm, temp) %>%
	# filter(temp==18 | temp==22)
# summary(fvfmdatafilter)
# write.csv(fvfmdatafilter, file = "fvfmdatafilter.csv")


#after compiling, read back in all data:
data <- read.csv("AllData.csv")
summary(data)

names(data)
#put all traits on the same scale
dat_scaled <- data %>% mutate(
    R.scaled = scale(abs_Rcorr, center = T, scale = T), 
    P.scaled = scale(Pgross_holo, center = T, scale = T), 
    FvFm.scaled = scale(FvFm, center = T, scale = T), 
    RedChannel.scaled = scale(AvgRedChannel_reanalysis, center = T, scale = T)) %>%
  dplyr::select(Genotype_holo, Origin_holo, TempC_holo, Color_holo, R.scaled, P.scaled, FvFm.scaled, RedChannel.scaled)
names(dat_scaled)

#Separate by VA and RI populations
VA_data <- subset(dat_scaled, Origin_holo=="VA")
RI_data <- subset(dat_scaled, Origin_holo=="RI")

######## Build VA G Matrix
ntraits <- 4
ntraits

prior <- list(R = list(V = diag(ntraits), nu = ntraits-0.98),
              G = list(G1 = list(V = diag(ntraits), nu = ntraits-0.98, alpha.mu = rep(0, ntraits), alpha.V=diag(ntraits)*1000)))

set.seed(1)
VA.model <- MCMCglmm(cbind(P.scaled, R.scaled, FvFm.scaled, RedChannel.scaled)~TempC_holo,
              random=~us(trait): Genotype_holo,
              rcov=~us(trait):units,
              data=VA_data,
              prior=prior,
              family=rep("gaussian",4),
              nitt=25000,thin=20,burnin=10000)
summary(VA.model)
            # post.mean l-95% CI u-95% CI eff.samp pMCMC
# (Intercept)  -0.35559 -2.02849  1.09082      750 0.640
# TempC_holo    0.01097 -0.06404  0.08646      750 0.741

plot(VA.model)

summary(VA.model)$Gcovariances #G-matrix entries with 95% credible intervals
# "Significant" if credible interval doesn't include zero
signif <- summary(VA.model)$Gcovariances[,2]*summary(VA.model)$Gcovariances[,3]>0   
table(signif)
# FALSE  TRUE 
   # 12     4 
summary(VA.model)$Gcovariances[signif=="TRUE",]
                                                            # post.mean     l-95% CI  u-95% CI eff.samp
# traitP.scaled:traitP.scaled.Genotype_holo                   0.2967578 2.296246e-04 0.6839144 750.0000
# traitR.scaled:traitR.scaled.Genotype_holo                   0.1023639 1.369075e-07 0.3365503 750.0000
# traitFvFm.scaled:traitFvFm.scaled.Genotype_holo             0.3241452 1.876100e-05 0.8858626 750.0000
# traitRedChannel.scaled:traitRedChannel.scaled.Genotype_holo 1.3917533 5.469316e-01 2.7197390 592.3525

# Construct G-Matrix
gmatVA <- matrix(summary(VA.model)$Gcovariances,nrow=4,ncol=4)
dimnames(gmatVA) <- list(c("Pgross", "Rcorr", "FvFm", "Color"), c("Pgross", "Rcorr", "FvFm", "Color"))
round(gmatVA,4) 
                         # Pgross_holo abs_Rcorr    FvFm AvgRedChannel_reanalysis
# Pgross_holo                   0.2968   -0.0585  0.1753                  -0.4532
# abs_Rcorr                    -0.0585    0.1024 -0.0890                   0.1615
# FvFm                          0.1753   -0.0890  0.3241                  -0.4224
# AvgRedChannel_reanalysis     -0.4532    0.1615 -0.4224                   1.3918

eigen(gmatVA)$values # all positive eigevalues
# [1] 1.72849524 0.18912778 0.12716873 0.07022846

#make heatmap
colf <- colorRampPalette(c("lightyellow","gold","coral","red"))
colf2 <- colorRampPalette(c("skyblue","lightyellow"))
cols <- c(colf2(10),colf(95))
cols <- colf(100)
pheatmap(gmatVA,color=cols,treeheight_col=4,treeheight_row=4, main="VA G Matrix")

######## Build RI G Matrix
ntraits <- 4
ntraits

prior <- list(R = list(V = diag(ntraits), nu = ntraits-0.98),
              G = list(G1 = list(V = diag(ntraits), nu = ntraits-0.98, alpha.mu = rep(0, ntraits), alpha.V=diag(ntraits)*1000)))

set.seed(1)
RI.model <- MCMCglmm(cbind(P.scaled, R.scaled, FvFm.scaled, RedChannel.scaled)~TempC_holo,
              random=~us(trait): Genotype_holo,
              rcov=~us(trait):units,
              data=RI_data,
              prior=prior,
              family=rep("gaussian",4),
              nitt=25000,thin=20,burnin=10000)
summary(RI.model)
            # post.mean  l-95% CI  u-95% CI eff.samp pMCMC
# (Intercept)  0.147732 -2.058096  2.470727      750 0.899
# TempC_holo  -0.005392 -0.127685  0.100817      750 0.939

plot(RI.model)

summary(RI.model)$Gcovariances #G-matrix entries with 95% credible intervals
# "Significant" if credible interval doesn't include zero
signif <- summary(RI.model)$Gcovariances[,2]*summary(RI.model)$Gcovariances[,3]>0   
table(signif)
# FALSE  TRUE 
   # 12     4 
summary(RI.model)$Gcovariances[signif=="TRUE",]
                                                            # post.mean     l-95% CI  u-95% CI eff.samp
# traitP.scaled:traitP.scaled.Genotype_holo                   0.6029728 1.716538e-04 1.5972196 750.0000
# traitR.scaled:traitR.scaled.Genotype_holo                   0.3373048 3.637771e-06 0.9984210 750.0000
# traitFvFm.scaled:traitFvFm.scaled.Genotype_holo             0.2004272 1.748888e-08 0.5619874 750.0000
# traitRedChannel.scaled:traitRedChannel.scaled.Genotype_holo 0.8964751 2.449049e-01 1.6819281 612.6676

# Construct G-Matrix
gmatRI <- matrix(summary(RI.model)$Gcovariances,nrow=4,ncol=4)
dimnames(gmatRI) <- list(c("Pgross", "Rcorr", "FvFm", "Color"), c("Pgross", "Rcorr", "FvFm", "Color"))
round(gmatRI,4) 
                         # Pgross_holo abs_Rcorr    FvFm AvgRedChannel_reanalysis
# Pgross_holo                   0.6030   -0.0349  0.1898                  -0.4596
# abs_Rcorr                    -0.0349    0.3373 -0.0730                   0.0707
# FvFm                          0.1898   -0.0730  0.2004                  -0.2275
# AvgRedChannel_reanalysis     -0.4596    0.0707 -0.2275                   0.8965

eigen(gmatRI)$values # all positive eigevalues
# [1] 1.3198965 0.3411749 0.2688150 0.1072934

# make heatmap
colf <- colorRampPalette(c("lightyellow","gold","coral","red"))
colf2 <- colorRampPalette(c("skyblue","lightyellow"))
cols <- c(colf2(10),colf(95))
cols <- colf(100)
pheatmap(gmatRI,color=cols,treeheight_col=4,treeheight_row=4)
#pheatmap(gmatRI, main = "RI G Matrix")

### Compare VA and RI G matrices

# this test utilizes the Lande equation (delta_z = GB). The mean correlation resulting from this test is how correlated the responses of the two populations are to a selection pressure (Melo et al. 2015 "EvolQG - An R package for evolutionary quantitative genetics"). In other words, this is a statistic of how often two populations respond similarly (in the same direction) to the same selctive pressure (B). Below, our two populations respond significantly the same in response to the same selective pressure.
# values range between -1 (matrices have opposite structure) to 1 (matrices share the same structure), 0 means the matrices have distinct structure. 
# Steppan et al 2002 (TREE), Aguirre et al 2014 (Heredity), and Mezey et al 2003 are other papers to consider when thinking about G matrices and comparing them across populations
RandomSkewers(gmatRI, gmatVA)
   # correlation    probability correlation_sd 
     # 0.8772386      0.0300000      0.1600074 


### Since G-matrices are more or less the same, all corals from both locations in a single G-matrix, fit multi-trait model with fixed effects of temperature*location and random effect of genet

ntraits <- 4
ntraits

prior <- list(R = list(V = diag(ntraits), nu = ntraits-0.98),
              G = list(G1 = list(V = diag(ntraits), nu = ntraits-0.98, alpha.mu = rep(0, ntraits), alpha.V=diag(ntraits)*1000)))

set.seed(1)
AllSamples.model <- MCMCglmm(cbind(P.scaled, R.scaled, FvFm.scaled, RedChannel.scaled)~TempC_holo*Origin_holo,
              random=~us(trait): Genotype_holo,
              rcov=~us(trait):units,
              data=dat_scaled,
              prior=prior,
              family=rep("gaussian",4),
              nitt=25000,thin=20,burnin=10000)
summary(AllSamples.model)
                         # post.mean l-95% CI u-95% CI eff.samp pMCMC
# (Intercept)               -0.44227 -1.87981  1.01494    675.3 0.552
# TempC_holo                 0.02573 -0.05017  0.09597    673.9 0.517
# Origin_holoVA              0.59984 -0.72185  1.88342    750.0 0.384
# TempC_holo:Origin_holoVA  -0.03706 -0.10060  0.03033    750.0 0.272

plot(AllSamples.model)

summary(AllSamples.model)$Gcovariances #G-matrix entries with 95% credible intervals
# "Significant" if credible interval doesn't include zero
signif <- summary(AllSamples.model)$Gcovariances[,2]*summary(AllSamples.model)$Gcovariances[,3]>0   
table(signif)
# FALSE  TRUE 
   # 6     10 
summary(AllSamples.model)$Gcovariances[signif=="TRUE",]
                                                             # post.mean      l-95% CI   u-95% CI eff.samp
# traitP.scaled:traitP.scaled.Genotype_holo                    0.3698456  5.937816e-02  0.7247436 818.6552
# traitFvFm.scaled:traitP.scaled.Genotype_holo                 0.2405283  4.732664e-02  0.4736368 750.0000
# traitRedChannel.scaled:traitP.scaled.Genotype_holo          -0.4910866 -9.086046e-01 -0.2030878 750.0000
# traitR.scaled:traitR.scaled.Genotype_holo                    0.1915306  6.335007e-08  0.4996225 750.0000
# traitP.scaled:traitFvFm.scaled.Genotype_holo                 0.2405283  4.732664e-02  0.4736368 750.0000
# traitFvFm.scaled:traitFvFm.scaled.Genotype_holo              0.2749240  1.542310e-04  0.5992919 750.0000
# traitRedChannel.scaled:titFvFm.scaled.Genotype_holo       -0.4209912 -8.232711e-01 -0.1126043 750.0000
# traitP.scaled:traitRedChannel.scaled.Genotype_holo          -0.4910866 -9.086046e-01 -0.2030878 750.0000
# traitFvFm.scaled:traitRedChannel.scaled.Genotype_holo       -0.4209912 -8.232711e-01 -0.1126043 750.0000
# traitRedChannel.scaled:traitRedChannel.scaled.Genotype_holo  1.0722500  5.145158e-01  1.6740529 750.0000

# Construct G-Matrix
gmatALL <- matrix(summary(AllSamples.model)$Gcovariances,nrow=4,ncol=4)
dimnames(gmatALL) <- list(c("Pgross", "Rcorr", "FvFm", "Color"), c("Pgross", "Rcorr", "FvFm", "Color"))

round(gmatALL,4) 
        # Pgross   Rcorr    FvFm   Color
# Pgross  0.3698 -0.0920  0.2405 -0.4911
# Rcorr  -0.0920  0.1915 -0.1112  0.1657
# FvFm    0.2405 -0.1112  0.2749 -0.4210
# Color  -0.4911  0.1657 -0.4210  1.0723

eigen(gmatALL)$values # all positive eigevalues
# [1] 1.54889602 0.17320993 0.12021238 0.06623193

# make heatmap
colf <- colorRampPalette(c("lightyellow","gold","coral","red"))
colf2 <- colorRampPalette(c("skyblue","lightyellow"))
cols <- c(colf2(10),colf(95))
cols <- colf(100)
GMatAll <- pheatmap(gmatALL,color=cols,treeheight_col=4,treeheight_row=4)
ggsave(file="GMatAllPops.pdf", GMatAll, width = 4, height = 4, units = c("in"), useDingbats=FALSE)
#pheatmap(gmatALL, main = "All Samples G Matrix")

