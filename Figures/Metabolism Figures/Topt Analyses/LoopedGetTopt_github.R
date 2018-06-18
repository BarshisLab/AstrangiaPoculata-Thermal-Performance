#Written by Daniel Padfield, updated by Hannah Aichelman to analyze Astrangia data
#Script adapted from D Padfield's github here: https://padpadpadpad.github.io/post/bootstrapping-non-linear-regressions-with-purrr/
#Last updated March 3, 2018

#devtools::install_github('padpadpadpad/nls.multstart')
#devtools::install_github('thomasp85/patchwork')

# load packages
library(nls.multstart) 
library(patchwork) 
library(ggplot2)
library(broom)
library(purrr)
library(dplyr)
library(tidyr)
library(nlstools)
library(modelr)
library(readr)

setwd("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/TPCFitting_Padfield_Analysis")

#Using same randomly sampled Corrected Respiration data used for plotting
holo.data<-read.csv("~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respirometry_Experiments/Respirometry_ANOVA/CorrR_RandSample.csv")

glimpse(holo.data)

#Visualize all ln_R data in a couple different ways:
plot(holo.data$TempC_holo, holo.data$ln_R)

Fig.all.R <-  ggplot(holo.data, aes(x=TempC_holo,y=ln_R,color=colorstate)) + #set up plot information
  #geom_errorbar(aes(x=temp, ymax=mean+sd, ymin=mean-sd), colour="black", width=.1, position = position_dodge(width = 0.2)) + #add standard error bars about the mean
  geom_point(aes(shape=colorstate), position = position_dodge(width = 0.2), size=4) + #plot points
  scale_shape_manual(values=c(17,2,16,1)) + #sets point shape manually
  #geom_line(aes(), position = position_dodge(width = 0.2), size = 0.5) + #add lines
  xlab("Temperature (°C)") + #label x axis
    #scale_x_discrete(breaks=c(6,9,12,15,18,22,26,29,32)) +
  ylab(bquote('ln_Respiration ('*mu~ 'mol' ~O[2] ~cm^-2~h^-1*')')) + #label y axis
  ylim(-10, 2.5)+
  theme_bw() + #Set the background color
  theme(axis.text.x = element_text(color = 'black', size=10), #Set the text angle
        axis.line = element_line(color = 'black'), #Set the axes color
        panel.border = element_blank(), #Set the border
        panel.grid.major = element_blank(), #Set the major gridlines
        panel.grid.minor = element_blank(), #Set the minor gridlines
        plot.background=element_blank(),  #Set the plot background
        legend.key = element_rect(colour='grey95'),
        legend.background = element_rect(colour='black'),
        legend.position=c(.09,.2)) + #set legend location
  ggtitle("Respiration - Holobiont") + #add a main title
  theme(plot.title = element_text(face = 'bold', 
                                  size = 14, 
                                  hjust = 0)) #set title attributes
Fig.all.R #view plot


# write function for sharpe schoolfield model
schoolfield_high <- function(lnc, E, Eh, Th, temp, Tc) {
  Tc <- 273.15 + Tc
  k <- 8.62e-5
  boltzmann.term <- lnc + log(exp(E/k*(1/Tc - 1/temp)))
  inactivation.term <- log(1/(1 + exp(Eh/k*(1/Th - 1/temp))))
  return(boltzmann.term + inactivation.term)
}


#######
# for ln_R data:
for (i in levels(holo.data$Genotype_holo)){
  d_1 <- subset(holo.data, Genotype_holo == i) 
  
  print(ggplot(d_1, aes(TempK_holo - 273.15, ln_R)) +
    geom_point(col = 'green4') +
    ylab('log Metabolic rate') +
    xlab('Assay temperature (ºC)') +
    theme_bw(base_size = 12, base_family = 'Helvetica'))
  quartz.save(paste0(i,"_LogMetabolicRate.pdf"), type="pdf",device=dev.cur())
  
  # run nls_multstart
  #Tc in degrees C, Th can't be lower than 0C (273.15K)
  fit <- nls_multstart(ln_R ~ schoolfield_high(lnc, E, Eh, Th, temp = TempK_holo, Tc = 18),
                       data = d_1,
                       iter = 1000,
                       start_lower = c(lnc=-10, E=0.1, Eh=0.2,Th=273.15),
                       start_upper = c(lnc=10, E=2, Eh=5, Th=330),
                       supp_errors = 'Y',
                       na.action = na.omit,
                       lower = c(lnc = -10, E = 0, Eh = 0, Th = 0))
  
  fit
  
  fit_boots <- d_1 %>% 
    modelr::bootstrap(., n = 200, id = 'boot_num') %>%
    group_by(boot_num) %>%
    mutate(., fit = map(strap, ~nls_multstart(ln_R ~ schoolfield_high(lnc, E, Eh, Th, temp = TempK_holo, Tc = 18),
                                              data = data.frame(.),
                                              iter = 100,
                                              start_lower = c(lnc = -10, E = 0.1, Eh = 0.2, Th = 273.15),
                                              start_upper = c(lnc = 10, E = 2, Eh = 5, Th = 330),
                                              lower = c(lnc=-10, E=0, Eh=0, Th=0),
                                              na.action = na.omit,
                                              supp_errors = 'Y')
    ))
  
  fit_boots
  
  # calculate confidence intervals of predictions
  new_data <- d_1 %>%
    do(data.frame(TempK_holo = seq(min(.$TempK_holo, na.rm = TRUE), max(.$TempK_holo, na.rm = TRUE), length.out = 200), stringsAsFactors = FALSE))
  
  # make predictions
  preds <- augment(fit, newdata = new_data)
  preds <- fit_boots %>%
    unnest(fit %>% map(augment, newdata = new_data)) %>%
    group_by(TempK_holo) %>%
    summarise(., lwr_CI = quantile(.fitted, 0.025),
              upr_CI = quantile(.fitted, 0.975)) %>%
    ungroup() %>%
    merge(., select(preds, TempK_holo, .fitted), by = 'TempK_holo')
  
  # plot
  print(ggplot() +
    geom_point(aes(TempK_holo - 273.15, ln_R), d_1) +
    geom_line(aes(TempK_holo - 273.15, .fitted), preds) +
    geom_ribbon(aes(TempK_holo - 273.15, ymin = lwr_CI, ymax = upr_CI), alpha = 0.2, preds) +
    theme_bw() +
    ggtitle('Bootstrapping using new datasets'))
  quartz.save(paste0(i,"_BootstrappingOriginal.pdf"), type="pdf",device=dev.cur())
  
  ### Starting here with the alternative bootstrap approach following Thomas (2012) in Science
  # get preds for each point
  preds <- augment(fit) %>%
    select(., - c(X.weights., .resid)) %>%
    rename(., fitted = .fitted)
  
  # get sigma (sd) for each fit
  info <- glance(fit) %>%
    select(sigma) %>%
    # create new columns needed to create adjusted sigma
    # number of samples - no NAs can be present here - will need to be changed if doing multiple curves
    mutate(., n_samps = nrow(d_1),
           # chi_sq number
           chi_sq_num = rchisq(1, n_samps - 1),
           new_sd = sigma^2*sqrt((n_samps -1)/chi_sq_num))
  
  # merge two dataframe together
  boots <- merge(preds, info, all = TRUE)
  glimpse(boots)
  
  # we can now simulate new datasets based on these columns and the above equation
  # define number of bootstraps
  nboot <- 1000
  
  # create new replicate dataframes
  boots <- mutate(boots, n = 1:n()) %>%
    # creates nboot replicates of the dataset
    slice(rep(1:n(), times = nboot)) %>%
    # add column for bootnum
    mutate(., boot_num = rep(1:nboot, each = n()/nboot),
           # new n col
           n = 1:n()) %>%
    # nest the dataset per row
    nest(., -n) %>%
    # create a new simulated y value based on the predicted value and the new sd
    mutate(., sim_y = map_dbl(data, ~ rnorm(n = 1, mean = .x$fitted, sd = .x$new_sd))) %>%
    unnest()
  
  boots <- boots %>%
    group_by(boot_num) %>%
    nest() %>%
    mutate(fit = purrr::map(data, ~ nls_multstart(sim_y ~ schoolfield_high(lnc, E, Eh, Th, temp = TempK_holo, Tc = 18),
                                                  data = data.frame(.),
                                                  iter = 100,
                                                  start_lower = c(lnc = -10, E = 0.1, Eh = 0.2, Th = 273.15),
                                                  start_upper = c(lnc = 10, E = 2, Eh = 5, Th = 330),
                                                  supp_errors = 'Y',
                                                  na.action = na.omit,
                                                  lower = c(lnc = -10, E = 0, Eh = 0, Th = 0))))
  
  # calculate confidence intervals of predictions
  new_data <- d_1 %>%
    do(data.frame(TempK_holo = seq(min(.$TempK_holo, na.rm = TRUE), max(.$TempK_holo, na.rm = TRUE), length.out = 200), stringsAsFactors = FALSE))
  
  # make predictions
  preds <- augment(fit, newdata = new_data)
  preds <- boots %>%
    unnest(fit %>% map(augment, newdata = new_data)) %>%
    group_by(TempK_holo) %>%
    summarise(., lwr_CI = quantile(.fitted, 0.025),
              upr_CI = quantile(.fitted, 0.975)) %>%
    ungroup() %>%
    merge(., select(preds, TempK_holo, .fitted), by = 'TempK_holo')
#write.csv(preds, "Preds_RIbrown_CorrR.csv")  #to write out the preds to take over to the TPC figures
  # plot
  print(ggplot() +
    geom_point(aes(TempK_holo - 273.15, ln_R), d_1) +
    geom_line(aes(TempK_holo - 273.15, .fitted), preds) +
    geom_ribbon(aes(TempK_holo - 273.15, ymin = lwr_CI, ymax = upr_CI), alpha = 0.2, preds) +
    theme_bw() +
    ggtitle('Bootstrapping incorporating uncertainty of sigma'))
  quartz.save(paste0(i,"_BootstrappingSigma.pdf"), type="pdf",device=dev.cur())
  
  # And from these we can calculate confidence intervals for Topt and other parameters
  
  get_topt <- function(E, Th, Eh){
    return((Eh*Th)/(Eh + (8.62e-05 *Th*log((Eh/E) - 1))))
  }
  
  # get parameters, calculate Topt and get confidence intervals
  params <- boots %>%
    unnest(fit %>% map(tidy)) %>%
    select(boot_num, term, estimate) %>%
    spread(term, estimate) %>%
    mutate(Topt = get_topt(E, Th, Eh)) %>%
    gather(., 'term', 'estimate', c(2:ncol(.))) %>%
    group_by(., term) %>%
    summarise(., mu = mean(estimate, na.rm = TRUE),
              lwr_CI = quantile(estimate, 0.025, na.rm = TRUE),
              upr_CI = quantile(estimate, 0.975, na.rm = TRUE)) %>%
    ungroup()
  write_csv(params, path = paste0("/Users/hannahaichelman/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/TPCFitting_Padfield_Analysis/DataOutput/",i,"_ParameterEstimate.csv"))
 
  # And plot these estimates...
  print(ggplot(params, aes(term, y = mu, ymin = lwr_CI, ymax = upr_CI)) +
    geom_point(size = 3) +
    geom_linerange() +
    facet_wrap(~ term, scales = 'free') +
    theme_bw())
  quartz.save(paste0(i,"_ParameterEstimates.pdf"), type="pdf",device=dev.cur())
}

