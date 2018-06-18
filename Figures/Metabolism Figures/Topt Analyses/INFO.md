# Info on Topt calculations and plots

-LoopedGetTopt_github.R reads in the randomly sampled corrected respiration data (CorrR_RandSample.csv), loops through by individual, and calculates Topt for each individual. These estimates were used to create the bar plot of Topt in Figure 3 (using Topt_analysis.R).

-The population-level (RI-white, RI-brown, VA-white, VA-brown) estimates used to plot the estimates in Figure 3 (dotted lines) were calculated using the same LoopedGetTopt_github.R, but the script was changed so that it loops through by population instead of by individual. 
