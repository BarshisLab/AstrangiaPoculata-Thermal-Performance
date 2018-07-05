# Info on analyzing all P and R data and creating Figures 3, 4, and S3

Work Flow for Analyzing Respirometry Data

1. Collect raw oxygen evolution data from Oxy. This is saved as a .txt file. 

2. Correct the raw output files from the Oxy based on the spreadsheet provided by PreSens. This cleans up the file and corrects it using the calculations on the PreSens correction spreadsheet
** Save both the raw .txt and the _cleaned.txt files from the Oxy inside the folder designated for the date of the experiment, included here in the Raw Data folder

3. The photosynthesis and respiration files are separated and added to their own files inside a folder called /RAnalysis/. I use cp at the command line to move the _cleaned files over only

4. Run Segmented_Respirometry_olitotest.R, some parts of code will have to be changed depending whether you are analyzing a cold or heat ramp experiment. The only thing that is changing is the names of variables based on temperature. This Olito script produces corrected photosynthesis and respiration rates (corrected for time and surface area of fragment) and writes out the photosynthesis and respiration data to separate .csv files

5. Copied and pasted all individual ramp days into a master data sheet called: 
~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Respiration_master_holobiont.csv
~/ODU_MS/ODU_MS/ExperimentalTanks/ExperimentData/Photosynthesis_master_holobiont.csv

6. Still in Segmented_Respirometry_olitotest.R, bind together corrected data and important fragment identifying data, calculate gross photosynthesis as Pnet - Rdark. Created these databases separately, they are saved as Heat_Data.csv and Cold_Data.csv

7. Use Respirometry_ANOVA.R script to analyze and plot the rate data from the Segmented_Respirometry_olitotest.R. Because this script analyzes the heat and cold ramp data put together (6 - 32C), we first randomly sample the measurements made at 18C. The randomly sampled metabolic rate files are included here, in the Data For Figures folder. 

