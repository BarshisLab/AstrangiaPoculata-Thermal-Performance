# Info on creating Figure S1 
-Fig S1 is a plot of the temperature recorded by hobo loggers during all 16 thermal ramp experiments by group (heat holobiont, cold holobiont, heat skeleton and cold skeleton).

-First, use hobocleaner.py script to clean the raw hobo logger files (.txt). 

-Then, trim the cleaned files to only include the times during which an experiment was being run.

-Next, organize all cleaned.txt files by type of ramp experiment.

-Now, can use the FigS1_PlotRampHoboLoggers.R script to plot all temperature data by category of ramp and R/Illustrator to combine all plots into the figure. 
