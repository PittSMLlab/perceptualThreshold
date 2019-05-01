# perceptualThreshold experiment

# data organization
Original data from experiment is stored as .mat files (datlogs). mat files are NOT stored in Github because of git issues, but are available at ...
Datlog to .csv conversion is done through the processDatlog.m file in the /data/ folder. Resulting csv files ARE stored on each subject's /data/ subfolder. There are two types of files: the ones named 'BlockX.csv' store a table in csv format, and contain summary information about each of the 24 trials in each block. The files name 'BlockX_yyyy.csv' store some additional information about EACH stride (25 per trial) on each of the 24 trials of the block.

# analysis
First, add the /src/ folder to Matlab's path.
The monoLS toolbox is also needed for some analyses, find at github.com/pabloi/monoLS and add to path.
Basic analysis can be performed through the accuracyPlots.m function, called with the tables stored in csv files. 
Example: t=readtable('../data/AA01/Block1.csv'); accuracyPlots(t);
To generate plots for all subjects, run genAccPlots.m


