%%
close all
addpath('../../../ext/altmany-export_fig-b1a7288/')
fh=figure('Units','pixels','InnerPosition',[100 100 3*300 1*150]);
%% Static Protocol Schematic
ax=axes('Position',[.1 .5 .85 .4]);
staticProtocolSchematic
%% Perceptual task
ax=axes('Position',[.1 .1 .4 .4]);
taskDescription
%% Block profile
ax=axes('Position',[.55 .1 .4 .4]);
staticBlocks