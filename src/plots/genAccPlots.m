%% Generate accuracy plots for each subject
%%
loadAllDataIntoTable
%%
fh=accPlots(superSuperT);
saveFig(fh,'../fig/all/',['accuracyAll'],0)
%%

fh=rtPlots(superSuperT);
saveFig(fh,'../fig/all/',['rtAll'],0)

%%
fh=ssPlots(superSuperT);
%saveFig(fh,'../fig/all/',['ssAll'],0)