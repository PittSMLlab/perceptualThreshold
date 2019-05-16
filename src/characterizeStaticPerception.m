%% Generate accuracy plots for each subject
%%
run loadAllDataIntoTable.m
addpath(genpath('../../monoLS'))
%%
fh=accPlots(superSuperT);
%saveFig(fh,'../../fig/all/',['accuracyAll'],0)
%%
%All subjects:
fh=rtPlots(superSuperT);
%w/o subject 2:
fh=rtPlots(superSuperT(superSuperT.subID~=2,:));
%saveFig(fh,'../../fig/all/',['rtAll'],0)

%%
fh=ssPlots(superSuperT);
%saveFig(fh,'../../fig/all/',['ssAll'],0)