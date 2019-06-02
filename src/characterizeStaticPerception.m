%% Generate accuracy plots for each subject
%%
clear all
addpath(genpath('../'))
run loadAllDataIntoTable.m
addpath(genpath('../../monoLS'))
%%
[f1]=accPlots(superSuperT);
%saveFig(f1,'../fig/allstatic/',['accuracyAll'],0)
%saveFig(f2,'../fig/allstatic/',['accuracySubjectAndBlockEffects'],0)
%%
%All subjects:
[fh,f2]=rtPlots(superSuperT);
%saveFig(fh,'../fig/allstatic/',['rtAll'],0)
%saveFig(f2,'../fig/allstatic/',['EZdd'],0)
%w/o subject 2:
%fh=rtPlots(superSuperT(superSuperT.subID~=2,:));
%saveFig(fh,'../../fig/all/',['rtAll'],0)

%%
fh=ssPlots(superSuperT);
%saveFig(fh,'../../fig/all/',['ssAll'],0)

%% EZ modeling
[f1,f2] = EZplots(superSuperT);