%% Generate accuracy plots for each subject
%%
clear all
addpath(genpath('../'))
run loadAllDataIntoTable.m
addpath(genpath('../../monoLS'))
%%
[f1,f2]=accPlots(superSuperT);
saveFig(f1,'../fig/allstatic/',['accuracyAll'],0)
saveFig(f2,'../fig/allstatic/',['accuracySubjectAndBlockEffects'],0)
%%
%All subjects:
[fh,f2]=rtPlots(superSuperT);
%saveFig(fh,'../fig/allstatic/',['rtAll'],0)
saveFig(f2,'../fig/allstatic/',['EZdd'],0)
%w/o subject 2:
%fh=rtPlots(superSuperT(superSuperT.subID~=2,:));
%saveFig(fh,'../../fig/all/',['rtAll'],0)

%%
fh=ssPlots(superSuperT);
%saveFig(fh,'../../fig/all/',['ssAll'],0)

%% Do EZ modeling:
x=superSuperT.pertSize;
y=superSuperT.initialResponse;
y=superSuperT.initialResponse==-1+.5*isnan(superSuperT.initialResponse);

[params, predictedY, Likelihood] = fitGenPsycho(x,y);
params
figure; hold on; G=findgroups(x); plot(splitapply(@(z) mean(z),x,G),splitapply(@(z) mean(z),y,G),'o'); hold on; plot(sort(x),sort(predictedY));