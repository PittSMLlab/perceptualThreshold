%% Generate accuracy plots for each subject
%%
clear all
addpath(genpath('../'))
run loadAllDataIntoTable.m
addpath(genpath('../../monoLS'))
%%
[f1]=accPlots(superSuperT);
ph=findobj(f1,'Type','Axes');
set(ph,'FontSize',10);
for i=1:length(ph)
    ph(i).XAxis.FontSize=8;
    ph(i).YAxis.FontSize=8;
end
%saveFig_(f1,'../fig/allstatic/',['acc'],0)
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
ph=findobj(fh,'Type','Axes');
set(ph,'FontSize',10);
for i=1:length(ph)
    ph(i).XAxis.FontSize=8;
    ph(i).YAxis.FontSize=8;
end
saveFig_(fh,'../fig/allstatic/',['ssAll'],0)

%% EZ modeling
[f1,f2] = EZplots(superSuperT);
extendedPanelWidth(f1,.1)
ph=findobj(f1,'Type','Axes');
set(ph,'FontSize',10);
for i=1:length(ph)
    ph(i).XAxis.FontSize=8;
    ph(i).YAxis.FontSize=8;
end
saveFig_(f1,'../fig/allstatic/',['EZfit'],0)