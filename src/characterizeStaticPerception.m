%% Generate accuracy plots for each subject
%%
clear all
addpath(genpath('../'))
run loadAllDataIntoTable.m
addpath(genpath('../../monoLS'))
addpath('../../../ext/altmany-export_fig-b1a7288/')
%% Accuracy plots
close all
[f1,f2]=accPlots(superSuperT);
ph=findobj(f1,'Type','Axes');
set(ph,'FontSize',10);
for i=1:length(ph)
    ph(i).XAxis.FontSize=10;
    ph(i).YAxis.FontSize=10;
    ph(i).Position([2,4])=[.16 .74];
end
%saveFig_(f1,'../fig/allstatic/',['acc'],0)
figure(f1)
%export_fig ../fig/allstatic/acc.eps -eps -c[0 5 0 5] -transparent -m2 -r600 -depsc
close(f2)
export_fig ../fig/allstatic/acc.png -png -c[0 5 0 5] -transparent -r600


%%
close all 
[f1,f2]=accPlots(superSuperT);
ph=findobj(f2,'Type','Axes');
set(ph,'FontSize',10);
for i=1:length(ph)
    ph(i).XAxis.FontSize=10;
    ph(i).YAxis.FontSize=10;
    ph(i).Position([2,4])=[.16 .74];
end
close(f1)
%saveFig_(f2,'../fig/allstatic/',['subjectEffects'],0)
figure(f2)
%export_fig ../fig/allstatic/subjectEffects.eps -eps -c[0 5 0 5] -transparent -r600 -depsc
export_fig ../fig/allstatic/subjectEffects.png -png -c[0 5 0 5] -transparent -r600
%%
%All subjects:
[fh,f2]=rtPlots(superSuperT);
%saveFig(fh,'../fig/allstatic/',['rtAll'],0)
%saveFig(f2,'../fig/allstatic/',['EZdd'],0)
%w/o subject 2:
%fh=rtPlots(superSuperT(superSuperT.subID~=2,:));
%saveFig(fh,'../../fig/all/',['rtAll'],0)

%%
close all
fh=ssPlots(superSuperT);
ph=findobj(fh,'Type','Axes');
set(ph,'FontSize',10);
for i=1:length(ph)
    ph(i).XAxis.FontSize=10;
    ph(i).YAxis.FontSize=10;
    ph(i).Position([2,4])=[.16 .74];
end
%saveFig_(fh,'../fig/allstatic/',['ssAll'],0)
export_fig ../fig/allstatic/ssAll.png -png -c[0 5 0 5] -transparent -r600

%% EZ modeling
[f1,f2] = EZplots(superSuperT);
extendedPanelWidth(f1,.1)
ph=findobj(f1,'Type','Axes');
set(ph,'FontSize',10);
for i=1:length(ph)
    ph(i).XAxis.FontSize=10;
    ph(i).YAxis.FontSize=10;
end
close(f2)
%saveFig_(f1,'../fig/allstatic/',['EZfit'],0)
export_fig ../fig/allstatic/EZfit.png -png  -transparent -r600
cc=get(ph(1),'ColorOrder');
ll=findobj(f1,'Type','Line','Color',cc(1,:));
delete(ll)
export_fig ../fig/allstatic/EZfitBis.png -png  -transparent -r600