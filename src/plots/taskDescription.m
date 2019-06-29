addpath('../../../ext/altmany-export_fig-b1a7288/')

%% Task figure
%figure('Units','Normalized','OuterPosition',[.5 .5 .5 .5])
f1=figure('Units','pixels','InnerPosition',[100 100 1.5*300 1*300]);
v=[zeros(1,5),ones(1,2),nan(1,26),zeros(1,6)];
patch([7 32 32 7],[-1 -1 2 2],.85*ones(1,3),'FaceAlpha',.65,'EdgeColor','none');
hold on
p1=plot(v,'LineWidth',2,'DisplayName','set speeds');
set(gca,'YLim',[-.2 1.5],'XLim',[1 38])
x1=    [0.1438
    0.1982
    0.2347
    0.2530
    0.3100
    0.3442
    0.3794
    0.3897
    0.4274
    0.4540
    0.4625
    0.4725
    0.5292
    0.5673
    0.6280
    0.6934
    0.7217]';
x2=    [   -0.0845
   -0.0714
   -0.0597
   -0.0517
   -0.0417
   -0.0060
    0.0169
    0.0199
    0.0216
    0.0489
    0.0620
    0.0828
    0.1178
    0.1250
    0.1271
    0.1333
    0.1705]';
r1=[nan(1,6),ones(1,6),1-[.04, .08 .1 x1]];
r2=[nan(1,6),ones(1,6),1-[-.05 -.09 -.11 x2]];
p2=plot(r1,'Color',.5*ones(1,3),'LineWidth',1,'DisplayName','sample responses');
plot([32:34],r1(end)*[1,.5,.0],'Color',p1.Color,'LineStyle','--','LineWidth',2)
p3=plot(r2,'Color',.5*ones(1,3),'LineWidth',1);
plot([32:34],r2(end)*[1,.5,.0],'Color',p1.Color,'LineStyle','--','LineWidth',2)
plot([12 12],[-.2 1.5],'k--')
set(gca,'ColorOrderIndex',2)
p4=scatter(7,1,30,'filled','MarkerEdgeColor','w','DisplayName','Probe size');
p5=plot(32,r1(end),'o','MarkerFaceColor','k','MarkerEdgeColor','none','DisplayName','reported PSE');
p5=plot(32,r2(end),'o','MarkerFaceColor','k','MarkerEdgeColor','none','DisplayName','reported PSE');
ylabel('belt speed diff. (mm/s)')
set(gca,'XTick',[7,32],'XtickLabel',{'start cue','end cue','return to normal'},'YTick',0,'FontName','OpenSans')

%Annotate
h=1.7;
w=42;
annotation('textarrow',.03+[15.6 14.5]/40,.2/h+[1.15 1.05]/h,'String',{'incorrect';' response'})
annotation('textarrow',.03+[15.5 14.55]/40,.2/h+[.75 .9]/h,'String',{'correct  '; '  response'})
annotation('doublearrow',[11 15]/w,.2/h+[.3 .3]/h)
text(6.7,.335,{'reaction';'  time'},'FontName','OpenSans')
annotation('doublearrow',[8.3 8.3]/w,(.2+[.15 .95])/h)
text(2.9,.3,'probe size','Rotation',90,'FontName','OpenSans')
%annotation('doublearrow',[32 32]/w,(.2+[.15 r1(end)+.08])/h)
%annotation('doublearrow',[31 31]/w,(.2+[.15 r2(end)+.03])/h)
annotation('textarrow',[30 32.4]/w,([.8 r1(end)+.32 ])/h,'String',{'reported';' PSE'})
annotation('arrow',[30.5,32.3]/w,([.87 r2(end)+.18])/h)
%text(28.8,.1,{'Reported PSE'},'Rotation',90,'FontName','OpenSans')
title('perceptual task')

lg=legend([p1 p2],'Location','NorthEast');
lg.Position(1)=lg.Position(1)+.1;
set(gca,'FontSize',10)
ph=gca;
ph.XAxis.FontSize=12;
ph.YAxis.FontSize=12;
%saveFig_(gcf,'../../fig/','PerceptutalTask',0)
export_fig ../../fig/perceptualTask.png -png -c[0 5 0 5] -transparent -r600