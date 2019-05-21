%% Task figure
figure('Units','Normalized','OuterPosition',[.7 .6 .3 .4])
v=[zeros(1,5),ones(1,2),nan(1,26),zeros(1,4)];
patch([7 32 32 7],[-1 -1 2 2],.6*ones(1,3),'FaceAlpha',.2,'EdgeColor','none');
hold on
p1=plot(v,'LineWidth',2,'DisplayName','Set speeds');
set(gca,'YLim',[-.2 1.5],'XLim',[1 37])
r1=[nan(1,6),ones(1,6),1-[.04, .08 .1 .1+cumsum(.07*abs(rand(1,17)))]];
r2=[nan(1,6),ones(1,6),1-[-.05 -.09 -.11 -.11+cumsum(.05*abs(rand(1,17)))]];
p2=plot(r1,'Color',.5*ones(1,3),'LineWidth',1,'DisplayName','Sample responses');
plot([32:34],r1(end)*[1,.5,.0],'Color',p1.Color,'LineStyle','--','LineWidth',2)
p3=plot(r2,'Color',.5*ones(1,3),'LineWidth',1);
plot([32:34],r2(end)*[1,.5,.0],'Color',p1.Color,'LineStyle','--','LineWidth',2)
plot([12 12],[-.2 1.5],'k--')
p4=plot(7,1,'o','MarkerFaceColor','r','MarkerEdgeColor','none','DisplayName','Initial Condition');
p5=plot(32,r1(end),'o','MarkerFaceColor','k','MarkerEdgeColor','none','DisplayName','Reported PSE');
p5=plot(32,r2(end),'o','MarkerFaceColor','k','MarkerEdgeColor','none','DisplayName','Reported PSE');
ylabel('Belt speed diff. (mm/s')
set(gca,'XTick',[7,32],'XtickLabel',{'Start cue','End cue','Return to normal'})

%Annotate
h=1.7;
annotation('textarrow',.02+[15.6 14.5]/40,.2/h+[1.2 1.05]/h,'String',{'Incorrect';' response'})
annotation('textarrow',.02+[15.65 14.55]/40,.2/h+[.75 .9]/h,'String',{'Correct  '; '  response'})
annotation('doublearrow',[11 15]/41,.2/h+[.3 .3]/h)
text(7.3,.32,{'Reaction';'    time'})
title('Perceptual task')

saveFig(gcf,'.','PerceptutalTask',0)
legend([p1 p2 p4 p5],'Location','NorthEast')