f1=figure('Units','pixels','InnerPosition',[100 100 1.5*300 1*300]);
hold on
x=[-400:400];
y=1-1./(1+exp(x/75));
plot(x,y,'k','LineWidth',2)
ax=gca;
ax.XAxis.TickValues=[0];
ax.XAxis.TickLabels={'PSE'};
ax.XAxis.Label.String='belt speed difference (mm/s)';
ax.YAxis.Label.String='choice ratio';
ax.YAxis.TickValues=[0 .5 1];
ax.YAxis.TickLabels={'100% right','50/50','100% left'};
ax.YAxis.TickLabelRotation=90;
hold on
plot([-400 0],[.5 .5],'k--')
plot([0 0],[0 .5],'k--')
plot(0,.5,'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',7)

set(ax,'FontName','OpenSans')
export_fig ../../fig/PSEdescription.png -png -c[0 5 0 5] -transparent -r600