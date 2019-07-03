%% Dynamic Protocol Schematic Option 2

fh=figure('Units','pixels','InnerPosition',[100 100 3*300 1*150]);
% out=[zeros(1,10),NaN(1,10),zeros(1,40),NaN(1,40)];
% plot(out)
%patch([0 10 10 0],[-1 -1 1 1],.6*ones(1,3),'FaceAlpha',.2,'EdgeColor','k','LineStyle', '--');
set(gca,'YLim',[-1 1],'XLim',[0 140])
set(gca,'Position',[.05 .1 .9 .7])
siz=20;

patch([0 19.9 19.9 0],[-1 -1 1 1],0.75*ones(1,3),'FaceAlpha',1,'EdgeColor','k','LineWidth', 2);
patch([20.1 39.9 39.9 20.1],[-1 -1 1 1],0.75*ones(1,3),'FaceAlpha',1,'EdgeColor','k','LineWidth', 2, 'LineStyle', '--');
patch([40.1 59.9 59.9 40.1],[-1 -1 1 1],ones(1,3),'FaceAlpha',1,'EdgeColor','k','LineWidth', 2);
patch([60.1 79.9 79.9 60.1],[-1 -1 1 1],ones(1,3),'FaceAlpha',1,'EdgeColor','k','LineWidth', 2);
patch([80.1 100 100 80.1],[-1 -1 1 1],ones(1,3),'FaceAlpha',1,'EdgeColor','k','LineWidth', 2);
patch([99.1 119.9 119.9 99.1],[-1 -1 1 1],0.1*ones(1,3),'FaceAlpha',1,'EdgeColor','k','LineWidth', 2);
patch([120.1 140 140 120.1],[-1 -1 1 1],ones(1,3),'FaceAlpha',1,'EdgeColor','k','LineWidth', 2);

%Annotations

x=siz./2;

text(x,0.4,{'Familiarization'},'FontName','OpenSans', 'HorizontalAlignment', 'center','FontWeight', 'bold');
text(x,0.1,{'(with visual'},'FontName','OpenSans','FontSize',8, 'HorizontalAlignment', 'center','FontWeight', 'bold');
text(x,-0.2,{'feedback)'},'FontName','OpenSans','FontSize',8, 'HorizontalAlignment', 'center','FontWeight', 'bold');
%text(0.5,-0.3,{'(Visual Feedback)'},'FontName','OpenSans', 'FontSize',8);

x=x+siz;

text(x,0.4,{'Familiarization'},'FontName','OpenSans', 'HorizontalAlignment', 'center','FontWeight', 'bold');
text(x,0.1,{'(with visual'},'FontName','OpenSans','FontSize',8, 'HorizontalAlignment', 'center','FontWeight', 'bold');
text(x,-0.2,{'feedback'},'FontName','OpenSans','FontSize',8, 'HorizontalAlignment', 'center','FontWeight', 'bold');
text(x,-0.5,{'-as needed)'},'FontName','OpenSans','FontSize',8, 'HorizontalAlignment', 'center','FontWeight', 'bold');
%text(0.5,-0.3,{'(Visual Feedback)'},'FontName','OpenSans', 'FontSize',8);


x=x+siz;

text(x,0.4,{'Baseline'},'FontName','OpenSans', 'HorizontalAlignment', 'center','FontWeight', 'bold');
text(x,0.1,{'Long'},'FontName','OpenSans','HorizontalAlignment', 'center','FontWeight', 'bold');

%text(x,-0.2,{'(as needed)'},'FontName','OpenSans', 'FontSize',8,'HorizontalAlignment', 'center','FontWeight', 'bold');

x=x+siz;

text(x,0.4,{'Baseline'},'FontName','OpenSans', 'HorizontalAlignment', 'center','FontWeight', 'bold');
text(x,0.1,{'Short'},'FontName','OpenSans', 'HorizontalAlignment', 'center','FontWeight', 'bold');
text(x,-0.2,{'1'},'FontName','OpenSans', 'HorizontalAlignment', 'center','FontWeight', 'bold');

x=x+siz;

text(x,0.4,{'Baseline'},'FontName','OpenSans', 'HorizontalAlignment', 'center','FontWeight', 'bold');
text(x,0.1,{'Short'},'FontName','OpenSans', 'HorizontalAlignment', 'center','FontWeight', 'bold');
text(x,-0.2,{'2'},'FontName','OpenSans', 'HorizontalAlignment', 'center','FontWeight', 'bold');

x=x+siz;

text(x,0.4,{'Adaptation'},'FontName','OpenSans', 'Color', 'w','HorizontalAlignment', 'center','FontWeight', 'bold');
text(x,0.1,{'(Inter-task'},'FontName','OpenSans', 'FontSize',8,'Color', 'w','HorizontalAlignment', 'center','FontWeight', 'bold');
text(x,-0.2,{'walking at'},'FontName','OpenSans', 'FontSize',8, 'Color', 'w','HorizontalAlignment', 'center','FontWeight', 'bold');
text(x,-0.5,{'\DeltaV=500)'},'FontName','OpenSans', 'FontSize',8, 'Color', 'w','HorizontalAlignment', 'center','FontWeight', 'bold');

x=x+siz;

text(x,0.4,{'Washout'},'FontName','OpenSans', 'HorizontalAlignment', 'center','FontWeight', 'bold');
%text(x,0.1,{'Adaptation'},'FontName','OpenSans', 'HorizontalAlignment', 'center','FontWeight', 'bold');

% Title, Ticks, etc

title('Perceptual Dynamics Protocol')

set(gca,'FontSize',10)
ph=gca;
set(gca,'xtick',[])
set(gca,'ytick',[])

export_fig ../../fig/dynamicProtocol.png -png -c[0 5 0 5] -transparent -r600