addpath('../../../ext/altmany-export_fig-b1a7288/')
%%
%%
fh=figure('Units','pixels','InnerPosition',[100 100 2*300 .8*300]);
ax=axes('Position',[.1 .15 .85 .8])
plot(100*[1 1],[-1 1],'k','LineWidth',2)
hold on
plot([100 300],[1 1],'k--')
plot([100 300],-1*[1 1],'k--')

ax.XAxis.Limits=[0 300];
ax.YAxis.Limits=[-1.5 1.5];
ax.YAxis.TickValues=[];
ax.YAxis.TickLabels={};
ax.XAxis.FontSize=12;
text(200,1.15,'left choice barrier','FontName','OpenSans','FontSize',12)
text(200,-1.15,'right choice barrier','FontName','OpenSans','FontSize',12)
annotation('doublearrow',[.1 .375],[1 1]*.52)
text(5,.15,'non-decision time','FontName','OpenSans','FontSize',12)

%Add some random walks:
v=.01;
s=.1;
for i=1:10
    x=zeros(200,1);
    for t=1:200
        x(t+1)=x(t)+v+s*randn;
    end
   idx=find(abs(x)>1,1,'first');
       cc=.8*ones(1,3);
   plot([0:idx-1]+100,[x(1:idx-1); sign(x(idx))],'Color',cc);
end

load randomWalk.mat
x=x(36:end);
x(36)=0;
idx=find(abs(x)>1,1,'first');
cc=.2*ones(1,3);
plot([0:idx-1]+100,[x(1:idx-1); sign(x(idx))],'Color',cc,'LineWidth',1);
plot(100+idx-1,1,'ko','MarkerFaceColor',zeros(1,3),'MarkerSize',5)

ax.XAxis.TickValues=[0 100+idx];
ax.XAxis.TickLabels={'task start','response time'};
%ax.YAxis.Label.String='evidence accumulation';
ylabel({'evidence';'accumulation'});
ax.YAxis.Label.FontSize=12;
%%
export_fig ../../fig/ddm.png -png -c[0 5 0 5] -transparent -r600