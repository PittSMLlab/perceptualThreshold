load ../../data/dynamicProfiles.mat
%fh=figure('Units','Normalized','OuterPosition',[.5 .5 .5 .5]);
%Baseline tasks:
pS=[100;-200;0;100;0;200;-100;200;0;-100;0;-200];
sstb=[find(isnan(vLb(1:end-1)) & ~isnan(vLb(2:end)))]; %Start location of each task;
eetb=[find(~isnan(vLb(1:end-1)) & isnan(vLb(2:end)))]; %Start location of each task;
%Adaptation tasks:
pSa=[200;200;400;200;400;200;400;200;400;200;400;200;400;200;400;200;400;200;200];
ssta=[find(isnan(vL(1:end-1)) & ~isnan(vL(2:end)))]; %Start location of each task;
eeta=[find(~isnan(vL(1:end-1)) & isnan(vL(2:end)))]; %Start location of each task;

%Post-adapt tasks:
pSp=[100;0;100;200;0;100;0;100;-100;0;100;-100;0];
sstp=[find(isnan(vLp(1:end-1)) & ~isnan(vLp(2:end)))]; %Start location of each task;
eetp=[find(~isnan(vLp(1:end-1)) & isnan(vLp(2:end)))]; %Start location of each task;

%All together:
breakSize=20;
sst=[sstb; sstb+length(vLb)+breakSize; ssta+2*length(vLb)+2*breakSize; sstp+2*length(vLb)+length(vL)+3*breakSize];
eet=[eetb; eetb+length(vLb)+breakSize; eeta+2*length(vLb)+2*breakSize; eetp+2*length(vLb)+length(vL)+3*breakSize];
pps=[pS; pS; pSa; pSp];

v=[vRb-vLb; nan(breakSize,1); vRb-vLb; nan(breakSize,1); vR-vL; nan(breakSize,1); vRp-vLp];
p1=plot(v,'LineWidth',2,'DisplayName','Set speeds');
hold on
p2=scatter(sst,pps,20,'filled','MarkerEdgeColor','w','DisplayName','Probe start');
for i=1:length(pps)
   pp=patch([sst(i) eet(i) eet(i) sst(i)],500*[-1 -1 1 1],.85*ones(1,3),'EdgeColor','none','FaceAlpha',.6,'DisplayName','Task')
    uistack(pp,'bottom')
end
set(gca,'XLim',[0 max(eet)+10],'YLim',[-300 500])
legend([p1 p2 pp],'Location','South')
ylabel('Belt speed difference (mm/s)')
xlabel('Stride cycles')
title('Perceptual assessment protocol')
set(gca,'FontSize',10)
ph=gca;
ph.XAxis.FontSize=8;
ph.YAxis.FontSize=8;
%saveFig_(fh,'../../fig/','dynamicProtocol',0)