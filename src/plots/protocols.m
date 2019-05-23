%% Static protocol:
pS=[-250;100;25;-50;300;125;-75;-10;150;75;-200;-125;-25;0;-300;10;50;-100;-350;0;200;-150;250;350];
sT=[22.7219950407743;71.8679968267679;122.065994516015;173.413996025920;224.091995134950;273.719999939203;324.947996065021;375.670000165701;426.549999788404;477.708996087313;528.312000632286;578.844999521971;629.546001926064;679.892992973328;731.772001460195;781.585993990302;832.899998500943;883.257993310690;934.044994041324;984.847997501493;1036.44699342549;1086.55499964952;1137.80599907041;1185.73399297893];
eT=[47.3539978265762;97.4669929593801;147.932993620634;199.308001622558;249.031992629170;299.819999188185;350.938993692398;401.860996708274;452.616997808218;503.435993939638;553.886001929641;604.755993559957;655.289002507925;706.236996129155;756.805993616581;807.684001699090;858.911997824907;909.421993419528;959.116996452212;1011.56699359417;1061.84999383986;1112.45999895036;1161.85199618340;1211.59799471498];

%
fh=figure('Units','Normalized','OuterPosition',[.5 .5 .5 .5]);
hold on
sst=round(sT/1.2); %Deriving stride numbers from data. Task was 25 strides to respond, 25 strides of tied, and repeat.
eet=round(eT/1.2);
v=zeros(1,max(eet)+10);
for i=1:length(sst)
   v(sst(i)-1:sst(i))=pS(i);
   v(sst(i)+1:eet(i))=nan;
end
p1=plot(v,'LineWidth',2,'DisplayName','Set speeds');
p2=plot(sst,pS,'o','MarkerFaceColor','red','MarkerEdgeColor','none','DisplayName','Probe start');
for i=1:length(pS)
   pp=patch([sst(i) eet(i) eet(i) sst(i)],400*[-1 -1 1 1],.85*ones(1,3),'FaceAlpha',.9,'EdgeColor','none','FaceAlpha',.6,'DisplayName','Task')
    uistack(pp,'bottom')
end
set(gca,'XLim',[0 max(eet)+10])
legend([p1 pp p2],'Location','North')
ylabel('Belt speed difference (mm/s)')
xlabel('Stride cycles')
title('Perceptual assessment protocol')
saveFig_(fh,'../../fig/','staticProtocol',0)

%% Dynamic protocol:
load ../../data/dynamicProfiles.mat
fh=figure('Units','Normalized','OuterPosition',[.5 .5 .5 .5]);
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
p2=plot(sst,pps,'o','MarkerFaceColor','red','MarkerEdgeColor','none','DisplayName','Probe start');
for i=1:length(pps)
   pp=patch([sst(i) eet(i) eet(i) sst(i)],500*[-1 -1 1 1],.85*ones(1,3),'EdgeColor','none','FaceAlpha',.6,'DisplayName','Task')
    uistack(pp,'bottom')
end
set(gca,'XLim',[0 max(eet)+10],'YLim',[-300 500])
legend([p1 p2 pp],'Location','South')
ylabel('Belt speed difference (mm/s)')
xlabel('Stride cycles')
title('Perceptual assessment protocol')
saveFig_(fh,'../../fig/','dynamicProtocol',0)