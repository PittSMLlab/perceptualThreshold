addpath('../../../ext/altmany-export_fig-b1a7288/')
%%
pS=[-250;100;25;-50;300;125;-75;-10;150;75;-200;-125;-25;0;-300;10;50;-100;-350;0;200;-150;250;350];
sT=[22.7219950407743;71.8679968267679;122.065994516015;173.413996025920;224.091995134950;273.719999939203;324.947996065021;375.670000165701;426.549999788404;477.708996087313;528.312000632286;578.844999521971;629.546001926064;679.892992973328;731.772001460195;781.585993990302;832.899998500943;883.257993310690;934.044994041324;984.847997501493;1036.44699342549;1086.55499964952;1137.80599907041;1185.73399297893];
eT=[47.3539978265762;97.4669929593801;147.932993620634;199.308001622558;249.031992629170;299.819999188185;350.938993692398;401.860996708274;452.616997808218;503.435993939638;553.886001929641;604.755993559957;655.289002507925;706.236996129155;756.805993616581;807.684001699090;858.911997824907;909.421993419528;959.116996452212;1011.56699359417;1061.84999383986;1112.45999895036;1161.85199618340;1211.59799471498];

%
fh=figure('Units','pixels','InnerPosition',[100 100 1.5*300 1*300]);
hold on
sst=round(sT/1.2); %Deriving stride numbers from data. Task was 25 strides to respond, 25 strides of tied, and repeat.
eet=round(eT/1.2);
v=zeros(1,max(eet)+10);
for i=1:length(sst)
   v(sst(i)-1:sst(i))=pS(i);
   v(sst(i)+1:eet(i))=nan;
end
p1=plot(v,'LineWidth',2,'DisplayName','set speeds');
p2=scatter(sst,pS,20,'filled','MarkerEdgeColor','w','DisplayName','probe start');
for i=1:length(pS)
   pp=patch([sst(i) eet(i) eet(i) sst(i)],400*[-1 -1 1 1],.85*ones(1,3),'FaceAlpha',.9,'EdgeColor','none','DisplayName','task')
    uistack(pp,'bottom')
end
set(gca,'XLim',[0 max(eet)+10])
lg=legend([p1 pp p2],'Location','North');
lg.Position(1)=lg.Position(1)+.07;
ylabel('belt speed difference (mm/s)')
xlabel('stride cycles')
title('test block')
set(gca,'FontSize',10)
ph=gca;
ph.XAxis.FontSize=12;
ph.YAxis.FontSize=12;
ph.Position(1)=ph.Position(1)+.03;
%saveFig_(fh,'../../fig/','staticProtocol',0)
export_fig ../../fig/staticBlockMLMC.png -png -c[0 5 0 5] -transparent -r600