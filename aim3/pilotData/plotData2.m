%% Load data (adaptation)
%load('20180409T154404_PerceptualDynamics_Baseline.mat')
clear filename
%Pilot 3:
filename{1}='20180409T155049_PerceptualDynamics_Adaptation.mat';
filename{2}='20180409T161056_PerceptualDynamics_PostAdaptation.mat';

%Pilot 2:
%filename{1}='20180404T145759_PerceptualDynamics_Adaptation.mat';
%filename(2)=[];
%Pilot 1:
%filename{1}='20180404T140614_PerceptualDynamics_Adaptation.mat';
%filename{2}='20180404T142714_PerceptualDynamics_PostAdaptation.mat';
%filename(2)=[];
%% 
figure('Name',filename{1})
t0=0;
for k=1:length(filename)
    load(filename{k})
[trialData,strideData]=datlogSummarize(datlog);

%%
subplot(2,4,1:4)
hold on
cc=get(gca,'ColorOrder');
l1=plot(datlog.TreadmillCommands.read(:,4)+t0,datlog.TreadmillCommands.read(:,1)-datlog.TreadmillCommands.read(:,2),'DisplayName','Actual \Delta v','Color',cc(1,:));
pp=unique(trialData.pertSize);
if k==1
    A=mode(datlog.TreadmillCommands.read(:,1)-datlog.TreadmillCommands.read(:,2)); %Estimated speed diff amplitude:
end
%Adding expected progression:
t=[0:max(datlog.TreadmillCommands.read(:,4))];
f=.25;
p1=@(t) f*A*((2-k)+exp(-t/150)*(-1)^k);
l2=plot(t+t0,p1(t),'r','DisplayName','Expected perceptual SS');
f=.5;
p2=@(t) f*A*((2-k)+exp(-t/150)*(-1)^k);
l3=plot(t+t0,p2(t),'k','DisplayName','Adj. Expected perceptual SS');
    
for i=1:length(pp)
    subplot(2,4,1:4)
    hold on
    y1=trialData.pertSize(trialData.pertSize==pp(i));
    y2=trialData.lastSpeedDiff(trialData.pertSize==pp(i));
    x1=trialData.startTime(trialData.pertSize==pp(i));
    x2=trialData.endTime(trialData.pertSize==pp(i));
    
    s1=scatter(x2+t0,y2,50,'filled','DisplayName','Task init');
    s2=scatter(x1+t0,y1,50,'MarkerEdgeColor',s1.CData,'DisplayName', 'Task end');
    %plot(x2,y2,'Color',s1.CData)
    %plot([x1 x2]',[y1 y2]','k')
    xlabel('Time (s)')
    title('RAW')
    ylabel('Belt-speed diff (mm/s)')
    
    %add corrected amounts
    subplot(2,4,6+2*(k-1))
    hold on
    bar(x1+t0,y2-y1,'BarWidth',.4,'EdgeColor','none','FaceColor',s1.CData,'FaceAlpha',.5)
    title('Corrected amount')
    e1=p1(x1)-y1;
    e2=p2(x1)-y1;
    %plot(x1+t0,.9*(e1),'r')
    plot(x1+t0,.9*(e2),'Color',s1.CData,'LineWidth',2)
    xlabel('Time (s)')
    ylabel('Change in \Delta v (mm/s)')

    %add reaction times
    y2=trialData.reactionTime(trialData.pertSize==pp(i));
    subplot(2,4,5+2*(k-1))
    hold on
    scatter(x1,y2,'filled','MarkerFaceColor',s1.CData);
    %s2=scatter(x1,y1,'MarkerEdgeColor',s1.CData);
    %plot(x1,y2,'Color',s1.CData)
    %plot([x1 x2]',[y1 y2]','k')
    xlabel('Time (s)')
    title(['RAW times ' filename{k}(36:end-4)])
    ylabel('Reaction time (s)')
    %plot(x1,1./(.00325*abs(e1)+0.1269),'r')
    plot(x1,1./(.0023*abs(e2)+0.1269),'Color',s1.CData,'LineWidth',2)
    
end
t0=t(end);
subplot(2,4,1:4)
legend([l1 l2 l3 s1 s2])
end