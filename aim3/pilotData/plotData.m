%% Load data (adaptation)
%load('20180409T154404_PerceptualDynamics_Baseline.mat')
load('20180409T155049_PerceptualDynamics_Adaptation.mat')
load('20180409T161056_PerceptualDynamics_PostAdaptation.mat')

load('20180404T145759_PerceptualDynamics_Adaptation.mat')
%% Summarize datlog
[trialData,strideData]=datlogSummarize(datlog);

%% 
figure
pp=unique(trialData.pertSize);
for i=1:length(pp)
    y1=trialData.pertSize(trialData.pertSize==pp(i));
    y2=trialData.lastSpeedDiff(trialData.pertSize==pp(i));
    x1=trialData.startTime(trialData.pertSize==pp(i));
    x2=trialData.endTime(trialData.pertSize==pp(i));
    subplot(2,2,1)
    hold on
    s1=scatter(x2,y2,'filled');
    s2=scatter(x1,y1,'MarkerEdgeColor',s1.CData);
    plot(x2,y2,'Color',s1.CData)
    plot([x1 x2]',[y1 y2]','k')
    xlabel('Time (s)')
    title('RAW')
    ylabel('Belt-speed diff (mm/s)')

    subplot(2,2,2)
    hold on
    offset=pp(i);
    s1=scatter(x2,y2-offset,'filled','DisplayName', 'End data');
    s2=scatter(x1,y1-offset,'MarkerEdgeColor',s1.CData,'DisplayName','Init. data');
    plot(x2,y2-offset,'Color',s1.CData)
    plot([x1 x2]',[y1 y2]'-offset,'k')
    xlabel('Time (s)')
    title('Full align')
    ylabel('Belt-speed diff (mm/s)')
    legend([s1 s2],'Location','SouthEast')
    
        subplot(2,2,3)
    hold on
    offset=.6*pp(i);
    s1=scatter(x2,y2-offset,'filled');
    s2=scatter(x1,y1-offset,'MarkerEdgeColor',s1.CData);
    plot(x2,y2-offset,'Color',s1.CData)
    plot([x1 x2]',[y1 y2]'-offset,'k')
    xlabel('Time (s)')
    title('60% align')
    ylabel('Belt-speed diff (mm/s)')
    
    y2=trialData.reactionTime(trialData.pertSize==pp(i));
    subplot(2,2,4)
    hold on
    s2=scatter(x2,y2,'filled','MarkerFaceColor',s1.CData);
    %s2=scatter(x1,y1,'MarkerEdgeColor',s1.CData);
    plot(x2,y2,'Color',s1.CData)
    %plot([x1 x2]',[y1 y2]','k')
    xlabel('Time (s)')
    title('RAW times')
    ylabel('Reaction time (s)')
end
