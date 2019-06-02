function [fh,fh2]=trialsSummaryPlots(trialData)
colorOff=3;
%%
B=findgroups(trialData.pertSize); %pertSize>0 means vR>vL
pp=unique(trialData.pertSize);
cmap=parula(length(pp)); %parula, jet, redbluecmap
cmap=cmap*.8;
 
%% Define quantities of interest
%nullTrials=trialData.pertSize==0;
%correctResponses=trialData.initialResponse==-sign(trialData.pertSize) & ~nullTrials; %Negative response means LEFT IS SLOW (RIGHT IS FAST) choice
%nonResponse=isnan(trialData.initialResponse) & ~nullTrials;
%incorrectResponses=trialData.initialResponse==sign(trialData.pertSize) & ~nullTrials;

%Adding prev perturbation to table:
trialData.prevSize=[0;trialData.pertSize(1:end-1)]; 
trialData.prevSize(trialData.pertSize==-250*((-1).^trialData.blockNo))=0; %Assigning NaN to previous perturbation for first trial in each block (+250 in even blocks, -250 in odd)
trialData.prevFinalSpeed=[0;trialData.lastSpeedDiff(1:end-1)].*[0;diff(trialData.subID)==0].*[0;diff(trialData.blockNo)==0];

%Creating binary response variable(s):
trialData.leftResponse=trialData.initialResponse==-1;
trialData.rightResponse=trialData.initialResponse==1;
trialData.noResponse=isnan(trialData.initialResponse);
trialData.nullTrials=trialData.pertSize==0;
trialData.correctResponses=trialData.initialResponse==-sign(trialData.pertSize) & ~trialData.nullTrials;
trialData.incorrectResponses=trialData.initialResponse==sign(trialData.pertSize) & ~trialData.nullTrials;

%Creating a (modified) subject ID field: %See comment at end about ANOVA
aux=trialData.subID;
%aux(aux==1)=10;
%aux(aux==4)=1;
trialData.ID=categorical(aux); 

%% First figure: global stats
fh=figure('Units','Normalized','OuterPosition',[.5 0 .5 1]);
Q=6;

%% First 3 plots: some summary stats
subplot(1,3,1)
c=[sum(trialData.pertSize<0) sum(trialData.nullTrials) sum(trialData.pertSize>0) ]/numel(trialData.nullTrials);
bar(3,c(:,1),'FaceColor',cmap(end-colorOff+1,:),'EdgeColor','none')
hold on
bar(2,c(:,2),'FaceColor',cmap(round(size(cmap,1)/2),:),'EdgeColor','none')
bar(1,c(:,3),'FaceColor',cmap(colorOff,:),'EdgeColor','none')
set(gca,'XTick',1:3,'XTickLabel',{'vL>vR','NT','vL<vR'},'YLim',[0 .5])
title('Trial decomposition')
grid on

subplot(1,3,2)
d=[sum(trialData.initialResponse==1 & trialData.nullTrials) sum(trialData.initialResponse==-1 & trialData.nullTrials)  sum(trialData.noResponse & trialData.nullTrials)]/sum(trialData.nullTrials);
bar(1,d(:,1),'FaceColor',cmap(colorOff,:),'EdgeColor','none')
hold on
bar(2,d(:,2),'FaceColor',cmap(end-colorOff+1,:),'EdgeColor','none')
bar(3,d(:,3),'FaceColor',cmap(round(size(cmap,1)/2),:),'EdgeColor','none')
hold on
set(gca,'XTick',1:3,'XTickLabel',{'>','<','NR'})
title('Null trial responses')
grid on

subplot(1,3,3)
d=[sum(trialData.correctResponses) sum(trialData.incorrectResponses)  sum(trialData.noResponse)]/sum(~trialData.nullTrials);
dN=[sum(trialData.noResponse) sum(trialData.nullTrials)]/numel(trialData.nullTrials);
d2=[sum(trialData.correctResponses & (trialData.pertSize<0))/sum(trialData.pertSize<0)];
d3=sum(trialData.correctResponses & (trialData.pertSize>0))/sum(trialData.pertSize>0);
bar(1:3,d,'FaceColor',.4*ones(1,3),'EdgeColor','none')
hold on
b=bar([6],[d2]','FaceColor',cmap(end-colorOff+1,:),'EdgeColor','none');
b=bar([5],[d3]','FaceColor',cmap(colorOff,:),'EdgeColor','none');
set(gca,'XTick',1:6,'XTickLabel',{'OK','BAD','NR','','>','<'})
grid on
title('Non-null trial accuracy')
end
