%%
addpath(genpath('../src/datlogManipulation/'))
%%
dataDir='../data/';
subList=dir([dataDir 'AB*']);
%% Load fast & slow baseline trials
for j=1:length(subList)
        b1=readtable([dataDir subList(j).name '/Baseline1.csv']);
        b2=readtable([dataDir subList(j).name '/Baseline2.csv']);
        if j>1
           allT=[allT;b1;b2];
        else
            allT=[b1;b2];
        end

     bs=readtable([dataDir subList(j).name '/BaselineSlow.csv']);
    if j>1
       allT1=[allT1;bs];
    else
        allT1=bs;
    end
end

%% Plot baseline performance
figure;
subplot(2,2,1) %Accuracy
B=findgroups(allT.pertSize); %pertSize>0 means vR>vL
B1=findgroups(allT1.pertSize); %pertSize>0 means vR>vL
pp=unique(allT.pertSize);
pp1=unique(allT1.pertSize);
S=splitapply(@(x) (sum(x==-1)+.5*sum(isnan(x)))/length(x),allT.initialResponse,B); %Counting LEFT IS SLOW choices plus HALF of no response
scatter(pp,S,50,pp,'filled')
hold on
S1=splitapply(@(x) (sum(x==-1)+.5*sum(isnan(x)))/length(x),allT1.initialResponse,B1); %Counting LEFT IS SLOW choices plus HALF of no response
scatter(pp1,S1,20,.4*ones(1,3),'filled')
grid on
title('<- Choices')

subplot(2,2,2) %Clickrates
S=splitapply(@(x) mean(x),allT.Lclicks-allT.Rclicks,B); %Counting LEFT IS SLOW choices plus HALF of no response
scatter(pp,S,50,pp,'filled')
grid on
title('Net clicks')

subplot(2,2,3) %Clickrates
S=splitapply(@(x) mean(x),allT.lastSpeedDiff,B); %Counting LEFT IS SLOW choices plus HALF of no response
scatter(pp,S,50,pp,'filled')
grid on
title('Last diff (mm/s)')

%% Get adapt trials
for j=1:length(subList)
    td=readtable([dataDir subList(j).name '/Adaptation.csv']);
    if j>1
       allTA=[allTA;td];
    else
        allTA=td;
    end
end
%% Plot adaptation
figure;
vars={'lastSpeedDiff','netClicks','reactionTime'};
allTA.netClicks=allTA.Lclicks-allTA.Rclicks;
allT.netClicks=allT.Lclicks-allT.Rclicks;
for k=1:3
    v=allTA.(vars{k});
    vB=allT.(vars{k});

    subplot(3,1,k)
    hold on
    aux=mod([1:length(allTA.startTime)]-1,19);
    scatter(aux(allTA.pertSize==200),v(allTA.pertSize==200),50,[1,0,0])
    scatter(aux(allTA.pertSize==400),v(allTA.pertSize==400),50,[0,0,1])
    for i=0:18
        scatter(i,mean(v(aux==i)),100,[-1,0,1]*(mean(allTA.pertSize(aux==i))/200-1) +[1,0,0],'filled')
    end
    scatter(-2,mean(vB(allT.pertSize==200)),100,[.7,.2,0],'filled')
    scatter(-2*ones(sum(allT.pertSize==200),1),vB(allT.pertSize==200),50,[.7,.2,0])
    axis tight
    pp=patch([.8 17.8 17.8 .8],[min(v)*[1 1] max(v)*[1 1]],.6*ones(1,3),'FaceAlpha',.3,'EdgeColor','none');
    uistack(pp,'bottom')
    grid on
    title(vars{k})
end

%%
subplot(3,1,2)
%scatter(allTA.startTime(allTA.pertSize==200),200,50,[1,0,0],'filled')
hold on
scatter(aux(allTA.pertSize==200),allTA.Lclicks(allTA.pertSize==200)-allTA.Rclicks(allTA.pertSize==200),50,[1,0,0])
%scatter(allTA.startTime(allTA.pertSize==400),400,50,[0,0,1],'filled')
scatter(aux(allTA.pertSize==400),allTA.Lclicks(allTA.pertSize==400)-allTA.Rclicks(allTA.pertSize==400),50,[0,0,1])
for i=0:18
    scatter(i,mean(allTA.Lclicks(aux==i)-allTA.Rclicks(aux==i)),100,[-1,0,1]*(mean(allTA.pertSize(aux==i))/200-1) +[1,0,0],'filled')
end
grid on
title('Net clicks')

%
subplot(3,1,3)
hold on
scatter(aux(allTA.pertSize==200),allTA.reactionTime(allTA.pertSize==200),50,[1,0,0])
scatter(aux(allTA.pertSize==400),allTA.reactionTime(allTA.pertSize==400),50,[0,0,1])
for i=0:18
    scatter(i,nanmean(allTA.reactionTime(aux==i)),100,[-1,0,1]*(mean(allTA.pertSize(aux==i))/200-1) +[1,0,0],'filled')
end
title('Reaction times')
grid on