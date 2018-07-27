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
%% Get post trials
for j=1:length(subList)
    td=readtable([dataDir subList(j).name '/PostAdaptation.csv']);
    if j>1
       allTP=[allTP;td];
    else
        allTP=td;
    end
end

%% Plot adaptation
func=@(x) nanmedian(x);
figure;
vars={'lastSpeedDiff','netClicks','reactionTime'};
allTA.netClicks=allTA.Lclicks-allTA.Rclicks;
allT.netClicks=allT.Lclicks-allT.Rclicks;
allTP.netClicks=allTP.Lclicks-allTP.Rclicks;
colors=get(gca,'ColorOrder');
for k=1:3
    v=allTA.(vars{k});
    vB=allT.(vars{k});
    vP=allTP.(vars{k});

    subplot(3,1,k)
    hold on
    %Adapt:
    aux=mod([1:length(allTA.startTime)]-1,19);
    pp=unique(allTA.pertSize);
    for j=1:length(pp)
    scatter(aux(allTA.pertSize==pp(j)),v(allTA.pertSize==pp(j)),50,colors(pp(j)/100+3,:))
    end
    %scatter(aux(allTA.pertSize==400),v(allTA.pertSize==400),50,[0,0,1])
    for i=0:18
        scatter(i,func(v(aux==i)),100,colors(allTA.pertSize(i+1)/100 + 3,:),'filled')
    end
    
    %Base:
    pp=unique(allT.pertSize);
    for j=1:length(pp)
        scatter(-10+j,func(vB(allT.pertSize==pp(j))),50,colors(pp(j)/100+3,:),'filled')
        scatter((-10+j)*ones(sum(allT.pertSize==pp(j)),1),(vB(allT.pertSize==pp(j))),50,colors(pp(j)/100+3,:))
    end
    %scatter(-2,func(vB(allT.pertSize==200)),100,[.7,.2,0],'filled')
    %scatter(-2*ones(sum(allT.pertSize==200),1),vB(allT.pertSize==200),50,[.7,.2,0])
    
    %Paraphernalia
    axis tight
    pp=patch([.8 17.8 17.8 .8],[min(v)*[1 1] max(v)*[1 1]],.6*ones(1,3),'FaceAlpha',.3,'EdgeColor','none');
    uistack(pp,'bottom')
    grid on
    title(vars{k})
    
    % Add post adaptation:
    aux=mod([1:length(allTP.startTime)]-1,13)+19;
    pp=unique(allTP.pertSize);
    pp=pp(1:4);
    for j=1:length(pp)
        scatter(aux(allTP.pertSize==pp(j)),vP(allTP.pertSize==pp(j)),50,colors(pp(j)/100+3,:))
    end
    for i=19:13+19
        scatter(i,func(vP(aux==i)),100,colors(allTP.pertSize(i-18)/100 + 3,:),'filled')
    end
    %for i=0:12
    %    scatter(i,func(vP(aux==i)),100,colors(allTP.pertSize(i+1)/100 + 3,:),'filled')
    %end
end


%% Plot all click rates vs. reaction times
figure;
subplot(2,1,1)
hold on
for ps=[-200:100:200]
    aux=allT.pertSize==ps;
s=scatter((allT.reactionTime(aux)),(abs(allT.Lclicks(aux)-allT.Rclicks(aux))),20,'filled');
scatter(nanmean((allT.reactionTime(aux))),mean(abs(allT.Lclicks(aux)-allT.Rclicks(aux))),100,s.CData,'filled');
end

for ps=[200,400]
    aux=allTA.pertSize==ps;
s=scatter((allTA.reactionTime(aux)),(abs(allTA.Lclicks(aux)-allTA.Rclicks(aux))),20);
scatter(nanmean((allTA.reactionTime(aux))),mean(abs(allTA.Lclicks(aux)-allTA.Rclicks(aux))),100,s.CData);
end
