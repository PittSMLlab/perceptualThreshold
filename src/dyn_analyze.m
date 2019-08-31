%%
%addpath(genpath('../src/datlogManipulation/'))
%%
 addpath(genpath('../'))
dataDir='../data/';
subList=dir([dataDir 'AB*']);
%% Load fast & slow baseline trials
sL=1:10;
%sL=[1:7,9,10]; %Excluding subj 8
for j=sL
        b1=readtable([dataDir subList(j).name '/Baseline1.csv']); %Short response time task (8 strides!)
        b1.blockNo=ones(size(b1,1),1);
        b1.subID=j*ones(size(b1,1),1);
        b2=readtable([dataDir subList(j).name '/Baseline2.csv']); 
        b2.blockNo=2*ones(size(b2,1),1);
        b2.subID=j*ones(size(b1,1),1);
        if j>1
           allT=[allT;b1;b2];
        else
           allT=[b1;b2];
        end

     bs=readtable([dataDir subList(j).name '/BaselineSlow.csv']);  %Long-response time task (24 strides)
    if j>1
        allT1=[allT1;bs];
    else
        allT1=bs;
    end
end

%% Plot baseline performance
accPlots(allT)
%rtPlots(allT)

%% Alt:
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
xlabel('Initial speed diff. (mm/s)')
ylabel('Left choice (%)')

subplot(2,2,2) %Clickrates
S=splitapply(@(x) mean(x),allT.Lclicks-allT.Rclicks,B); 
scatter(pp,S,50,pp,'filled') %Fast baseline
hold on
S1=splitapply(@(x) mean(x),allT1.Lclicks-allT1.Rclicks,B1);  %Does this mean that subjects are better at the short task? Is it the urgency? Or they being used to the task?
scatter(pp1,S1,20,.4*ones(1,3),'filled') %Slow baseline
grid on
title('Net clicks')

subplot(2,2,3) %last speed diff reached
S=splitapply(@(x) mean(x),allT.lastSpeedDiff,B); 
scatter(pp,S,50,pp,'filled')
hold on
S1=splitapply(@(x) mean(x),allT1.lastSpeedDiff,B1); 
scatter(pp1,S1,20,.4*ones(1,3),'filled')
grid on
title('Last diff (mm/s)')

subplot(2,2,4) %last speed diff reached
nonlin=@(x,alpha) sign(x).*abs(x).^alpha;
nonlininv=@(y,alpha) sign(y).*abs(y).^(1/alpha);
S=splitapply(@(x) mean(x),-allT.lastSpeedDiff+allT.pertSize,B); 
scatter(pp,S,50,pp,'filled')
hold on
S1=splitapply(@(x) mean(x),-allT1.lastSpeedDiff+allT1.pertSize,B1); 
scatter(pp1,S1,20,.4*ones(1,3),'filled')
hold on;
plot([-200 200],[-200 200],'k--')
%p=fitPsycho(pp,S);
alpha=1;%1.4;
k1=110*(1/75).^alpha; %projection constant. I presume from baseline data that subjects that correct 150mm/s had a PSE actually of 200mm/s away from the starting point
plot(k1*nonlin([-200:200],alpha),[-200:1:200],'r')
%plot([-200:200],300*(psycho([0,85],[-200:200])-.5),'r')
grid on
title('Total correction (mm/s)')
xlabel('Init diff (mm/s)')
ylabel('Corrected amount')

%% Get adapt trials
for j=sL
    td=readtable([dataDir subList(j).name '/Adaptation.csv']);
    if j>1
       allTA=[allTA;td];
    else
        allTA=td;
    end
end
%% Get post trials
for j=sL
    td=readtable([dataDir subList(j).name '/PostAdaptation.csv']);
    if j>1
       allTP=[allTP;td];
    else
        allTP=td;
    end
end

%% Plot adaptation
%func=@(x) nanmedian(x);
func=@(x) nanmean(x);
allTA.leftResponse=(allTA.initialResponse==-1) + .5*isnan(allTA.initialResponse); %No response is coded as 50/50
allTP.leftResponse=(allTP.initialResponse==-1) + .5*isnan(allTP.initialResponse); %No response is coded as 50/50
allT.leftResponse=(allT.initialResponse==-1)+.5*isnan(allT.initialResponse);

fh=figure('Units','Normalized','OuterPosition',[.5 .2 .5 .8]);

vars={'lastSpeedDiff','projectedPSE','leftResponse','reactionTime'}; %Variables to plot
%Add: projection to estimate PSE, assuming that subjects undershoot the
%real 'PSE' target. The problem with this approach is that the true PSE
%target is unknown, so the projection needs to be made from the actual
%correction, which is noisier.
allTA.projectedPSE=allTA.pertSize+k1*nonlin(allTA.lastSpeedDiff-allTA.pertSize,alpha);
allTP.projectedPSE=allTP.pertSize+k1*nonlin(allTP.lastSpeedDiff-allTP.pertSize,alpha);
allT.projectedPSE=allT.pertSize+k1*nonlin(allT.lastSpeedDiff-allT.pertSize,alpha);

u=[0 85];
allTA.projectedPSE2=allTA.pertSize+invpsycho(u,(allTA.lastSpeedDiff-allTA.pertSize)/300+.5);
allTP.projectedPSE2=allTP.pertSize+invpsycho(u,(allTP.lastSpeedDiff-allTP.pertSize)/300+.5);
allT.projectedPSE2=allT.pertSize+invpsycho(u,(allT.lastSpeedDiff-allT.pertSize)/300+.5);

%Add: projection to estimate PSE, assuming that subjects prefer to go back
%to what they were just doing

%Add: accuracy of responses (avg. only): can it be used to infer the
%population PSE?

%Add: infer initial distance to PSE from reaction times?

%Get stride coordinates:
load dynamicProfiles.mat %Storage for dynamic profiles: need to change to csv
firstResponseStrideA=find(isnan(vL(2:end)) & ~isnan(vL(1:end-1)));
firstResponseStrideP=find(isnan(vLp(2:end)) & ~isnan(vLp(1:end-1)))+length(vL);

allTA.netClicks=allTA.Lclicks-allTA.Rclicks;
allT.netClicks=allT.Lclicks-allT.Rclicks;
allTP.netClicks=allTP.Lclicks-allTP.Rclicks;
colors=get(gca,'ColorOrder');
%colors=zeros(size(colors));
mK=4;
%allT=allT(~isnan(allT.initialResponse),:);
%allTA=allTA(~isnan(allTA.initialResponse),:);
%allTP=allTP(~isnan(allTP.initialResponse),:);
names={'Reported PSE','Estimated PSE','< response %','Reaction Time'};
for k=1:mK
    v=allTA.(vars{k});
    vB=allT.(vars{k});
    vP=allTP.(vars{k});

    subplot(mK,3,[2:3]+k*3-3)
    hold on
    %Adapt:
    aux=mod([1:length(allTA.startTime)]-1,19);
    pp=unique(allTA.pertSize);
    pp=200;
    for j=1:length(pp)
        x=aux(allTA.pertSize==pp(j));
        y=v(allTA.pertSize==pp(j));
        %scatter(x,y,10,.4*ones(1,3),'filled')
        x=reshape(x,numel(x)/length(sL),length(sL));
        y=reshape(y,numel(x)/length(sL),length(sL));
        %plot(x,y,'Color',.4*ones(1,3))
    end
    %scatter(aux(allTA.pertSize==400),v(allTA.pertSize==400),50,[0,0,1])
    for i=0:18 %19 repetitions of the task
        %if allTA.pertSize(i+1)==200
        if false %k==2
             ss(i+1)=scatter(firstResponseStrideA(i+1),func(v(aux==i)),50,'k','filled','DisplayName',['Probe =' num2str(allTA.pertSize(i+1))]);
        else
            ss(i+1)=scatter(firstResponseStrideA(i+1),func(v(aux==i)),50,colors(allTA.pertSize(i+1)/100 + 3,:),'filled','DisplayName',['Probe =' num2str(allTA.pertSize(i+1))]);
        end
        errorbar(firstResponseStrideA(i+1),func(v(aux==i)),nanstd(v(aux==i))/sqrt(sum((aux==i))),'Color','k')
        %end
    end
   
    %Paraphernalia
    axis tight
    aStart=find((vR-vL)==500,1,'first');
    aEnd=find((vR-vL)==500,1,'last');
    pp=patch([aStart aEnd aEnd aStart],[min([v;vP])*[1 1] max([v;vP])*[1 1]],.6*ones(1,3),'FaceAlpha',.3,'EdgeColor','none');
    uistack(pp,'bottom')
    grid on
    title(vars{k})
    
    % Add post adaptation:
    aux=mod([1:length(allTP.startTime)]-1,13)+19;
    pp=unique(allTP.pertSize);
    pp=pp(1:4);
    for j=1:length(pp)
        %scatter(aux(allTP.pertSize==pp(j)),vP(allTP.pertSize==pp(j)),10,.4*ones(1,3),'filled')
    end
    for i=19:12+19
        if false %k==2
            ss(i+1)=scatter(firstResponseStrideP(i-18),func(vP(aux==i)),50,'k','filled','DisplayName',['Probe =' num2str(allTP.pertSize(i-18))]);
        else
            ss(i+1)=scatter(firstResponseStrideP(i-18),func(vP(aux==i)),50,colors(allTP.pertSize(i-18)/100 + 3,:),'filled','DisplayName',['Probe =' num2str(allTP.pertSize(i-18))]);
        end
        %scatter(i,func(vP(aux==i)),50,'k','filled')
        errorbar(firstResponseStrideP(i-18),func(vP(aux==i)),nanstd(vP(aux==i))/sqrt(sum((aux==i))),'Color','k')
    end
    %for i=0:12
    %    scatter(i,func(vP(aux==i)),100,colors(allTP.pertSize(i+1)/100 + 3,:),'filled')
    %end
    set(gca,'XTick',[-8,10,25],'XTickLabel',{'Baseline','Adapt','Post'})
    if k==1
        legend(ss([1,3,20,21,28]))
    end
    title(names{k})
    if k<=2
        ylabel('\Delta v (mm/s)')
    end 
    
    subplot(mK,3,1+k*3-3)
    %Base:
    hold on
    pp=unique(allT.pertSize);
    for j=1:length(pp)
        if false %k==1
        scatter(pp(j),func(vB(allT.pertSize==pp(j))),50,colors(pp(j)/100+3,:),'filled')
        else
        scatter(pp(j),func(vB(allT.pertSize==pp(j))),50,'k','filled')
        end
        %scatter((-10+j)*ones(sum(allT.pertSize==pp(j)),1),(vB(allT.pertSize==pp(j))),10,.4*ones(1,3),'filled')
        errorbar(pp(j),func(vB(allT.pertSize==pp(j))),nanstd(vB(allT.pertSize==pp(j)))/sqrt(sum(allT.pertSize==pp(j))),'Color','k')
    end
    title('Baseline performance')
    
    set(gca,'XTick',[-200:100:200],'XTickLabel',{'-200','-100','0','100','200'},'XLim',[-250 250])
    grid on
    xlabel('Probe (mm/s)')
    if k==1
        %Adding the baseline fit curve
        x=[-200:100:200];
        y=x+nonlininv(-x,alpha)/k1;
        p1=plot(x,y,'Color',.5*ones(1,3),'LineWidth',2);
        uistack(p1,'bottom')
        
    end
    %scatter(-2,func(vB(allT.pertSize==200)),100,[.7,.2,0],'filled')
    %scatter(-2*ones(sum(allT.pertSize==200),1),vB(allT.pertSize==200),50,[.7,.2,0])
end

%saveFig(fh,'../fig/alldyn/','timecourse',0)
%% Find temporal decay during post-adapt
figure; 
vars={'projectedPSE','lastSpeedDiff','projectedPSE2'};
load dynamicProfiles.mat %Storage for dynamic profiles: need to change to csv
firstResponseStrideA=find(isnan(vL(2:end)) & ~isnan(vL(1:end-1)));
firstResponseStrideP=length(vL)+find(isnan(vLp(2:end)) & ~isnan(vLp(1:end-1)));
taskInd=[firstResponseStrideA; firstResponseStrideP];
U=[vR-vL; vRp-vLp];
U2=[zeros(45,1); 500*ones(905,1); zeros(645,1)]; %spped profile if we didnt have the tasks
pS=[allTA.pertSize; allTP.pertSize]; %Perturbations
pS=pS(1:32);
U3=zeros(1595,1);
U3(taskInd)=pS;
%x=[15,25+[15,50,80,110,140,170,230,290,350,410,470,530,590]]; %Strides of post in which each trial was conducted

for i=1:2
pse=allTP.(vars{i});
pse=reshape(pse,13,10);
pseA=allTA.(vars{i});
pseA=reshape(pseA,19,10);
pse=[pseA;pse];

pse=pse(:,[1:3,5:end]);%Excluding subject 4, who did not respond in most post-adapt trials
y=nanmedian(pse'); %Median across subjects

pse=pse(19:end,:);
subplot(3,1,i)
%plot(pse); 
%Fit linear model of order 1 to post-adapt data only:
YY=nan(1,1595);
YY(taskInd)=y;
y=y(951:end); %Post-adapt only

t1=find(~isnan(y),1,'first');
YY2=y(t1:end);
UU=zeros(1,size(YY2,2));
opts.fixB=0;
opts.fixD=0;
opts.fixC=[]; %This could be fixed, but makes no difference
opts.Nreps=3;
opts.fastFlag=false;
opts.refineTol=1e-3;
mdl=linsys.id(dset(UU,YY2),1,opts); %No input
th(1)=mdl.initCondPrior.state;
th(2)=-1./log(mdl.A)
%LS fit to decaying data:
%th=fminunc(@(t) sum((y-t(1)*exp(-x/t(2))).^2),[200,30]);

%Trying to track PSE with additional elements:
opts.indB=1;
opts.indD=[];
opts.fixB=[];
opts.fixD=0;
opts.fixC=[]; %Could be fixed, but also we can normalize post-hoc
opts.fixX0=0;
opts.fixP0=0;
opts.Nreps=10;
opts.fastFlag=false;
UU=[U2'; U2'-U3'];
UU=U2';
t1=find(~isnan(YY),1,'first');
YY3=YY(t1:end);
UU3=UU(:,t1:end);
mdl=linsys.id(dset(UU3,YY3),1,opts);

hold on; plot(x,y,'k','LineWidth',2)
 
title(vars{i})
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
