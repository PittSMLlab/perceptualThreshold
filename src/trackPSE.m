%%
%addpath(genpath('../src/datlogManipulation/'))
addpath('../../ext/altmany-export_fig-b1a7288/')
%%
clear all
 addpath(genpath('../'))
dataDir='../data/';
subList=dir([dataDir 'AB*']);
addpath(genpath('../../matlab-linsys/'))
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
             bs.blockNo=9*ones(size(bs,1),1);
        bs.subID=j*ones(size(bs,1),1);
    if j>1
        allT1=[allT1;bs];
    else
        allT1=bs;
    end
end
%% Get adapt trials
for j=sL
    td=readtable([dataDir subList(j).name '/Adaptation.csv']);
    td.blockNo=3*ones(size(td,1),1);
    td.subID=j*ones(size(td,1),1);
    if j>1
       allTA=[allTA;td];
    else
        allTA=td;
    end
end
%% Get post trials
for j=sL
    td=readtable([dataDir subList(j).name '/PostAdaptation.csv']);
    td.blockNo=4*ones(size(td,1),1);
    td.subID=j*ones(size(td,1),1);
    if j>1
       allTP=[allTP;td];
    else
        allTP=td;
    end
end

%% Speed profiles:
load ../data/dynamicProfiles.mat
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
breakSize=0;
sst=[sstb; sstb+length(vLb)+breakSize; ssta+2*length(vLb)+2*breakSize; sstp+2*length(vLb)+length(vL)+3*breakSize];
eet=[eetb; eetb+length(vLb)+breakSize; eeta+2*length(vLb)+2*breakSize; eetp+2*length(vLb)+length(vL)+3*breakSize];
pps=[pS; pS; pSa; pSp];

v=[vRb-vLb; nan(breakSize,1); vRb-vLb; nan(breakSize,1); vR-vL; nan(breakSize,1); vRp-vLp]; %This will be the input

%% Some pre processing
allData=[allT;allTA;allTP];
startTime=allData.date*1e5+allData.startTime;
[~,idx]=sort(startTime,'ascend');
allData=allData(idx,:);
baseData=[allT;allT1];
baseData.leftResponse=baseData.initialResponse==-1;
trueResponses=baseData(~isnan(baseData.reactionTime),:);
psychom=fitglm(trueResponses,'leftResponse~pertSize','Distribution','binomial')
sigma=1/psychom.Coefficients.Estimate(2) %Empirically estimated sigma for accuracy response
bias=-psychom.Coefficients.Estimate(1)/psychom.Coefficients.Estimate(2)
baseData=[allT];
baseData.leftResponse=baseData.initialResponse==-1;
trueResponses=baseData(~isnan(baseData.reactionTime),:);
rsp=fitglm(trueResponses,'lastSpeedDiff~pertSize')
err=trueResponses.lastSpeedDiff-(rsp.Coefficients.Estimate(1)+rsp.Coefficients.Estimate(2)*trueResponses.pertSize);
sig=std(err)
bis=rsp.Coefficients.Estimate(2)


%%
%From baseline data: drop non-responses(?)
%From remaining trials, estimate a distribution p(left| \Delta V)
%We'll assume that p(left|\DeltaV) really is p(left|\Delta V - PSE), only
input=[v];
range=[-500:2:500];
%range=[-300:6:300];
for dataUsed=2%1:2

switch dataUsed
    case 1
        obs=discretizeObs(allData.initialResponse==1,2,[0,1]);
        bias=0;
        %sigma=75/1.1; %Realistic based on group-level analysis of responses. These numbers are from the static task. Do they match the dynamic one?
        p=1./(1+exp((range+bias)/sigma));
        pObsGivenState=[p;1-p];
        O=@(u) [0;1]+[1;-1].*1./(1+exp((range+bias-u)/sigma));
        oRange=[0 1];
    case 2
        obs=discretizeObs(allData.lastSpeedDiff,length(range),[min(range) max(range)]); %Same range
        O=@(u) exp(-(range'-range).^2/(2*sig^2)); %Presuming that lastSpeedDiff is an unbiased estimator of PSE
        O=@(u) exp(-(range'-range-bis*(u-range)).^2/(2*sig^2)); % 35% undercorrection, derived from base 
        oRange=range;
end

obsTimes=find(isnan(input(2:end)) & ~isnan(input(1:end-1)));
obsTimes=repmat(obsTimes,10,1);
%Remove nan responses:
obs=obs(~isnan(allData.reactionTime));
obsTimes=obsTimes(~isnan(allData.reactionTime));
perSize=allData.pertSize(~isnan(allData.reactionTime));


N=length(input);
%We'll assume some arbitrary transition matrices p(x_{k+1}|x_{k},u_{k}), which
s=7;
T=@(u) exp(-(range'-range).^2./(2*s^2));
pStateInitial=ones(size(range'))/numel(range);
N=length(input);
%Run a hidden-markov model inference to get p(x_k| obs)
%Inference:
[pPredicted, pUpdated, pSmoothed] = HMMnonStationaryInferenceAlt(obs,obsTimes,input,O,T,pStateInitial);

%Compute viterbi sequence:
%[optSeq,logL]=nonStatViterbi(obs,T,O,pStateInitial,input,obsTimes);
%Viz:
[fh] = vizHMMInference(pSmoothed,T(0),O(0),obs,obsTimes,range,oRange,1:N);
%[fh] = vizHMMInference(pUpdated,T(0),O(0),obs,obsTimes,range,oRange,1:N);
fh.Units='Pixels';
fh.InnerPosition=[500 300 300*2 2*300];
%Add adapt/post separation\
ph=findobj(fh,'Type','Axes');
axes(ph(2))
hold on
plot(888*[1 1],500*[1 -1],'k')
%plot(858*[1 1],500*[1 -1],'Color',.5*ones(1,3)) %Last measurement of base
%plot(903*[1 1],500*[1 -1],'Color',.5*ones(1,3)) %First of A
plot(1793*[1 1],500*[1 -1],'k')
%plot(1763*[1 1],500*[1 -1],'Color',.5*ones(1,3)) %Last measurement of base
%plot(1808*[1 1],500*[1 -1],'Color',.5*ones(1,3)) %First of A
[~,ml]=max(pSmoothed);
plot(1:N,range(ml),'k') %Ml states
set(gca,'XTick',sort(unique(obsTimes)),'XTickLabel',{},'YTick',[-450:150:450])
grid on
if dataUsed==2 %Add viterbi sequence
%plot(1:N,range(optSeq),'r')
end
axes(ph(1))
hold off
fg=findgroups(obsTimes);
aux=splitapply(@(x) nanmean(x),oRange(obs)',fg);
aux2=splitapply(@(x) nanmean(x),obsTimes,fg);
ps=splitapply(@(x) nanmean(x),perSize,fg);
ss=scatter(aux2(ps==400),aux(ps==400),'filled','DisplayName','\Delta V=400 mm/s');
   ss.CData=.7*ss.CData;
   hold on
   scatter(aux2(ps==200),aux(ps==200),'filled','DisplayName','\Delta V=200 mm/s','CData',ss.CData*.5 +.5*[1 1 1]);
   scatter(aux2(ps==100),aux(ps==100),'filled','DisplayName','\Delta V=100 mm/s','CData',ss.CData*.15 +.85*[1 1 1],'MarkerEdgeColor',ss.CData)
   ph(1).Position(3)=ph(2).Position(3);
   ph(1).XAxis.Limits=ph(2).XAxis.Limits;
   title('observations')
      legend('Location','NorthWest','AutoUpdate','off')
if dataUsed==1
   ph(3).YAxis.Label.String='choice';
   ph(3).YAxis.TickValues=[0 1];
   ph(1).YAxis.Limits=[0 1];
   ph(3).YAxis.TickLabels={'right','left'};
   ylabel({'% of right'; 'side choices'})
   axes(ph(1))
  rr=rectangle('Position',[885 0 900 1],'LineWidth',2);
  uistack(rr,'bottom')
else
   ph(3).YAxis.Label.String='reported PSE (mm/s)';
   ph(1).YAxis.Label.String={'avg. reported'; 'PSE (mm/s)'};
   ph(1).YAxis.Limits=[-50 450];
   %for k=[400,200,100]; %Add errorbars
   % errorbar(aux2(ps==k),aux(ps==k),);
   %end
   axes(ph(1))
   rr=rectangle('Position',[885 -50 900 500],'LineWidth',2);
   uistack(rr,'bottom')
end
plot(888*[1 1],500*[1 -1],'k')
%plot(858*[1 1],500*[1 -1],'Color',.5*ones(1,3)) %Last measurement of base
%plot(903*[1 1],500*[1 -1],'Color',.5*ones(1,3)) %First of A
plot(1793*[1 1],500*[1 -1],'k')
ph(2).YAxis.Limits=[-100 450];
ph(1).XAxis.TickValues=ph(2).XAxis.TickValues;
ph(1).XTickLabels={};
ph(3).Title.String='emission probabilities';
ph(4).Title.String='transition probabilities';
ph(3).XAxis.Label.String='PSE - \Delta V (mm/s)';
ph(4).XAxis.Label.String='current PSE (mm/s)';
ph(4).YAxis.Label.String='next PSE (mm/s)';
ph(2).YAxis.Label.String={'estimated';'PSE (mm/s)'};
ph(1).XAxis.Label.String='strides';
ph(2).Title.String='PSE estimates (probability density)';
ll=findobj(ph(2),'Type','Line','Color',[0 0 0]);
set(ll(end),'DisplayName','MLE')
axes(ph(2))
legend(ll(end))
ph(1).Title.String='avg. observation';
ph(1).Position(1)=.12;
ph(1).Position(4)=.22;
ph(2).Position(4)=.22;
ph(2).Position(1)=.12;
ph(3).Position(1)=.6;
ph(1).Position(3)=ph(2).Position(3);
vizRange=300*[-1 1];
ph(4).XAxis.Limits=vizRange;
ph(3).XAxis.Limits=vizRange;
ph(4).YAxis.Limits=vizRange;
ph(4).CLim=[0 .15];
ph(1).Position(2)=.05;
for jj=1:4
    ph(jj).Title.FontWeight='normal';
end
axes(ph(2))
rectangle('Position',[885 -100 900 550],'LineWidth',2)
set(ph,'FontSize',10,'FontName','OpenSans')
%set(ph,'CLim',[0 .3]) 
if dataUsed==1
    export_fig ../fig/trackPSE_2AFC.png -png -c[0 5 0 5] -transparent -r600
else
    export_fig ../fig/trackPSE_speed.png -png -c[0 5 0 5] -transparent -r600
end
end

%% The real deal:
%Get table with responses, generate observation vectors

%Generate appropriate obs and transition matrices

%Estimate optimal probabilities


%% Find temporal decay during post-adapt

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
