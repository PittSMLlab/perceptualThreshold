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

%% Laundry list:
%From baseline data: drop non-responses(?)
%From remaining trials, estimate a distribution p(left| \Delta V)
%We'll assume that p(left|\DeltaV) really is p(left|\Delta V - PSE), only
input=[v];
obsTimes=find(isnan(input(2:end)) & ~isnan(input(1:end-1)));
allData=[allT;allTA;allTP];
startTime=allData.date*1e5+allData.startTime;
[~,idx]=sort(startTime,'ascend');
allData=allData(idx,:);
obs=discretizeObs(allData.initialResponse==1,2,[0,1]);
obsTimes=repmat(obsTimes,10,1);

N=length(input);
bias=0;
sigma=75/1.1; %Realistic based on group-level analysis of responses
range=[-500:10:500];
p=1./(1+exp((range+bias)/sigma));
pObsGivenState=[p;1-p];
O=@(u) [0;1]+[1;-1].*1./(1+exp((range+bias-u)/sigma));
%that in baseline PSE=0
%We'll assume some arbitrary transition matrices p(x_{k+1}|x_{k},u_{k}), which
s=20;
T=exp((range'-range)/(2*s^2));
pStateInitial=ones(size(range'))/numel(range);
%Run a hidden-markov model inference to get p(x_k| obs)
%Inference:
[pPredicted, pUpdated, pSmoothed] = HMMnonStationaryInferenceAlt(obs,obsTimes,input,O,T,pStateInitial);

%Compute viterbi sequence:
%[optSeq,logL]=nonStatViterbi(observations,pStateGivenPrev,pObsGivenState,pStateInitial,input,observationTimes);
%Viz:
[fh] = vizHMMInference(pSmoothed,T,O(0),obs,obsTimes,range,[0 1],1:N);

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
