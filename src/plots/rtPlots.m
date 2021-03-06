function [fh,f2]=rtPlots(trialData,goodOnly)
if nargin<2 || isempty(goodOnly)
    goodOnly=0;
end


%%
%Creating binary response variable(s):
trialData.leftResponse=trialData.initialResponse==-1;
trialData.rightResponse=trialData.initialResponse==1;
trialData.noResponse=isnan(trialData.initialResponse);

%
trialData.absPertSize=abs(trialData.pertSize);
trialData.pertSign=sign(trialData.pertSize);

%
nullTrials=trialData.pertSize==0;
trialData.correctResponse=trialData.initialResponse==-sign(trialData.pertSize) & ~nullTrials; 
%
aux=trialData.subID;
aux(aux==1)=10;
aux(aux==9)=1;% To use subject 9 as reference
aux(aux==10)=9;
trialData.ID=categorical(aux);

%Define some aux vars:
correctResponses=trialData.initialResponse==-sign(trialData.pertSize) & ~nullTrials; %Negative response means LEFT IS SLOW (RIGHT IS FAST) choice
incorrectResponses=trialData.initialResponse==sign(trialData.pertSize) & ~nullTrials;

B2=findgroups(abs(trialData.pertSize)); %pertSize>0 means vR>vL
pp2=unique(abs(trialData.pertSize));
Q=4;
%% reaction times vs. pert size
rt=trialData.reactionTime;
fun=@nanmean;
%fun=@nanmedian;
fh=figure('Units','Normalized','OuterPosition',[.5 0 .5 1]);

subplot(2,Q,1:2)
RT=splitapply(fun,rt,B2); %Sort by absolute pert size
eRT=splitapply(@(x) nanstd(x)/sqrt(sum(~isnan(x))),rt,B2); 
ss=scatter(pp2,RT,50,.4*ones(1,3),'filled','DisplayName','All data');
hold on

grid on
ylabel('Reaction time [RT] (s)') 
xlabel('Abs. speed difference (mm/s)')
hold on 
%Adding vR>vL and vL<vR overlapping
%scatter(pp(pp>0),S(pp>0),20,cmap(end,:),'filled')
%scatter(abs(pp(pp<0)),S(pp<0),20,cmap(1,:),'filled')
set(gca,'XLim',[-10 380])
%set(gca,'YScale','log')

%Add" reaction times for correct vs. incorrect responses separately:
rtC=rt(correctResponses);
BC=B2(correctResponses)-1; %Subtracting 1 to eliminate the first group, which corresponds to null trials. splitapply fails for empty groups
S2=splitapply(fun,rtC,BC); %Sort by absolute pert size
E2=splitapply(@(x) std(x)/sqrt(length(x)),rtC,BC); 
%scatter(pp2(2:end)+10,S2,20,zeros(1,3),'filled')
%errorbar(pp2(2:end)*1.1,S2,E2,'k') 
rtI=rt(incorrectResponses);
BI=B2(incorrectResponses)-1;
S2=splitapply(fun,rtI,BI); %Sort by absolute pert size
E2=splitapply(@(x) std(x)/sqrt(length(x)),rtI,BI); 
%scatter(pp2(2:end)+10,S2,20,zeros(1,3))
%errorbar(pp2(2:end)*1.05,S2,E2,'k')  %Very large bars

% Add modeling:


%These models allow for a time offset while still being
%inversely proportional to perturbation size. However, the correct-based
%offset comes from a few subjects (2 and 3) who are both bad at the task
%and slow to respond.
trialData.invAbsPertSize=1./max(trialData.absPertSize,35);
X=trialData(~trialData.noResponse & ~nullTrials ,:);
mm2=fitglm(X,'reactionTime~invAbsPertSize*pertSign+blockNo','Distribution','poisson','Link','identity','DispersionFlag',true); %Logisitc regression, excluding null responses
mm3=fitglm(X,'reactionTime~invAbsPertSize','Distribution','poisson','Link','identity','DispersionFlag',true); %Logisitc regression, excluding null responses

%These models ignored the correctness of responses:
%mm2=fitglm(X,'reactionTime~absPertSize*pertSign+blockNo','Distribution','poisson','Link','reciprocal','DispersionFlag',true); %Logisitc regression, excluding null responses
%mm3=fitglm(X,'reactionTime~absPertSize','Distribution','poisson','Link','reciprocal','DispersionFlag',true); %Logisitc regression, excluding null responses
yp=mm3.predict(X(X.correctResponse,:));
xc=X.absPertSize(X.correctResponse,:);
[~,idx]=sort(xc);
pm=plot(xc(idx),yp(idx),'k','LineWidth',2,'DisplayName','Model fit (group)');
uistack(pm,'bottom')
%yp=mm3.predict(X(~X.correctResponse,:));
%xc=X.absPertSize(~X.correctResponse,:);
%[~,idx]=sort(xc);
%plot(xc(idx),yp(idx),'k--','LineWidth',2)

%Add individual models and plots:
Nsubs=unique(trialData.ID);
for i=1:length(Nsubs)
   X2=X(X.ID==Nsubs(i),:);
   %mm6=fitglm(X2,'reactionTime~absPertSize','Distribution','poisson','Link','reciprocal','DispersionFlag',true); %Logisitc regression, excluding null responses
   mm6=fitglm(X2,'reactionTime~invAbsPertSize','Distribution','poisson','Link','identity','DispersionFlag',true); %Logisitc regression, excluding null responses
   y3=mm6.predict(trialData(1:24,:));
   ph=plot(sort(trialData.absPertSize(1:24),'ascend'),sort(y3,'descend'),'Color',.7*ones(1,3),'DisplayName','Individual fits');
   text(370,min(y3),num2str(i),'FontSize',6)
   uistack(ph,'bottom')
%   S7=splitapply(@nanmedian,rt(trialData.ID==Nsubs(i)),B2(trialData.ID==Nsubs(i))); %Sort by absolute pert size
%scatter(pp2,S7,10,.7*ones(1,3),'filled')
%plot(pp2,S7,'Color',.7*ones(1,3))
end


e=errorbar(pp2,RT,eRT,'k');
e.LineStyle='none';
legend([ss,pm,ph])
%mm2=fitglm(X,'reactionTime~absPertSize*pertSign+blockNo+correctResponse*absPertSize','Distribution','poisson','Link','reciprocal','DispersionFlag',true); %Logisitc regression, excluding null responses
%mm3=fitglm(X,'reactionTime~absPertSize+correctResponse:absPertSize','Distribution','poisson','Link','reciprocal','DispersionFlag',true); %Logisitc regression, excluding null responses
%mm4=fitglm(X,'reactionTime~invAbsPertSize+correctResponse*invAbsPertSize','Distribution','poisson','Link','identity','DispersionFlag',true); %Logisitc regression, excluding null responses
text(400,8,removeTags(evalc('mm2.disp')),'FontSize',7,'Clipping','off')
text(400,1,removeTags(evalc('mm3.disp')),'FontSize',7,'Clipping','off')
%text(400,.6,removeTags(evalc('mm4.disp')),'FontSize',7,'Clipping','off')
%TO INCLUDE: prediction of accuracy can be improved bvy considering
%response time, possibly because it measures trials in which subjects were
%distracted, which increases RT and reduces accuracy:
%mm=fitglm(X,'correctResponse~absPertSize+reactionTime:absPertSize','Link','logit','Distribution','binomial');
%predictedCorrect=mm.predict(X)>.5;
%accuracy=sum(predictedCorrect==X.correctResponse)/length(predictedCorrect)
%mm2=fitglm(X,'correctResponse~absPertSize','Link','logit','Distribution','binomial');
%predictedCorrect=mm2.predict(X)>.5;
%accuracy=sum(predictedCorrect==X.correctResponse)/length(predictedCorrect)

%% Modeling individuals:
X=trialData(~trialData.noResponse & ~nullTrials,:);
mmAll2=fitglm(X,'reactionTime~invAbsPertSize*ID','Distribution','poisson','Link','identity','DispersionFlag',true) %Logisitc regression, excluding null responses
%mmAll2=fitglm(X,'reactionTime~absPertSize*ID','Distribution','poisson','Link','reciprocal','DispersionFlag',true) %Logisitc regression, excluding null responses
%text(20,.55,removeTags(evalc('mmAll.disp')),'FontSize',9,'Clipping','off')
%text(20,-.35,removeTags(evalc('mmAll2.disp')),'FontSize',9,'Clipping','off')

%% Add individual mean RT vs. mean acc:
subplot(2,Q,5:6)
%acc vs. RT for each subject & pertSize
Nsubs=unique(trialData.ID);
fun=@nanmean;
B2=findgroups(X.absPertSize);
pp=unique(X.absPertSize);
hold on
for i=1:length(Nsubs)
    mrt(i)=nanmean(X.reactionTime(X.ID==Nsubs(i)));
    N(i)=sum(X.ID==Nsubs(i));
    ert(i)=nanstd(X.reactionTime(X.ID==Nsubs(i)))/sqrt(N(i));
    macc(i)=nanmean(X.correctResponse(X.ID==Nsubs(i)));
   text(mrt(i)+.1,macc(i)-.01,num2str(i),'FontSize',6)
end
scatter(mrt,macc,50,.4*ones(1,3),'filled')
errorbar(mrt,macc,macc.*(1-macc)./sqrt(N),'Vertical','Color','k','LineStyle','none')
errorbar(mrt,macc,ert,'Horizontal','Color','k','LineStyle','none')
xlabel('Mean RT (s)')
ylabel('Mean accuracy (%)')
title('RT vs accuracy across subjects')
text(6,.8,removeTags(evalc('mmAll2.disp')),'FontSize',7,'Clipping','off')

%% Add some comparison to EZ-DD model (group level):
trialData.nullTrials=trialData.pertSize==0;
trialData.correctResponses=trialData.initialResponse==-sign(trialData.pertSize) & ~trialData.nullTrials;
%Compute accuracy: and RT as function of pertSize
rt=trialData.reactionTime;
fun=@nanmean;
%fun=@nanmedian;
f2=figure('Units','Normalized','OuterPosition',[.5 0 .5 .3]);


B2=findgroups(abs(trialData.pertSize)); %pertSize>0 means vR>vL
pp2=unique(abs(trialData.pertSize));
Acc=splitapply(@(x) nansum(x)/sum(~isnan(x)),trialData.correctResponses.*abs(trialData.initialResponse),B2); %Not counting NR trials at all
Acc(1)=.5; %No definition of accuracy for null-trials
eAcc=sqrt(Acc.*(1-Acc))./sqrt(splitapply(@(x) sum(~isnan(x)),trialData.correctResponses.*abs(trialData.initialResponse),B2));
empDiff= -1./log(Acc.^-1 -1);%Estimate of s^/av
rt=rt(trialData.correctResponses | trialData.nullTrials);
B2=B2(trialData.correctResponses | trialData.nullTrials);
RT=splitapply(fun,rt,B2); %Sort by absolute pert size
vRT=splitapply(@(x) nanvar(x),rt,B2);
eRT=splitapply(@(x) nanstd(x)/sqrt(sum(~isnan(x))),rt,B2);

%regardless of how we map stimulus intensity to drift rate!
diff=[.01:.01:5]; %difficulty values= s^2/av in Wagenmakers EZ model notation
Toff=0;
s=[.28];
for i=1:length(s)
    y=-1./diff; %Inverse difficulty
    a=1; %Can be arbitrarily defined (only s^2/a^2 matters)
    a_v= a^2/s(i)^2 .*diff;%a/v
    %v=s(i)^2./(a.*diff);
    Pc=1./(1+exp(y)); 
    MDT=Toff+(.5*a_v).*(1-exp(y))./(1+exp(y)); 
    %VRT=.5*diff.*(a./v).^2 .*(2.*y.*exp(y)-exp(2*y)+1).*Pc.^2; %Equivalent
    VRT=.5*(a_v).^2 .*(diff.*(2*Pc-1)-2*Pc.*(1-Pc));
    subplot(1,3,2); hold on; plot(diff,MDT,'DisplayName',['s=' num2str(s(i))],'LineWidth',2); 
    subplot(1,3,1); hold on; plot(diff,Pc,'DisplayName',['s=' num2str(s(i))],'LineWidth',2); 
    subplot(1,3,3); hold on; plot(MDT,Pc,'LineWidth',2);
    %subplot(2,4,8); hold on; plot(VRT,Pc,'LineWidth',2);
    %subplot(2,2,3); hold on; plot(diff,VRT,'LineWidth',2);
end
diffFun=@(x) .3.*(350./x).^.75;
subplot(1,3,2);
xlabel('Difficulty'); 
ylabel('Mean Reaction Time (MRT)')
cc=pp2;
cc=.4*ones(1,3);
scatter(diffFun(pp2),RT,50,cc,'filled')
errorbar(diffFun(pp2),RT,eRT,'Vertical','Color','k','LineStyle','none')
%plot([0 5],(Toff+m(i)/2) *[1 1],'k--')
%legend({'EZ model fit','Group data','ste','MDT(Diff=\infty)'},'Location','SouthEast'); 
subplot(1,3,1);
xlabel('Difficulty'); 
ylabel('p(correct)')
scatter(diffFun(pp2),Acc,50,cc,'filled')
errorbar(diffFun(pp2),Acc,eAcc,'Vertical','Color','k','LineStyle','none')
legend({'DD model fit','Group data','ste'},'Location','NorthEast'); 
%subplot(2,2,3)
%scatter(diffFun(pp2),vRT,50,cc,'filled')
%xlabel('Difficulty')
%ylabel('Variance of RT')
subplot(1,3,3);
scatter(RT,Acc,50,cc,'filled')
errorbar(RT,Acc,eRT,'Horizontal','Color','k','LineStyle','none')
errorbar(RT,Acc,eAcc,'Vertical','Color','k','LineStyle','none')
%legend({'EZ model fit','Data','ste'},'Location','SouthWest')
xlabel('MRT')
ylabel('p(correct)')
%subplot(2,4,8)
%scatter(vRT,Acc,50,cc,'filled')
%xlabel('VRT')
%ylabel('p(correct)')

%% Group fits do not work well, but individual fits may:
%Step 1: fit a psychometric curve to each subject
%Step 2: use the param fits to infer difficulty for each belt speed difference
%Step 3: use VRT, RT from all correct trials somehow to infer the best
%value for s^2 for that subject
%Idea: express VRT, MDT as a function of difficulty and (a/s). [CAN BE
%DONE], use empirical Accuracy & difficulty (or psychometric accuracy and
%difficulty?) to infer optimal value of a/s given the subject's VRT -> Does
%this work? Too little data per subject!
end

