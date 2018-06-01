function fh=accuracyPlots(trialData,goodOnly)
if nargin<2 || isempty(goodOnly)
    goodOnly=0;
end
%%
B=findgroups(trialData.pertSize); %pertSize>0 means vR>vL
pp=unique(trialData.pertSize);
cmap=parula(length(pp)); %parula, jet, redbluecmap
cmap=cmap*.8;
 
fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
Q=6;
%%
nullTrials=trialData.pertSize==0;
correctResponses=trialData.initialResponse==-sign(trialData.pertSize) & ~nullTrials; %Negative response means LEFT IS SLOW (RIGHT IS FAST) choice
nonResponse=isnan(trialData.initialResponse) & ~nullTrials;
incorrectResponses=trialData.initialResponse==sign(trialData.pertSize) & ~nullTrials;
B2=findgroups(abs(trialData.pertSize)); %pertSize>0 means vR>vL
pp2=unique(abs(trialData.pertSize));
S2=splitapply(@(x) (sum(x))/length(x),correctResponses+.5*nonResponse,B2); %Counting LEFT IS SLOW choices + half of NR
S2(1)=NaN;
colorOff=3;
%% Third plot: reaction times
rt=trialData.reactionTime;
%rt=trialData.reactionStride;
%fun=@nanmean;
fun=@nanmedian;
%subplot(2,Q,1)
%scatter(trialData.pertSize,rt,10,zeros(1,3))
%hold on
%scatter(pp,S3,50,pp)
S4=splitapply(fun,rt,B);
% scatter(pp,S4,70,pp,'filled')
% grid on
% ylabel('Reaction log-time (s))') 
% xlabel('vL>vR      PERTURBATION         vL<vR')
% legend({'Indiv. trials','Median'},'Location','NorthWest')
% set(gca,'YScale','log')

subplot(2,Q,1:2)
S=splitapply(fun,rt,B); 
S2=splitapply(fun,rt,B2); 
scatter(pp2,S2,80,.4*ones(1,3),'filled')
grid on
ylabel('Reaction log-time (s)') 
xlabel('ABS SPEED PERTURBATION (mm/s)')
hold on 
%Adding vR>vL and vL<vR overlapping
scatter(pp(pp>0),S(pp>0),20,cmap(end,:),'filled')
scatter(abs(pp(pp<0)),S(pp<0),20,cmap(1,:),'filled')
%TODO: add STD
set(gca,'YScale','log')

%% Fifth: Reaction time vs. accuracy
subplot(2,Q,3:4)
M=min(20,round(size(trialData,1)/10));
S=splitapply(@nanmedian,trialData.reactionTime,B);
T=splitapply(@(x) (sum(x))/length(x),correctResponses,B);
T(12)=NaN;
[x,idx]=sort(trialData.reactionTime,'ascend');
y=correctResponses(idx);
y2=monoLS(y(~isnan(x)));
plot(x(~isnan(x)),y2,'k')
axis([.5 10 .4 1])
y=conv(y,ones(M,1)/M,'same');
y(1:floor(M/2))=NaN;
y(end-ceil(M/2):end)=NaN;
%plot(x,y,'k')
hold on
scatter(S,T,70,pp,'filled')
grid on
set(gca,'XScale','log')
xlabel('Reaction log-time (s)')
ylabel('% CORRECT')

%subplot(2,Q,4)
S=splitapply(@nanmedian,trialData.reactionTime,B2);
T=splitapply(@(x) (sum(x))/length(x),correctResponses,B2);
T(1)=NaN;
%Postivie:
rt=trialData.reactionTime(trialData.pertSize>0);
[x,idx]=sort(rt,'ascend');
y=correctResponses(trialData.pertSize>0);
y=y(idx);
y2=monoLS(y(~isnan(x)));
y=conv(y,ones(M,1)/M,'same');
y(1:floor(M/2))=NaN;
y(end-ceil(M/2):end)=NaN;
%plot(x,y,'Color',cmap(end,:))
%plot(x(~isnan(x)),y2,'Color',cmap(end,:))
hold on
%Negative:
rt=trialData.reactionTime(trialData.pertSize<0);
[x,idx]=sort(rt,'ascend');
y=correctResponses(trialData.pertSize<0);
y=y(idx);
y2=monoLS(y(~isnan(x)));
y=conv(y,ones(M,1)/M,'same');
y(1:floor(M/2))=NaN;
y(end-ceil(M/2):end)=NaN;
%plot(x,y,'Color',cmap(1,:))
%plot(x(~isnan(x)),y2,'Color',cmap(1,:))
%All:
scatter(S,T,70,.4*ones(1,3),'filled')
grid on
set(gca,'XScale','log')
xlabel('Reaction log-time (s)')
ylabel('% CORRECT')
hold on
[x,i]=sort(trialData.reactionTime);
y=correctResponses(i);
y2=monoLS(y(~isnan(x)));
%plot(x(~isnan(x)),y2,'k')
axis([.5 10 .4 1])

%%
%Adding prev perturbation to table:
trialData.prevSize=[0;trialData.pertSize(1:end-1)]; 
trialData.prevSize(trialData.pertSize==-250*((-1).^trialData.blockNo))=0; %Assigning NaN to previous perturbation for first trial in each block (+250 in even blocks, -250 in odd)
trialData.lastSpeedDiff=[0;trialData.lastSpeedDiff(1:end-1)];
trialData.lastSpeedDiff(trialData.pertSize==-250*((-1).^trialData.blockNo))=0; %Assigning NaN to previous perturbation for first trial in each block (+250 in even blocks, -250 in odd)

%Creating binary response variable(s):
trialData.leftResponse=trialData.initialResponse==-1;
trialData.rightResponse=trialData.initialResponse==1;
trialData.noResponse=isnan(trialData.initialResponse);

%
trialData.invAbsPertSize=1./abs(trialData.pertSize);
trialData.absPertSize=abs(trialData.pertSize);
trialData.pertSign=sign(trialData.pertSize);
%
trialData.correctResponse=trialData.initialResponse==-sign(trialData.pertSize) & ~nullTrials; 
%
aux=trialData.subID;
aux(aux==1)=10;
aux(aux==9)=1;
trialData.ID=categorical(aux);
%% Modeling rt as a function of abs(pertSize), sign, etc
X=trialData(~trialData.noResponse & ~nullTrials,:);
mm=fitglm(X,'reactionTime~invAbsPertSize*pertSign+blockNo+correctResponse*invAbsPertSize','Distribution','poisson','Link','log','DispersionFlag',true); %Logisitc regression, excluding null responses
mm1=fitglm(X,'reactionTime~invAbsPertSize*pertSign+blockNo+correctResponse*invAbsPertSize','Distribution','poisson','Link','identity','DispersionFlag',true); %Logisitc regression, excluding null responses
mm2=fitglm(X,'reactionTime~absPertSize*pertSign+blockNo+correctResponse*absPertSize','Distribution','poisson','Link','reciprocal','DispersionFlag',true); %Logisitc regression, excluding null responses
mm3=fitglm(X,'reactionTime~absPertSize+correctResponse:absPertSize','Distribution','poisson','Link','reciprocal','DispersionFlag',true); %Logisitc regression, excluding null responses
mmAll=fitglm(X,'reactionTime~absPertSize*pertSign+blockNo+ID*absPertSize*correctResponse','Distribution','poisson','Link','reciprocal','DispersionFlag',true); %Logisitc regression, excluding null responses
mmAll2=fitglm(X,'reactionTime~absPertSize+correctResponse:absPertSize+ID:absPertSize','Distribution','poisson','Link','reciprocal','DispersionFlag',true); %Logisitc regression, excluding null responses
%text(20/2e4,.1,removeTags(evalc('mm.disp')),'FontSize',9,'Clipping','off')
text(40/1e2,.1,removeTags(evalc('mm2.disp')),'FontSize',9,'Clipping','off')
text(40/1e2,-.4,removeTags(evalc('mm3.disp')),'FontSize',9,'Clipping','off')
text(40/1e4,.1,removeTags(evalc('mm1.disp')),'FontSize',9,'Clipping','off')
text(20,.55,removeTags(evalc('mmAll.disp')),'FontSize',9,'Clipping','off')
text(20,-.35,removeTags(evalc('mmAll2.disp')),'FontSize',9,'Clipping','off')
end

