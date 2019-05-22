function [fh,fh2]=accPlots(trialData)
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
subplot(5,3,1)
c=[sum(trialData.pertSize<0) sum(trialData.nullTrials) sum(trialData.pertSize>0) ]/numel(trialData.nullTrials);
bar(3,c(:,1),'FaceColor',cmap(end-colorOff+1,:),'EdgeColor','none')
hold on
bar(2,c(:,2),'FaceColor',cmap(round(size(cmap,1)/2),:),'EdgeColor','none')
bar(1,c(:,3),'FaceColor',cmap(colorOff,:),'EdgeColor','none')
set(gca,'XTick',1:3,'XTickLabel',{'vL>vR','NT','vL<vR'},'YLim',[0 .5])
title('Trial decomposition')
grid on

subplot(5,3,2)
d=[sum(trialData.initialResponse==1 & trialData.nullTrials) sum(trialData.initialResponse==-1 & trialData.nullTrials)  sum(trialData.noResponse & trialData.nullTrials)]/sum(trialData.nullTrials);
bar(1,d(:,1),'FaceColor',cmap(colorOff,:),'EdgeColor','none')
hold on
bar(2,d(:,2),'FaceColor',cmap(end-colorOff+1,:),'EdgeColor','none')
bar(3,d(:,3),'FaceColor',cmap(round(size(cmap,1)/2),:),'EdgeColor','none')
hold on
set(gca,'XTick',1:3,'XTickLabel',{'>','<','NR'})
title('Null trial responses')
grid on

subplot(5,3,3)
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


%% Fourth plot: % <- choices as function of perturbation size
subplot(5,2,[3,5])
%S=splitapply(@(x) (sum(x==-1)+.5*sum(isnan(x)))/length(x),trialData.initialResponse,B); %Counting LEFT IS SLOW choices plus HALF of no response
S=splitapply(@(x) sum(x==-1)/sum(~isnan(x)),trialData.initialResponse,B); %Not counting NR responses
scatter(pp,S,50,zeros(1,3),'filled')
grid on
ylabel('% ''<'' (left is slow) responses') 
xlabel('vL>vR     probe (mm/s)      vL<vR')
%TODO: add error shaded area (std)
axis([-400 400 0 1])

%Add fits:
hold on
%y=monoLS(S);
%plot(pp,y,'k','LineWidth',1)
%[z,peval,L] = fitPsycho(trialData.pertSize,trialData.initialResponse==-1 + .5*isnan(trialData.initialResponse)); %MLE fit
%plot(pp,unique(peval),'k','LineWidth',2)

%Add modeling:
trialData.sqrtPertSize=sign(trialData.pertSize).*sqrt(abs(trialData.pertSize));
X=trialData(~trialData.noResponse,:); %Table of trials with NR removed
%mm=fitglm(X,'leftResponse~ID*pertSize+pertSize*prevSize+blockNo*pertSize+lastSpeedDiff','Distribution','binomial'); %Logisitc regression, excluding null responses
%No subj model:
mm1=fitglm(X,'leftResponse~pertSize*prevSize+pertSize*prevFinalSpeed+blockNo*pertSize','Distribution','binomial');
mm3=fitglm(X,'leftResponse~pertSize+prevFinalSpeed-1','Distribution','binomial'); %Dropping non-sig terms
mm3.plotPartialDependence('pertSize')
text(-400,-.6,removeTags(evalc('mm1.disp')),'FontSize',6,'Clipping','off')
text(-400,-1.2,regexprep(removeTags(evalc('mm3.disp')),'model:\n','model: \n'),'FontSize',6,'Clipping','off')
Nsubs=unique(trialData.ID);
for i=1:length(Nsubs)
   mm=fitglm(X(X.ID==Nsubs(i),:),'leftResponse~pertSize+prevFinalSpeed-1','Distribution','binomial');
   mm.plotPartialDependence('pertSize')
   %y3=mm.predict(trialData(1:24,:));
   %plot(sort(trialData.pertSize(1:24)),sort(y3),'Color',.4*ones(1,3));
end
ll=findobj(gca,'Type','Line');
set(ll(1:end-1),'Color',.4*ones(1,3));
set(ll(end),'Color','k','LineWidth',2);
legend({'Data','Group model','Individuals'},'Location','SouthEast')
uistack(ll(1:end-1),'bottom')
title(['Choice vs. probe size'])
ylabel('% ''<'' (left is slow) responses') 
xlabel('vL>vR     probe (mm/s)      vL<vR')
set(gca,'XLim',[-350 350])
%% Fifth plot: accuracy

%subplot(5,Q,[3*Q+[1:Q/2], 4*Q+[1:Q/2]])%Accuracy as function of perturbation size
subplot(5,2,[4,6])
B2=findgroups(abs(trialData.pertSize)); %pertSize>0 means vR>vL
pp2=unique(abs(trialData.pertSize));
accuracy=splitapply(@(x) (sum(x))/length(x),trialData.correctResponses+.5*trialData.noResponse,B2); %Counting LEFT IS SLOW choices + half of NR
accuracy=splitapply(@(x) nansum(x)/sum(~isnan(x)),trialData.correctResponses.*abs(trialData.initialResponse),B2); %Not counting NR trials at all
accuracy(1)=NaN; %No definition of accuracy for null-trials
scatter(pp2,accuracy,50,0*ones(1,3),'filled')
grid on
ylabel('% correct') 
xlabel('(abs.) probe size (mm/s)')
axis([0 400 .5 1])
set(gca,'YTick',[.5:.05:1],'XTick',[0:50:350])
hold on 

%TODO: add STD

%Add stats:
Ntrials=splitapply(@(x) nansum(abs(x)),trialData.initialResponse,B2);
NRtrials=splitapply(@(x) sum(isnan(x)),trialData.initialResponse,B2);
p=binocdf(accuracy.*Ntrials,Ntrials,.5,'upper'); %Single tail, as we only care whether accuracy is above 50%, not below
%Add modeling:
trialData.absPertSize=abs(trialData.pertSize);
trialData.absPertDiffToPrev=abs(trialData.pertSize-trialData.prevSize);
trialData.absPertDiffToLastSS=abs(trialData.pertSize-trialData.prevFinalSpeed);
diffFun=@(x) .3.*(350./x).^.75;
trialData.sqrtAbsPertSize=1./diffFun(trialData.absPertSize);
X=trialData(~trialData.noResponse & ~trialData.nullTrials,:); %Table of trials with NR removed
mm1=fitglm(X,'correctResponses~absPertSize+absPertDiffToLastSS-1','Distribution','binomial');
mm3=fitglm(X,'correctResponses~sqrtAbsPertSize-1','Distribution','binomial'); %Dropping non-sig terms
mm2=fitglm(X,'correctResponses~absPertSize-1','Distribution','binomial'); %Dropping non-sig terms
y3=mm3.predict(trialData(1:24,:));
plot(sort(trialData(1:24,:).absPertSize),sort(y3),'k','LineWidth',2)
y3=mm2.predict(trialData(1:24,:));
%plot(sort(trialData(1:24,:).absPertSize),sort(y3),'k--','LineWidth',2)
Nsubs=unique(trialData.ID);
for i=1:length(Nsubs)
   mm=fitglm(X(X.ID==Nsubs(i),:),'correctResponses~sqrtAbsPertSize-1','Distribution','binomial');
   y3=mm.predict(trialData(1:24,:));
   plot(sort(trialData(1:24,:).absPertSize),sort(y3),'Color',.4*ones(1,3));
   %text(360,max(y3),num2str(i),'FontSize',6)
end
%mm2.plotPartialDependence('absPertSize')
mmAll=fitglm(X,'correctResponses~sqrtAbsPertSize+sqrtAbsPertSize:ID-1','Distribution','binomial') %Dropping non-sig terms
%Adding vR>vL and vL<vR overlapping
%scatter(pp(pp>0),S(pp>0),20,cmap(end,:),'filled')
%scatter(abs(pp(pp<0)),1-S(pp<0),20,cmap(1,:),'filled')


text(0,.25,removeTags(evalc('mm1.disp')),'FontSize',6,'Clipping','off')
text(0,.05,regexprep(removeTags(evalc('mm2.disp')),'model:\n','model: \n'),'FontSize',6,'Clipping','off')
text(0,-.15,regexprep(removeTags(evalc('mm3.disp')),'model:\n','model: \n'),'FontSize',6,'Clipping','off')

%title(['Accuracy vs. abs. perturbation size'])
ylabel('% correct') 
xlabel('abs. probe size (mm/s)')
set(gca,'XLim',[0 350])
title('Accuracy vs. probe size')
plot(.13*[1 1]*1050,[.5 1],'k--')
legend({'Group data','Group fit','Individual fits','Reported threshold'},'Location','SouthEast')

%% Second figure: block/learning effects & subject effects

fh2=figure('Units','Normalized','OuterPosition',[0 .5 .5 .5]);
%% Third plot: compare accuracy in odd/even and first/second blocks

subplot(1,3,1) %Block number comparison
for k=1:4 %Block number
trialData1=trialData(trialData.blockNo==k & trialData.pertSize~=0,:); %Non-null trials only
correctResponses2=trialData1.initialResponse==-sign(trialData1.pertSize); %Negative response means LEFT IS SLOW (RIGHT IS FAST) choice
m(k)=mean(correctResponses2);
s(k)=std(correctResponses2)/sqrt(numel(correctResponses2));
end
bar(1:4,m,'EdgeColor','none')
hold on 
errorbar(1:4,m,s,'Color','k','LineStyle','none','LineWidth',2)
grid on
axis([0 5 .5 1])
xlabel('Block No')
set(gca,'XTick',1:4)
title('Block effect')
ylabel('% Correct')

subplot(1,3,2) %Individual subject x block comparison
allM=nan(4,9);
allS=nan(4,9);
for j=1:9 %Subjects
    m=nan(4,1);
    s=nan(4,1);
for k=1:4 %Block number
trialData1=trialData(trialData.blockNo==k & trialData.pertSize~=0 & trialData.subID==j,:); %Non-null trials only
correctResponses2=trialData1.initialResponse==-sign(trialData1.pertSize); %Negative response means LEFT IS SLOW (RIGHT IS FAST) choice
    if ~isempty(correctResponses2) %Some subjects do not have 4 blocks
    m(k)=mean(correctResponses2);
    allM(k,j)=m(k);
    s(k)=std(correctResponses2)/sqrt(numel(correctResponses2));
    allS(k,j)=s(k);
    end
end
hold on 
%errorbar([1:4]+.1*randn(1,4),m,s,'Color','k')
plot([1:4],m,'Color','k')
text(.8,m(1)+(j-5)*.005,num2str(j),'FontSize',7)
idx=find(~isnan(m),1,'last');
text(idx+.1,m(idx)+(j-5)*.005,num2str(j),'FontSize',7)
grid on
end
axis([0 5 .5 1])
xlabel('Block No')
set(gca,'XTick',1:4)
title('Subject x block')
%ylabel('% Correct')

subplot(1,3,3) %Individual subject comparison
clear m s
for j=1:9 %Subjects
trialData1=trialData(trialData.pertSize~=0 & trialData.subID==j,:); %Non-null trials only
correctResponses2=trialData1.initialResponse==-sign(trialData1.pertSize); %Negative response means LEFT IS SLOW (RIGHT IS FAST) choice
m(j)=mean(correctResponses2);
s(j)=std(correctResponses2)/sqrt(numel(correctResponses2));
end
hold on 
bar(1:9,m,'EdgeColor','none')
errorbar(1:9,m,s,'Color','k','LineStyle','none','LineWidth',2)
grid on
axis([0 10 .5 1])
xlabel('Subj. ID')
set(gca,'XTick',1:9)
%ylabel('% Correct')
title('Subject effect')

%% Multi-factor analysis of responses:
%Using generalized linear model with binomial distribution and logit
%link (equivalent to multinomial logistic regression)
%Subj effect
%Block effect
%pert-size effect 
%Prev pert-size effect


%Creating a (modified) subject ID field: %See comment at end about ANOVA
aux=trialData.subID;
aux(aux==1)=10;
aux(aux==4)=1;
trialData.ID=categorical(aux);

%
%trialData.blockNo=categorical(trialData.blockNo);

%X=X(~isnan(trialData.initialResponse) & trialData.pertSize~=-250*((-1).^trialData.blockNo),:); %Removing null responses, and first trial of block (+250 in even blocks, -250 in odd)
X=trialData(~trialData.noResponse,:);
mm=fitglm(X,'leftResponse~ID*pertSize+pertSize*prevSize+blockNo*pertSize+lastSpeedDiff','Distribution','binomial'); %Logisitc regression, excluding null responses
%X=trialData(:,{'pertSize'});
%X.blockNo=trialData.blockNo;
%X.subID=nominal(trialData.subID);
%X.prevSize=[nan;trialData.pertSize(1:end-1)]; %Need to fix so first trial of each block has nan
%X.absSize=abs(trialData.pertSize);
%X.leftResponse=trialData.initialResponse==-1;
%X.lastSpeedDiff=[nan;trialData.lastSpeedDiff(1:end-1)];
%X=X(~isnan(trialData.initialResponse) & trialData.subID~=2 & trialData.pertSize~=-250*((-1).^trialData.blockNo),:); 

%No subj model:
mm1=fitglm(X,'leftResponse~pertSize*prevSize+lastSpeedDiff+blockNo*pertSize','Distribution','binomial');
%Reduced model removing subj 2:
mm2=fitglm(trialData(~trialData.noResponse & trialData.ID~=trialData.ID(find(trialData.subID==2,1,'first')),:),'leftResponse~pertSize*prevSize+lastSpeedDiff+blockNo*pertSize','Distribution','binomial');
mm3=fitglm(trialData(~trialData.noResponse & trialData.ID~=trialData.ID(find(trialData.subID==2,1,'first')),:),'leftResponse~pertSize+prevSize-1','Distribution','binomial'); %Dropping non-sig terms
%subplot(2,Q,4:5)
%hold on
%mm2.plotPartialDependence('pertSize')
text(500,-.8,removeTags(evalc('mm.disp')),'FontSize',6,'Clipping','off')
text(-500,-.4,removeTags(evalc('mm1.disp')),'FontSize',6,'Clipping','off')
text(-500,-1,regexprep(removeTags(evalc('mm2.disp')),'model:\n','model: (no subj 2)\n'),'FontSize',6,'Clipping','off')
text(-500,-1.5,regexprep(removeTags(evalc('mm3.disp')),'model:\n','model: (no subj 2)\n'),'FontSize',6,'Clipping','off')

%text(-400,-1.4,removeTags(evalc('mm.anova')),'FontSize',6,'Clipping','off')
%I'd like to do something like an anova to find out which subjects are
%significantly different from the mean, in either bias or slope, but could
%find no way to do this in Matlab. I settle for setting the 'most average'
%subject as a reference (subject 1) and comparing everyone to him/her

