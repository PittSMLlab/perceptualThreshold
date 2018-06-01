function fh=accPlots(trialData)
colorOff=3;
%%
B=findgroups(trialData.pertSize); %pertSize>0 means vR>vL
pp=unique(trialData.pertSize);
cmap=parula(length(pp)); %parula, jet, redbluecmap
cmap=cmap*.8;
 
fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
Q=7;

%% Define quantities of interest
nullTrials=trialData.pertSize==0;
correctResponses=trialData.initialResponse==-sign(trialData.pertSize) & ~nullTrials; %Negative response means LEFT IS SLOW (RIGHT IS FAST) choice
nonResponse=isnan(trialData.initialResponse) & ~nullTrials;
incorrectResponses=trialData.initialResponse==sign(trialData.pertSize) & ~nullTrials;

%% First: some summary stats
subplot(2,Q,1)
c=[sum(trialData.pertSize<0) sum(nullTrials) sum(trialData.pertSize>0) ]/numel(nullTrials);
bar(1,c(:,1),'FaceColor',cmap(end-colorOff+1,:))
hold on
bar(2,c(:,2),'FaceColor',cmap(round(size(cmap,1)/2),:))
bar(3,c(:,3),'FaceColor',cmap(colorOff,:))
set(gca,'XTick',1:3,'XTickLabel',{'vL<vR','NT','vL>vR'})
title('Trial decomposition')
grid on

subplot(2,Q,2)
d=[sum(trialData.initialResponse==1 & nullTrials) sum(trialData.initialResponse==-1 & nullTrials)  sum(nonResponse & nullTrials)]/sum(nullTrials);
bar(1,d(:,1),'FaceColor',cmap(end-colorOff+1,:))
hold on
bar(2,d(:,2),'FaceColor',cmap(colorOff,:))
bar(3,d(:,3),'FaceColor',cmap(round(size(cmap,1)/2),:))
hold on
set(gca,'XTick',1:3,'XTickLabel',{'<','>','NR'})
title('Null trial responses')
grid on

subplot(2,Q,3)
d=[sum(correctResponses) sum(incorrectResponses)  sum(nonResponse)]/sum(~nullTrials);
dN=[sum(nonResponse) sum(nullTrials)]/numel(nullTrials);
d2=[sum(correctResponses & (trialData.pertSize<0))/sum(trialData.pertSize<0)];
d3=sum(correctResponses & (trialData.pertSize>0))/sum(trialData.pertSize>0);
bar(1:3,d,'FaceColor',.4*ones(1,3))
hold on
b=bar([5],[d2]','FaceColor',cmap(end-colorOff+1,:));
b=bar([6],[d3]','FaceColor',cmap(colorOff,:));
set(gca,'XTick',1:6,'XTickLabel',{'OK','BAD','NR','','<0','>0'})
grid on
title('Non-null trial accuracy')


%% Second plot: % <- choices as function of perturbation size
subplot(2,Q,4:5)
S=splitapply(@(x) (sum(x==-1)+.5*sum(isnan(x)))/length(x),trialData.initialResponse,B); %Counting LEFT IS SLOW choices plus HALF of no response
scatter(pp,S,50,pp,'filled')
grid on
ylabel('% ''<-'' (left is slow) responses') 
xlabel('vL>vR     PERTURBATION (mm/s)      vL<vR')
%TODO: add error shaded area (std)
axis([-400 400 0 1])

%Add fits:
hold on
y=monoLS(S);
plot(pp,y,'k')
[z,peval,L] = fitPsycho(trialData.pertSize,trialData.initialResponse==-1 + .5*isnan(trialData.initialResponse)); %MLE fit
plot(pp,unique(peval),'r')

subplot(2,Q,6:7)
B2=findgroups(abs(trialData.pertSize)); %pertSize>0 means vR>vL
pp2=unique(abs(trialData.pertSize));
S2=splitapply(@(x) (sum(x))/length(x),correctResponses+.5*nonResponse,B2); %Counting LEFT IS SLOW choices + half of NR
S2(1)=NaN;
scatter(pp2,S2,80,.4*ones(1,3),'filled')
grid on
ylabel('% CORRECT') 
xlabel('ABS SPEED PERTURBATION (mm/s)')
axis([0 400 .5 1])
hold on 
%Adding vR>vL and vL<vR overlapping
scatter(pp(pp>0),S(pp>0),20,cmap(end,:),'filled')
scatter(abs(pp(pp<0)),1-S(pp<0),20,cmap(1,:),'filled')
%TODO: add STD

%% Third plot: compare accuracy in odd/even and first/second blocks

subplot(2,Q,Q+1) %Block number comparison
for k=1:4 %Block number
trialData1=trialData(trialData.blockNo==k & trialData.pertSize~=0,:); %Non-null trials only
correctResponses2=trialData1.initialResponse==-sign(trialData1.pertSize); %Negative response means LEFT IS SLOW (RIGHT IS FAST) choice
m=mean(correctResponses2);
s=std(correctResponses2)/sqrt(numel(correctResponses2));
bar(k,m)
hold on 
errorbar(k,m,s)
grid on
end
axis([0 5 .5 1])
xlabel('Block No')
set(gca,'XTick',1:4)
title('Block effect')
ylabel('% Correct')

subplot(2,Q,Q+2) %Individual subject x block comparison
allM=nan(4,9);
allS=nan(4,9);
for j=1:9 %Subjects
    m=nan(4,1);
    s=nan(4,1);
for k=1:4 %Block number
trialData1=trialData(trialData.blockNo==k & trialData.pertSize~=0 & trialData.subID==j,:); %Non-null trials only
correctResponses2=trialData1.initialResponse==-sign(trialData1.pertSize); %Negative response means LEFT IS SLOW (RIGHT IS FAST) choice
if ~isempty(correctResponses2)
m(k)=mean(correctResponses2);
allM(k,j)=m(k);
s(k)=std(correctResponses2)/sqrt(numel(correctResponses2));
allS(k,j)=s(k);
end
end
hold on 
%errorbar([1:4]+.1*randn(1,4),m,s,'Color','k')
plot([1:4],m,'Color','k')
grid on
end
axis([0 5 .5 1])
xlabel('Block No')
set(gca,'XTick',1:4)
title('Subject x block')
ylabel('% Correct')

subplot(2,Q,Q+3) %Individual subject comparison
clear m s
for j=1:9 %Subjects
trialData1=trialData(trialData.pertSize~=0 & trialData.subID==j,:); %Non-null trials only
correctResponses2=trialData1.initialResponse==-sign(trialData1.pertSize); %Negative response means LEFT IS SLOW (RIGHT IS FAST) choice
m(j)=mean(correctResponses2);
s(j)=std(correctResponses2)/sqrt(numel(correctResponses2));
end
hold on 
bar(1:9,m)
errorbar(1:9,m,s,'Color','k')
grid on
axis([0 10 .5 1])
xlabel('Subj. ID')
set(gca,'XTick',1:9)
ylabel('% Correct')
title('Subject effect')

%% Multi-factor analysis of responses:
%Using generalized linear model with binomial distribution and logit
%link (equivalent to multinomial logistic regression)
%Subj effect
%Block effect
%pert-size effect 
%Prev pert-size effect

%Adding prev perturbation to table:
trialData.prevSize=[0;trialData.pertSize(1:end-1)]; 
trialData.prevSize(trialData.pertSize==-250*((-1).^trialData.blockNo))=0; %Assigning NaN to previous perturbation for first trial in each block (+250 in even blocks, -250 in odd)
trialData.lastSpeedDiff=[0;trialData.lastSpeedDiff(1:end-1)];
trialData.lastSpeedDiff(trialData.pertSize==-250*((-1).^trialData.blockNo))=0; %Assigning NaN to previous perturbation for first trial in each block (+250 in even blocks, -250 in odd)

%Creating binary response variable(s):
trialData.leftResponse=trialData.initialResponse==-1;
trialData.rightResponse=trialData.initialResponse==1;
trialData.noResponse=isnan(trialData.initialResponse);

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
subplot(2,Q,4:5)
hold on
mm2.plotPartialDependence('pertSize')
text(500,-.8,removeTags(evalc('mm.disp')),'FontSize',6,'Clipping','off')
text(-500,-.4,removeTags(evalc('mm1.disp')),'FontSize',6,'Clipping','off')
text(-500,-1,regexprep(removeTags(evalc('mm2.disp')),'model:\n','model: (no subj 2)\n'),'FontSize',6,'Clipping','off')
text(-500,-1.5,regexprep(removeTags(evalc('mm3.disp')),'model:\n','model: (no subj 2)\n'),'FontSize',6,'Clipping','off')

%text(-400,-1.4,removeTags(evalc('mm.anova')),'FontSize',6,'Clipping','off')
%I'd like to do something like an anova to find out which subjects are
%significantly different from the mean, in either bias or slope, but could
%find no way to do this in Matlab. I settle for setting the 'most average'
%subject as a reference (subject 1) and comparing everyone to him/her

