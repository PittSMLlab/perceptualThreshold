function [fh]=accPlots(trialData)

%% Define quantities of interest
%nullTrials=trialData.pertSize==0;
%correctResponses=trialData.initialResponse==-sign(trialData.pertSize) & ~nullTrials; %Negative response means LEFT IS SLOW (RIGHT IS FAST) choice
%nonResponse=isnan(trialData.initialResponse) & ~nullTrials;
%incorrectResponses=trialData.initialResponse==sign(trialData.pertSize) & ~nullTrials;

%Adding prev perturbation to table:
trialData.prevSize=[0;trialData.pertSize(1:end-1)]; 
trialData.isFirstInBlock=[1;diff(trialData.blockNo)~=0].*sign(trialData.pertSize);
trialData.prevSize(trialData.isFirstInBlock~=0)=0; %Assigning NaN to previous perturbation for first trial in each block
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

%% Remove 250 if needed (and -250 for balancing):
trialData=trialData(trialData.isFirstInBlock==0,:); %Test, remove the first trial in each block
%% Get probe sizes
B=findgroups(trialData.pertSize); %pertSize>0 means vR>vL
pp=unique(trialData.pertSize);
%% First figure: global stats
fh=figure('Units','Normalized','OuterPosition',[.4 .65 .6 .35]);
Q=6;

%% First plot: choices as function of probe size
subplot(1,3,1)
S=splitapply(@(x) sum(x==-1)/sum(~isnan(x)),trialData.initialResponse,B); %Not counting NR responses
E=splitapply(@(x) nanstd(x==-1)/sqrt(sum(~isnan(x))),trialData.initialResponse,B); %Not counting NR responses
ss=scatter(pp,S,30,zeros(1,3),'filled','MarkerEdgeColor','w');
grid on
ylabel('% ''<'' (left is slow) responses') 
xlabel('vL>vR     probe (mm/s)      vL<vR')
axis([-360 360 0 1])
X=trialData;
%Add fits:
hold on
errorbar(pp,S,E,'k','LineStyle','none')
X.pertSign=sign(X.pertSize);
%X.pertSize=100*sign(X.pertSize).*abs(X.pertSize/100).^.6;
%Attempt at fitting with a given alpha:
%alpha=.6;
%LL.Link=@(mu) log(mu./(1-mu)).^(1/alpha); %Need to change all this to consider powers of negative numbers
%LL.Inverse=@(y) exp(y.^alpha)./(1+exp(y.^alpha));
%LL.Derivative=@(mu) (1/alpha).* log(mu./(1-mu)).^(1/alpha -1) .* (1./(mu.*(1-mu)));
%simple model:
%mm0=fitglm(X,'leftResponse~pertSize','Distribution','binomial')
%mm0=fitglm(X,'leftResponse~pertSize-1','Distribution','binomial')
%mm0.plotPartialDependence('pertSize')
%Full monty:
%mm0=fitglm(X,'leftResponse~ID+pertSize*ID+pertSize:blockNo+pertSize:blockNo:ID+pertSize:pertSign+pertSize:pertSign:ID+prevSize*ID','Distribution','binomial')
%mm0=mm0.step('Upper','leftResponse~ID+pertSize*ID+pertSize:blockNo*ID+pertSize:pertSign*ID+prevSize*ID','Criterion','Deviance','PEnter',0,'PRemove',0.05,'Nsteps',Inf)
%Full monty w/o ID:
frml='leftResponse~pertSize+isFirstInBlock+pertSize:blockNo+pertSize:pertSign+prevSize';
mm0=fitglm(X,frml,'Distribution','binomial','Link','logit')
%Automated step-down to drop non-sig terms. By default uses a deviance criterion equivalent to LRT test under Wilk's approximation
mm0=mm0.step('Upper',frml,'Criterion','Deviance','PEnter',0,'PRemove',0.05,'Nsteps',Inf)
mm0.plotPartialDependence('pertSize')

%PEnter does not matter, because this is step-down strictly
%This drops the three-way pertSize:blockNo:ID interaction
%individual simple model:
Nsubs=unique(trialData.ID);
if mm0.Formula.HasIntercept
    hi='';
else
    hi='-1';
end
for i=1:length(Nsubs)
   mm{i}=fitglm(X(X.ID==Nsubs(i),:),'leftResponse~pertSize','Distribution','binomial');
   mm{i}=fitglm(X(X.ID==Nsubs(i),:),['leftResponse~' mm0.Formula.LinearPredictor hi],'Distribution','binomial');
   %insignificantBias=mm{i}.Coefficients.pValue(1)>.05;
   %if insignificantBias %Dropping intercept
   %    mm{i}=fitglm(X(X.ID==Nsubs(i),:),'leftResponse~pertSize-1','Distribution','binomial');
   %end
   mm{i}.plotPartialDependence('pertSize')
end

ll=findobj(gca,'Type','Line');
set(ll(1:end-1),'Color',.7*ones(1,3));
set(ll(end),'Color','k','LineWidth',2);
uistack(ll(1:end-1),'bottom')
uistack(ss,'top')
title(['Choice vs. probe size'])
ylabel('% ''<'' (left is slow) responses') 
xlabel('vL>vR         probe size (mm/s)          vL<vR')
set(gca,'XLim',[-350 350])

%% Second plot: same data, but folded to get accuracy estimates and thresholds
subplot(1,3,2)
hold on
trialData.correctResponses=double(trialData.correctResponses);
trialData.correctResponses(isnan(trialData.initialResponse))=nan;
B2=findgroups(abs(trialData.pertSize));
S2=splitapply(@(x) nansum(x)/sum(~isnan(x)),trialData.correctResponses,B2); %Not counting NR responses
S2(S2==0)=NaN;
ap=sort(unique(abs(pp)));
ss=scatter(ap,S2,30,zeros(1,3),'filled','MarkerEdgeColor','w');
E2=splitapply(@(x) nanstd(x==1)/sqrt(sum(~isnan(x))),trialData.correctResponses,B2); %Not counting NR responses
grid on
errorbar(ap,S2,E2,'k','LineStyle','none')
ylabel('% correct') 
xlabel(' |probe size| (mm/s)')
axis([0 360 .5 1])

%Run binomial tests and BH on accuracy results:
clear pval h
for i=2:length(S2)
    Ntrials=sum(~isnan(trialData.initialResponse(abs(trialData.pertSize)==ap(i))));
    correctTrials=Ntrials*S2(i);
    pval(i-1)=binocdf(correctTrials,Ntrials,.5,'upper');
end
h=BenjaminiHochberg(pval, .05, true); %Two-stage BKY procedure    
disp('Significance testing on accuracy:')
if all(h)
    [mp,mpi]=max(pval);
    disp(['All test were significant. Largest p=' num2str(mp) ' for probe size=' num2str(ap(mpi+1)) 'mm/s'])
else
    disp('Some non-sig tests!')
    disp(num2str(ap([0 h])))
end

xx=[0:360];
%Add group fit:
b=mm0.Coefficients.Estimate;
if ~mm0.Formula.HasIntercept
    y=1-1./(1+exp(b(1)*xx));
else %With bias, need to fold data
   %y=1-1./(1+exp(b(1)+b(2)*xx)); 
   y=.5*(1-1./(1+exp(b(1)+b(2)*xx))+1./(1+exp(b(1)+b(2)*-xx)));
end
plot(xx,y,'Color',zeros(1,3),'LineWidth',2)
thm=find(y>.75,1,'first');
%Add individual fits:
for i=1:length(Nsubs)
    XX=X(X.ID==Nsubs(i),:);
    b=mm{i}.Coefficients.Estimate;
    if ~mm{i}.Formula.HasIntercept
        y=1-1./(1+exp(b(1)*xx));
    else %Biased subject, need to fold data
        i
        y=.5*(1-1./(1+exp(b(1)+b(2)*xx))+1./(1+exp(b(1)+b(2)*-xx)));
    end
     plot(xx,y,'Color',.7*ones(1,3),'LineWidth',1)
     th(i)=find(y>.75,1,'first');
end
%legend({'Data','Group fit','Individuals fits'},'Location','SouthEast')
title('Accuracy and psychometric fits')
ll=findobj(gca,'Type','Line','LineWidth',1');
uistack(ll,'bottom')
uistack(ss,'top')
disp('-------------Soft threshold stats (psycho fit):------------')
disp(['Group=' num2str(thm) ', mean=' num2str(mean(th)) ', std=' num2str(std(th)) ', range=[' num2str(min(th)) ',' num2str(max(th)) ']']);
%% Same, but with monotonic fits instead of psychom
subplot(1,3,3)
hold on
ss=scatter(sort(unique(abs(pp))),S2,30,zeros(1,3),'filled','MarkerEdgeColor','w');
grid on
ylabel('% correct') 
xlabel(' |probe size| (mm/s)')
axis([0 360 .5 1])
errorbar(ap,S2,E2,'k','LineStyle','none')

xx=[0:360];
%Add group fit:
    S2(S2==0)=.5;
y=monoLS(S2,2,1,0,1,-1); %Monotonic, concave, function fitting
ap=sort(unique(abs(pp)));
y=interp1(ap,y,xx,'linear');
thm=(find(y>.75,1,'first'));
p1=plot(xx,y,'Color',zeros(1,3),'LineWidth',2);
%Add individual fits:
for i=1:length(Nsubs)
    XX=X(X.ID==Nsubs(i),:);
     %Non-param fit:
    B2=findgroups(abs(XX.pertSize));
    S2=splitapply(@(x) sum(x==1)/sum(~isnan(x)),XX.correctResponses,B2); %Not counting NR responses
    S2(S2==0)=.5;
    y=monoLS(S2,2,1,0,1,-1); %Monotonic, concave, function fitting
    y=interp1(ap,y,xx,'linear');
     p2=plot(xx,y,'Color',.7*ones(1,3),'LineWidth',1);
     th2(i)=(find(y>.75,1,'first'));
end
th=th2;
legend([ss p1 p2],{'Data','Group fit','Individuals fits'},'Location','SouthEast','box','off')
title('Accuracy and monotonic fits')
ll=findobj(gca,'Type','Line','LineWidth',1');
uistack(ll,'bottom')
uistack(ss,'top')
disp('-------------Soft threshold stats (mono fit):------------')
disp(['Group=' num2str(thm) ', mean=' num2str(mean(th)) ', std=' num2str(std(th)) ', range=[' num2str(min(th)) ',' num2str(max(th)) ']']);