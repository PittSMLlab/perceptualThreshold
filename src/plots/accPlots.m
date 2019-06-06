function [fh,f2]=accPlots(trialData)

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

%% Remove null and no response trials (for accurate counting of DOF in stats)
trialData=trialData(~trialData.noResponse,:);
%% Get probe sizes
B=findgroups(trialData.pertSize); %pertSize>0 means vR>vL
pp=unique(trialData.pertSize);
%% First figure: global stats
fh=figure('Units','pixels','InnerPosition',[100 100 3*300 1*300]);
%set(fh,'PaperUnits','inches','PaperPosition',[0 0 6 2]); %To print at 300dpi
sSize=40;
[cmap,unsignedMap]=probeColorMap(23);
%% First plot: choices as function of probe size
subplot(1,3,1)
hold on
set(gca,'Colormap',cmap);
S=splitapply(@(x) sum(x==-1)/sum(~isnan(x)),trialData.initialResponse,B); %Not counting NR responses
E=splitapply(@(x) nanstd(x==-1)/sqrt(sum(~isnan(x))),trialData.initialResponse,B); %Not counting NR responses
ss=scatter(pp,S,sSize,pp,'filled','MarkerEdgeColor','w');
grid on
ylabel('proportion of left choices') 
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
ylabel('proportion of left choices') 
xlabel('vL>vR         probe (mm/s)         vL<vR')
set(gca,'XLim',[-350 350])

%% Second plot: same data, but folded to get accuracy estimates and thresholds
subplot(1,3,2)
set(gca,'Colormap',unsignedMap);
hold on
trialData.correctResponses=double(trialData.correctResponses);
trialData.correctResponses(isnan(trialData.initialResponse))=nan;
B2=findgroups(abs(trialData.pertSize));
S2=splitapply(@(x) nansum(x)/sum(~isnan(x)),trialData.correctResponses,B2); %Not counting NR responses
S2(S2==0)=NaN;
ap=sort(unique(abs(pp)));
ss=scatter(ap,S2,sSize,ap,'filled','MarkerEdgeColor','w');
E2=splitapply(@(x) nanstd(x==1)/sqrt(sum(~isnan(x))),trialData.correctResponses,B2); %Not counting NR responses
grid on
errorbar(ap,S2,E2,'k','LineStyle','none')
ylabel('accuracy') 
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
        %i
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
set(gca,'Colormap',unsignedMap);
ss=scatter(sort(unique(abs(pp))),S2,sSize,ap,'filled','MarkerEdgeColor','w');
grid on
ylabel('accuracy') 
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
%% Extend panels 10% in width
extendedPanelWidth(fh,.1)

%%
k=1.0986;
X=trialData; %Excludes no response trials already
X.acc=X.correctResponses+.5*trialData.nullTrials;
biases=cellfun(@(x) x.Coefficients.Estimate(1),mm);
slopes=cellfun(@(x) x.Coefficients.Estimate(2),mm);
CI=cell2mat(cellfun(@(x) x.coefCI,mm,'UniformOutput',false));
biasCI=reshape(CI(1,:),2,9);
slopeCI=reshape(CI(2,:),2,9);
f2=figure('Units','pixels','InnerPosition',[100 100 3*300 1*300]);
subplot(1,3,2) %Bias
hold on
%alt: plot % left choices corresponding to bias
%biases=1-1./(1+exp(biases))-.5;
%biasCI=1-1./(1+exp(biasCI))-.5;
%gb=1-1./(1+exp(mm0.Coefficients.Estimate(1)))-.5;
%gCI=1-1./(1+exp(mm0.coefCI))-.5;
%ylabel('excess left choice at \Delta V=0')
%Alt: PSE
biases=-biases./slopes; %PSE
biasCI=-biasCI./slopes; %First order approx to uncertainty: presume that slopes are fixed
gb=-mm0.Coefficients.Estimate(1)/mm0.Coefficients.Estimate(2);
gCI=-mm0.coefCI./mm0.Coefficients.Estimate(2);
ylabel('PSE [-\beta_0/\beta_1] (mm/s)')

%Do plot:
bb=bar(biases,'FaceColor',.6*ones(1,3),'EdgeColor','none');
bar(11,gb,'FaceColor',.2*ones(1,3),'EdgeColor','none');
errorbar(11,gb,gb-gCI(1,1),gCI(1,2)-gb,'k','LineStyle','none','LineWidth',1);
ee=errorbar(1:9,biases,biases-biasCI(2,:),biasCI(1,:)-biases,'k','LineStyle','none','LineWidth',1);

title('bias')
xlabel('subject')
set(gca,'XTick',[1:9,11],'XTickLabel',{'1','2','3','4','5','6','7','8','9','Group'})

subplot(1,3,3) %Slope
hold on
%Alt: plot 1.1/beta_1
slopes=k./slopes;
slopeCI=k./slopeCI;
bb=bar(slopes,'FaceColor',.6*ones(1,3),'EdgeColor','none');
ee=errorbar(1:9,slopes,slopes-slopeCI(2,:),slopeCI(1,:)-slopes,'k','LineStyle','none','LineWidth',1);
gb=k/mm0.Coefficients.Estimate(2);
bar(11,gb,'FaceColor',.2*ones(1,3),'EdgeColor','none');
gCI=k./mm0.coefCI;
errorbar(11,gb,gb-gCI(2,2),gCI(2,1)-gb,'k','LineStyle','none','LineWidth',1);
title('probe size effect')
ylabel('\Delta V for 75% acc. [\approx 1.1/\beta_1] (mm/s)')
set(gca,'XTick',[1:9,11],'XTickLabel',{'1','2','3','4','5','6','7','8','9','Group'})
xlabel('subject')

subplot(1,3,1) %Avg. accuracy vs. avg. RT
hold on
G=findgroups(X.subID);
acc=splitapply(@nanmean,X.acc,G);
eacc=splitapply(@(x) nanstd(x)/sqrt(numel(x)),X.acc,G);
RT=splitapply(@nanmean,X.reactionTime,G);
eRT=splitapply(@(x) nanstd(x)/sqrt(numel(x)),X.reactionTime,G);
errorbar(RT,acc,eacc,'k','LineStyle','none')
ee=errorbar(RT,acc,eRT,'k','Horizontal','LineStyle','none','DisplayName','ste');
ss=scatter(RT,acc,sSize,.5*ones(1,3),'filled','MarkerEdgeColor','w','DisplayName','Individual subject data')
rta=nanmean(X.reactionTime);
acca=nanmean(X.acc);
scatter(rta,acca,sSize,.2*ones(1,3),'filled','MarkerEdgeColor','w')
eRTa=nanstd(X.reactionTime)/sqrt(sum(~isnan(X.reactionTime)));
eacca=nanstd(X.acc)/sqrt(sum(~isnan(X.reactionTime)));
errorbar(rta,acca,eRTa,'k','Horizontal','LineStyle','none','DisplayName','ste');
errorbar(rta,acca,eacca,'k','LineStyle','none')
text(RT+.05,acc-.02,num2str([1:9]'),'Fontsize',7,'FontName','OpenSans')
text(rta-.2,acca+.025,'G','Fontsize',7,'FontName','OpenSans')
ylabel('accuracy')
xlabel('mean RT (s)')
title('indiv. accuracy vs. RT')

extendedPanelWidth(f2,.1)