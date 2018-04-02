%% Simulate repetitions to see how likely it would be to observe large deviations from the expected behavior in a binary task

%% Case 1: Accuracy vs. reaction times
%load('../timeVsAccuracyUpToAA05.mat') % This was extracted from actual data in 17 blocks over 5 subjects
load('../timeVsAccuracy736664.5474.mat')
t=allAuxX; 
a=allAuxY; %Actual accuracy
p=allAuxZ; %Actual pert size

%Sorting by time:
[~,inds2]=sort(t);
t=t(inds2);
a=a(inds2);
p=p(inds2); 

%Generate a model for the accuracy given the time (we are taking time as given, but we could simulate that
%too if we had a time generation model)
pp=polyfit(log(t),a,1);

%Use the model to predict:
%If we assume uniform time sampling:
t2=[1:length(t)]'*20/length(t); 
%If we leave sampling untouched:
t2=t;
%Comment:
%If we don't sample uniformly then the whole thing gets weird, 
%as there is a much higher chance of finding outliers in the densely 
%sampled intervals than in the sparsely sampled intervals.
%However, in the actual data we do NOT sample uniformly, as the sampling
%times are random. Therefore there truly is a higher chance of finding
%outlier behavior for smaller times.
projectedP=log(t2)*pp(1)+pp(2);
th=.95;
projectedP(projectedP>th)=th; %Saturating out of range values
projectedP(projectedP<.5)=.5;



%Now, simulations: 
Nsim=1e5;
NN=25; %Binwidth
out=nan(Nsim,2);
out2=nan(Nsim,2);
meanA=zeros(size(a));
secMomA=zeros(size(a));
sp=nan(Nsim,1);
r=nan(Nsim,1);
pv=nan(Nsim,1);
pvs=nan(Nsim,1);
%auxNormFactor=conv([th*ones((NN-1)/2,1);projectedP;.5*ones((NN-1)/2,1)],ones(NN,1)/NN,'valid');
%auxNormFactor=sqrt(auxNormFactor.*(1-auxNormFactor));
auxNormFactor=sqrt(projectedP.*(1-projectedP));
auxNormFactor=conv([th*ones((NN-1)/2,1);auxNormFactor;.5*ones((NN-1)/2,1)],ones(NN,1)/NN,'valid');
for i=1:Nsim
   %For each time in t, generate a binary value according to the model probability:
   aux=rand(size(t2)); %Uniform random variables
   simA=aux<projectedP; %If aux is less than projectedP for that time point, assing 1, 0 otherwise. The probability of getting 1 is actually projectedP
   
   %Compute the running avg. of the simulated results:
   rSimA=conv([th*ones((NN-1)/2,1);simA;.5*ones((NN-1)/2,1)],ones(NN,1)/NN,'valid');
   
   %Find the largest deviation between rSimA and the projectedP
   [~, out(i,2)] = max(abs(rSimA-projectedP));
   out(i,1)=rSimA(out(i,2))-projectedP(out(i,2));
   
   [~, out2(i,2)] = max(abs(rSimA-projectedP)./auxNormFactor);
   out2(i,1)=(rSimA(out2(i,2))-projectedP(out2(i,2)))/auxNormFactor(out(i,2));
   
   
   %Compute the mean of simA
   meanA=meanA+simA/Nsim;
   secMomA=secMomA + simA.^2/Nsim;
   
   %Compute Pearson's r and Spearman's r
   [r(i),pv(i)]=corr(log(t2),simA); %Linear corr
   [sp(i),pvs(i)]=corr(log(t2),simA,'type','Spearman'); %Spearman rank corr
   
   %Alternative: find correlation between data and model, instead of
   %log(t), which may be better, since the model is not truly linear in
   %log(t) because of saturation effects.
   %[r(i),pv(i)]=corr(projectedP,simA); %Linear corr
   %[sp(i),pvs(i)]=corr(projectedP,simA,'type','Spearman'); %Spearman rank corr
end

stdA=sqrt((secMomA-meanA.^2)*Nsim/(Nsim-1)); %Unbiased estimator of stdA
meanR=mean(r);
stdR=std(r);
meanPV=mean(pv);
stdPV=std(pv);
meanSP=mean(sp);
stdSP=std(sp);
meanPVS=mean(pvs);
stdPVS=std(pvs);


%% Do some plots:
figure
subplot(2,1,1)
hold on
%[i]=discretize(t,[0:.5:30]);
%for k=1:max(i)
%bar(.5*k-.25,mean(a(i==k)))
%end
ppp2=plot(t2,projectedP,'k','LineWidth',2); %Model predictions
plot(t2,meanA,'o','MarkerFaceColor','k','MarkerEdgeColor','none')
p1=plot(t2,rSimA);
plot(t2,meanA+stdA,'k--')
plot(t2,meanA-stdA,'k--')
aa=conv([th*ones((NN-1)/2,1);a;.5*ones((NN-1)/2,1)],ones(NN,1)/NN,'valid'); %Smoothed responses
[a1,a2]=max(abs(projectedP-aa));
plot(t(a2),aa(a2),'o','MarkerFaceColor','r')
plot(t,aa,'k') %size NN running avg.
legend('Model regression','Avg. simulation result','One simulation','Mean+-1std from sims','','Observed value','Real data')
text(10,.85,['r=' num2str(meanR,2) '\pm' num2str(stdR,2) ' ,p=' num2str(meanPV,2) '\pm' num2str(stdPV,2)],'Color','k')
text(10,.75,['sp=' num2str(meanSP,2) '\pm' num2str(stdSP,2) ' ,p=' num2str(meanPVS,2) '\pm' num2str(stdPVS,2)],'Color','k')


subplot(2,4,5)
hold on
h=histogram(t2(out(:,2)),'Normalization','pdf');
aux2=find(h.Values>.95,1,'first');
legend('empiric PDF')
xlabel('Perturbation size at which the max dev occurred')
ylabel('Bin count')
title('Location of max dev')

subplot(2,4,6)
hold on
h=histogram(abs(out(:,1)),'Normalization','cdf','EdgeColor','none');
aux2=find(h.Values>.95,1,'first');
plot(h.BinEdges(aux2)*[1 1],[0 1],'k','LineWidth',2)
plot(a1*[1 1],[0 1],'r','LineWidth',2)
legend('empiric CDF','95% threshold','Observed value')
xlabel('Absolute deviation from model')
ylabel('Cumulative probability')
title('CDF of max dev. under null')

subplot(2,4,7)
hold on
h=histogram(t2(out2(:,2)),'Normalization','pdf');
xlabel('Perturbation size')
ylabel('Bin count')
title('Location of max dev (corrected by std)')
text(5,.2,'This should reflect the sampling distribution')
h2=histogram(t2,'Normalization','pdf');
uistack(h2,'bottom')
legend([h2 h],{'sampling PDF','empiric PDF'})

subplot(2,4,8)
hold on
h=histogram(abs(out2(:,1)),'Normalization','cdf','EdgeColor','none');
aux2=find(h.Values>.95,1,'first');
plot(h.BinEdges(aux2)*[1 1],[0 1],'k','LineWidth',2)
plot(a1*[1 1]/sqrt(a1*(1-a1)),[0 1],'r','LineWidth',2)
legend('empiric CDF','95% threshold','observed value')
xlabel('Abs. deviation, corrected by std. from model (z-score)')
title('CDF of max dev. (corrected by std) under null')
ylabel('Cumulative probability')


%% Case 2: accuracy as a function of perturbation size
load('../reactionsUpToAA05.mat')
pp=unique(pertSize);
%Discard NaN reactions:
pertSize=pertSize(~isnan(reactionTime));
reactionSign=reactionSign(~isnan(reactionTime));
reactionTime=reactionTime(~isnan(reactionTime));

%Generate a model:
[modelPsycho,~] = fitPsycho(pertSize(pertSize~=250),reactionSign(pertSize~=250)==-1,'MLEsat');
projectedP=psycho(modelPsycho,pertSize);
summaryP=nan(size(pp));
empiricP=nan(size(pp));
for j=1:length(pp) %Group results by perturbation size
    summaryP(j)=mean(projectedP(pertSize==pp(j)));
    empiricP(j)=mean(reactionSign(pertSize==pp(j))==1);
end

%Now, simulate:
Nsim=1e6;
out=nan(Nsim,2);
out2=nan(Nsim,2);
meanA=zeros(size(summaryP));
secMomA=zeros(size(summaryP));
sp=nan(Nsim,1);
r=nan(Nsim,1);
pv=nan(Nsim,1);
pvs=nan(Nsim,1);
for i=1:Nsim
   %For each repetition, generate a binary value according to the model probability:
   aux=rand(size(projectedP)); %Uniform random variables
   simA=aux<projectedP; %If aux is less than projectedP for that time point, assing 1, 0 otherwise. The probability of getting 1 is actually projectedP
   
   summaryA=nan(size(pp));
   for j=1:length(pp) %Group results by perturbation size
       summaryA(j)=mean(simA(pertSize==pp(j)));
   end
   %Find the largest deviation between summaryA and the summaryP
   [~, out(i,2)] = max(abs(summaryA-summaryP));
   out(i,1)=summaryA(out(i,2))-summaryP(out(i,2));
   
   [~, out2(i,2)] = max(abs(summaryA-summaryP)./sqrt(summaryP.*(1-summaryP)));
   out2(i,1)=(summaryA(out2(i,2))-summaryP(out2(i,2)))/sqrt(summaryP(out2(i,2)).*(1-summaryP(out2(i,2))));
   
   %Compute the mean of simA
   meanA=meanA+summaryA/Nsim;
   secMomA=secMomA + summaryA.^2/Nsim;
   
      %Compute Pearson's r and Spearman's r
   [r(i),pv(i)]=corr(projectedP,simA); %Linear corr
   [sp(i),pvs(i)]=corr(projectedP,simA,'type','Spearman'); %Spearman rank corr
   
end

stdA=sqrt((secMomA-meanA.^2)*Nsim/(Nsim-1)); %Unbiased estimator of stdA

meanR=mean(r);
stdR=std(r);
meanPV=mean(pv);
stdPV=std(pv);
meanSP=mean(sp);
stdSP=std(sp);
meanPVS=mean(pvs);
stdPVS=std(pvs);

%% Do some plots:
figure
subplot(2,1,1)
hold on
ppp2=plot(pp,summaryP,'k','LineWidth',2); %Model predictions
plot(pp,meanA,'o','MarkerFaceColor','k','MarkerEdgeColor','none')
p1=plot(pp,summaryA);
plot(pp,meanA+stdA,'k--')
plot(pp,meanA-stdA,'k--')
plot(250,.58,'o','MarkerFaceColor','r')
plot(pp,1-empiricP,'k')
legend('Model regression','Avg. simulation result','One simulation','Mean+-1std from sims','','Observed value')
text(150,.25,['r=' num2str(meanR,2) '\pm' num2str(stdR,2) ' ,p=' num2str(meanPV,2) '\pm' num2str(stdPV,2)],'Color','k')
text(150,.15,['sp=' num2str(meanSP,2) '\pm' num2str(stdSP,2) ' ,p=' num2str(meanPVS,2) '\pm' num2str(stdPVS,2)],'Color','k')

subplot(2,4,5)
hold on
h=histogram(pp(out(:,2)),'Normalization','count');
aux2=find(h.Values>.95,1,'first');
legend('empiric PDF')
xlabel('Perturbation size at which the max dev occurred')
ylabel('Bin count')
title('Location of max dev')

subplot(2,4,6)
hold on
h=histogram(abs(out(:,1)),'Normalization','cdf','EdgeColor','none');
aux2=find(h.Values>.95,1,'first');
plot(h.BinEdges(aux2)*[1 1],[0 1],'k','LineWidth',2)
plot(.32*[1 1],[0 1],'r','LineWidth',2)
legend('empiric CDF','95% threshold','Observed value')
xlabel('Absolute deviation from model')
ylabel('Cumulative probability')
title('CDF of max dev. under null')

subplot(2,4,7)
hold on
h=histogram(pp(out2(:,2)),'Normalization','count');
aux2=find(h.Values>.95,1,'first');
legend('empiric PDF')
xlabel('Perturbation size')
ylabel('Bin count')
title('Location of max dev (corrected by std)')

subplot(2,4,8)
hold on
h=histogram(abs(out2(:,1)),'Normalization','cdf','EdgeColor','none');
aux2=find(h.Values>.95,1,'first');
plot(h.BinEdges(aux2)*[1 1],[0 1],'k','LineWidth',2)
plot(.32*[1 1]/sqrt(summaryP(pp==250)*(1-summaryP(pp==250))),[0 1],'r','LineWidth',2)
legend('empiric CDF','95% threshold','observed value')
xlabel('Abs. deviation, corrected by std. from model (z-score)')
title('CDF of max dev. (corrected by std) under null')
ylabel('Cumulative probability')