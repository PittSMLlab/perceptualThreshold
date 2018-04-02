%% TEst random walk integrate and fire

%%
threshold=1;
simTime=1; %End time, in a.u.
NstepsPerTimeUnit=100;
M=5;
mus=[-1:.02:1]*threshold*M; 
sigmas=[1e-1,3e-1,5e-1]*threshold*M; %std of random walk step, per unit of time
Niter=5e2;
allEndTime=nan(length(mus),length(sigmas),Niter);
allEndSign=nan(length(mus),length(sigmas),Niter);
bias=0;
fixedNoise=0;
for i=1:length(mus)
    driftRate=mus(i);
    for j=1:length(sigmas)
        noiseRate=sigmas(j);
        [walkPaths,endStep,endSign] = simulateRandomWalk(threshold, driftRate, noiseRate, simTime, NstepsPerTimeUnit,Niter,bias,fixedNoise);
        allEndTime(i,j,:)=deal(endStep);
        allEndSign(i,j,:)=deal(endSign);
    end
end

%% Compute some stats:
meanTime=nanmean(allEndTime,3);
stdTime=nanstd(allEndTime,[],3);
proportionRight=nanmean(allEndSign,3)/2 +.5;
stdRight=nanstd(allEndSign,[],3)/2+.5;
proportionNR=mean(isnan(allEndSign),3);
stdNR=std(isnan(allEndSign),[],3);
%% Do some plots
%Get legends:
legSigmas={};
for i=1:length(sigmas)
    legSigmas{i}=['\sigma=' num2str(sigmas(i),2)];
end
figure('Name',['\theta=' num2str(th) ', T_{end}=' num2str(simTime) ', stepsPerSec=' num2str(NstepsPerTimeUnit)])
subplot(2,2,1)
hold on
plot(mus,meanTime,'LineWidth',2)
cc=get(gca,'ColorOrder');
for i=1:length(sigmas)
   patch([mus'; mus(end:-1:1)'],[meanTime(:,i); meanTime(end:-1:1,i)]+[stdTime(:,i); -stdTime(end:-1:1,i)],cc(i,:),'EdgeColor','none','FaceAlpha',.3) 
end
legend(legSigmas)
xlabel('Speed diff (a.u.)')
ylabel('Mean response time')
set(gca,'YScale','log')
axis tight

subplot(2,2,2)
hold on
plot(mus,proportionRight)
for i=1:length(sigmas)
    %p=fitPsychoGauss(mus,proportionRight(:,i),'MLE');
    %plot(mus,psychoGauss(p,mus),'k')
    %text(-M*.8, .8 -.1*i,['\mu=' num2str(p(1),2)],'Color',cc(i,:))
    %text(-M*.8, .75-.1*i,['\sigma_A=' num2str(p(2),2) '=' num2str(p(2)/sigmas(i)^2,2) '\sigma^2'],'Color',cc(i,:))
    proportionRight(isnan(proportionRight(:,i)),i)=.5;
    p=fitPsycho(mus,proportionRight(:,i),'MLE');
    plot(mus,psycho(p,mus),'k')
    text(M*.1, .6 -.1*i,['a=' num2str(p(1),2)],'Color',cc(i,:))
    text(M*.1, .55-.1*i,['b=' num2str(p(2),2) '=' num2str(p(2)/sigmas(i)^2,2) '\sigma^2'],'Color',cc(i,:))

end
text(M*-.8,.9,['p=(1+e^{(x-a)/b})^{-1}'])
ylabel('% of rightward choices')
xlabel('Speed diff (a.u.)')

subplot(2,2,3)
hold on
plot(mus,proportionNR)
ylabel('% of NR')
xlabel('Speed diff (a.u.)')

subplot(2,2,4)
[h,c]=hist(permute(bsxfun(@minus,allEndTime,nanmean(allEndTime,3)),[3,1,2]),30);
hold on; imagesc(h)
title('Distribution of reaction times around mean')
axis tight