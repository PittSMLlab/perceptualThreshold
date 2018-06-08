%% Idea: map the RT-accuracy of using different thresholds as a function of drift rate
%Further, if we assume the subjects' cost is some linear combination of
%error rate and reaction time, find the optimal threshold for each rate:
%does any of this resemble true subject behavior?
%Last: is there a way to define a time-varying threshold on a standard DD
%model that would attain these results?

%% Define diffusion model and simulate:
N=3e3;
t=[1:N]/N;
th=1;
a=10;
alpha=1;
beta=2;
thresholdCurve=th*(1+a*t.^alpha).*(1-t.^beta);
thresholdCurve=.1./(t+.15);
%thresholdCurve=1.5*max(1.7,.4./(2*t+.2)).*(1-t.^beta);
%thresholdCurve=max(2*(1-t),1);
thresholdCurve=1.5*(ones(size(t))-t.^3);
beta=100;
thresholdCurve=(min(1.5*(N*t*.01).^.5,100)).*max((1-(1*t).^beta),0);
%thresholdCurve=2.1.*max((1-(1*t).^beta),0);
Nsim=3e3;
drifts=[0 .1 .2 .5 1 2];
P=length(drifts);
[endTime,correctResponse]=simulateEZ(thresholdCurve,drifts,Nsim);
correctRate=reshape(mean(correctResponse),P,1);
correctedRate=reshape(sum(correctResponse)./sum(~isnan(endTime)),P,1); %Excluding non-responses

%% Plot model results
clockStep=.01;
figure;
subplot(4,1,2)
plot([1:N]*clockStep,thresholdCurve,'k')
hold on
plot([1:N]*clockStep,-thresholdCurve,'k')
for l=1:P
    subplot(4,1,2)
    hold on
    ppp=plot(endTime(~isnan(endTime(:,l,1)),l,1),-thresholdCurve(round(endTime(~isnan(endTime(:,l,1)),l,1)/clockStep))'.*(-1).^(correctResponse(~isnan(endTime(:,l,1)),l,1)==1),'o');
    plot(mean(endTime(correctResponse(:,l,1),l,1)),.2,'o','LineWidth',2,'Color',ppp.Color)
    plot(mean(endTime(~correctResponse(:,l,1),l,1)),-.2,'x','LineWidth',2,'Color',ppp.Color)
    subplot(4,1,1)
    hold on
    histogram(endTime(correctResponse(:,l,1),l,1),[0:.33:30],'FaceColor',ppp.Color,'FaceAlpha',.4,'EdgeColor','none','Normalization','probability')
    axis([.5 30 0 .35])
    set(gca,'XScale','log')
    subplot(4,1,3)
    hold on
    histogram(endTime(~correctResponse(:,l,1),l,1),[0:.33:30],'FaceColor',ppp.Color,'FaceAlpha',.5,'EdgeColor','none','Normalization','probability')
    %axis([0 15 0 .35])
    axis([.5 30 0 .35])
    set(gca,'XScale','log')
    
    subplot(4,4,13)
    hold on
    plot(mean(endTime(:,l,1)),mean(correctResponse(:,l,1)),'o','Color',ppp.Color)
    r=corr(endTime(:,l,1),correctResponse(:,l,1));
    plot(mean(endTime(:,l,1))+std(endTime(:,l,1))*[-1 1],mean(correctResponse(:,l,1))+r*std(correctResponse(:,l,1))*[-1 1],'Color',ppp.Color)
    grid on
    title('Accuracy vs. RT')
    subplot(4,4,14)
    hold on
    plot(drifts(l),mean(correctResponse(:,l,1)),'o','Color',ppp.Color)
    grid on
    title('Accuracy vs. drift')
    subplot(4,4,15)
    hold on
    plot(drifts(l),mean(endTime(correctResponse(:,l,1),l,1)),'o','Color',ppp.Color)
    plot(drifts(l),mean(endTime(~correctResponse(:,l,1),l,1)),'x','Color',ppp.Color)
    grid on
    title('RT vs. drift')
    
        subplot(4,4,16)
    hold on
    plot(drifts(l),std(endTime(correctResponse(:,l,1),l,1)),'o','Color',ppp.Color)
    plot(drifts(l),std(endTime(~correctResponse(:,l,1),l,1)),'x','Color',ppp.Color)
    grid on
    title('sRT vs. drift')
end

subplot(4,1,2)
hold on
X=table(correctResponse(:), endTime(:),'VariableNames',{'correctResponse','endTime'});
mm1=fitglm(X,'correctResponse~endTime','Distribution','binomial','Link','identity');
text(10,0,evalc('mm1.disp'),'FontSize',6,'Clipping','off')
subplot(4,4,13)
hold on
t=mean(endTime(:))+std(endTime(:))*[-1 1];
plot(t,t*mm1.Coefficients.Estimate(2)+mm1.Coefficients.Estimate(1),'k')
plot(mean(endTime(:)),mean(correctResponse(:)),'ko')