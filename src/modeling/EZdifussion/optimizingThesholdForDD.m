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
beta=10;
thresholdCurve=(min(1.7*sqrt(.001+N*t*.01),4)).*max((1-(1.1*t).^beta),0);
%thresholdCurve=.5*ones(size(t));
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
    axis([0 15 0 .35])
    subplot(4,1,3)
    hold on
    histogram(endTime(~correctResponse(:,l,1),l,1),[0:.33:30],'FaceColor',ppp.Color,'FaceAlpha',.5,'EdgeColor','none','Normalization','probability')
    axis([0 15 0 .35])
        subplot(4,4,13)
    hold on
    plot(mean(endTime(:,l,1)),mean(correctResponse(:,l,1)),'o','Color',ppp.Color)
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