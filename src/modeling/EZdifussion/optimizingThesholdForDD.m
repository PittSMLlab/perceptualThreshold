%% Idea: map the RT-accuracy of using different thresholds as a function of drift rate
%Further, if we assume the subjects' cost is some linear combination of
%error rate and reaction time, find the optimal threshold for each rate:
%does any of this resemble true subject behavior?
%Last: is there a way to define a time-varying threshold on a standard DD
%model that would attain these results?

%% Define diffusion model and simulate:
N=3e3;
t=[1:N]/N;
th=.88;
a=.7;
alpha=.09;
beta=10;
thresholdCurve=th*(1+a*t.^alpha).*(1-t.^beta);
Nsim=1e4;
drifts=[0 .1 .2 .5 .75 1];
[endTime,correctResponse]=simulateEZ(thresholdCurve,drifts,Nsim);
correctRate=reshape(mean(correctResponse),P,Q);
correctedRate=reshape(sum(correctResponse)./sum(~isnan(endTime)),P,Q); %Excluding non-responses

%% Plot model results
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
    axis([0 30 0 .25])
    subplot(4,1,3)
    hold on
    histogram(endTime(~correctResponse(:,l,1),l,1),[0:.33:30],'FaceColor',ppp.Color,'FaceAlpha',.5,'EdgeColor','none','Normalization','probability')
    axis([0 30 0 .25])
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
end