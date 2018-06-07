function [out]=getAccRT(th,a,alpha,beta)
N=3e3;
t=[1:N]/N;
thresholdCurve=th*(1+a*t.^alpha).*(1-t.^beta);
drifts=[0,.1,.2,.5,.75,1];
Nsim=1e4;
[endTime,correctResponse]=simulateEZ(thresholdCurve,drifts,Nsim);

P=length(drifts);
correctRate=reshape(mean(correctResponse),P,1);
%correctedRate=reshape(sum(correctResponse)./sum(~isnan(endTime)),P); %Excluding non-responses
meanRT=nanmean(endTime)';
stdRT=nanstd(endTime)';
out=[meanRT correctRate];
end