function [endTime,correctResponse]=simulateEZ(thresholdCurve,drifts,Nsim)
N=length(thresholdCurve);
if nargin<3 || isempty(Nsim)
    Nsim=1e4; %Number of simulations for each parameter pair
end
if nargin<2 || isempty(drifts)
    drifts=[0,.1,.2,.4,.75,1];
end
noises=1;
P=length(drifts);
correctResponse=false(Nsim,P);
endTime=nan(Nsim,P);
clockStep=.01;
s=sqrt(noises);
for j=1:P
    for k=1:Nsim %Number of sims
        x=zeros(N,1);
        m=drifts(j);%+.5*randn;
        for i=2:N
            x(i)=x(i-1)+m*clockStep+s*randn*sqrt(clockStep);
            if abs(x(i))>=thresholdCurve(i)
                endTime(k,j)=i*clockStep;
                correctResponse(k,j)=x(i)>0;
                break
            end
        end
    end
end
