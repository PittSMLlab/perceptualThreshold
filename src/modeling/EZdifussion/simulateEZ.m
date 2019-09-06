function [endTime,choice]=simulateEZ(thresholdCurve,drifts,Nsim,noises)
%Simulate random walks until first cross to given thrshold curve
%Drift rates and threshold curves are normalized to noise levels if noises
%(variance)
%not given. 
%If noises given, the parameters are somewhat redundant: thresholds, drifts
%and sqrt(noises) can be arbitrarily scaled by some constant and results
%should remain unchanged.
clockStep=1; %Assuming rates normalized to clock step
N=length(thresholdCurve);
if nargin<3 || isempty(Nsim)
    Nsim=1e4; %Number of simulations for each parameter pair
end
if nargin<2 || isempty(drifts)
    drifts=[0,.1,.2,.4,.75,1];
end
if nargin<4
    noises=ones(size(drifts)); %Presuming normalization to noise
end
P=length(drifts);
choice=nan(Nsim,P);
endTime=nan(Nsim,P);
for j=1:P
    for k=1:Nsim %Number of sims
        x=zeros(N,1);
        m=drifts(j);%+.5*randn;
        s=sqrt(noises(j));
        for i=2:N
            x(i)=x(i-1)+m*clockStep+s*randn*sqrt(clockStep);
            if abs(x(i))>=thresholdCurve(i)
                endTime(k,j)=(i-1)*clockStep;
                choice(k,j)=sign(x(i));
                break
            end
        end
    end
end
