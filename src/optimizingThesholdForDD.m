%% Idea: map the RT-accuracy of using different thresholds as a function of drift rate
%Further, if we assume the subjects' cost is some linear combination of
%error rate and reaction time, find the optimal threshold for each rate:
%does any of this resemble true subject behavior?
%Last: is there a way to define a time-varying threshold on a standard DD
%model that would attain these results?

%% Define diffusion model and simulate:
N=3e4; %Simulation steps
th=a/2; %Threshold, arbitrary scale
M=2e3; %Number of simulations for each parameter pair
drifts=[0,.05,.1,.2,.4,.8,1];
noises=1;
Q=length(noises); 
P=length(drifts);
%correctRate=zeros(P,Q);
correctResponse=false(M,P,Q);
endTime=nan(M,P,Q);
clockStep=.001;
thresholdCurve=th*ones(1,N);
%thresholdCurve=1.1*th*(1+(2*[1:N]/N).^.5);
for l=1:Q
    s=sqrt(noises(l));
    for j=1:P
        m=drifts(j);
        for k=1:M %Number of sims
            x=zeros(N,1);
            for i=2:N
                x(i)=x(i-1)+m*clockStep+s*randn*sqrt(clockStep);
                if abs(x(i))>=thresholdCurve(i)
                    endTime(k,j,l)=(i)*clockStep;
                    correctResponse(k,j,l)=x(i)>0;
                    %correctRate(j,l)=correctRate(j,l)+(x(i)>0)/M;
                    break
                end
            end
        end
    end
end
correctRate=reshape(mean(correctResponse),P,Q);
correctedRate=reshape(sum(correctResponse)./sum(~isnan(endTime)),P,Q); %Excluding non-responses

%% Plot model results
figure;
subplot(4,1,2)
plot([1:N]*clockStep,thresholdCurve,'k')
hold on
plot([1:N]*clockStep,-thresholdCurve,'k')
for l=1:2:P
    subplot(4,1,2)
    hold on
    ppp=plot(endTime(~isnan(endTime(:,l,1)),l,1),-thresholdCurve(round(endTime(~isnan(endTime(:,l,1)),l,1)/clockStep))'.*(-1).^(correctResponse(~isnan(endTime(:,l,1)),l,1)==1),'o');
    subplot(4,1,1)
    hold on
    histogram(endTime(correctResponse(:,l,1),l,1),[0:.33:30],'FaceColor',ppp.Color,'FaceAlpha',.4,'EdgeColor','none','Normalization','probability')
    axis([0 30 0 .15])
    subplot(4,1,3)
    hold on
    histogram(endTime(~correctResponse(:,l,1),l,1),[0:.33:30],'FaceColor',ppp.Color,'FaceAlpha',.5,'EdgeColor','none','Normalization','probability')
    axis([0 30 0 .15])
end