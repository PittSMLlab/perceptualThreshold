%Estimating how many binary choices we need to measure to get an accurate
%estimate of PSE shift

%Assumptions: 
%initial (before adaptation) PSE=0, exactly known
%p(left response) = 1/(1+exp((PSE-DeltaV)/sigma)
%sigma=100mm/s, value estimated from pilot data, assumed to be known exactly
%Subject 1: true PSE shift is 150mm/s
%Subject 2: true PSE shift is 300mm/s
%All measurements are conducted for DeltaV=200mm/s

%The design target is as follows: that two subjects independently estimated
%using N binary choices, but whose true PSE difference is 150mm/s, are not
%estimated to have the opposite ordering more than 5% of the time.
%Equivalently, that the difference in PSE estimates of these subjects is
%positive 95% of the time.
%Under some general conditions (probability distribution falls at least as
%fast as a gaussian and symmetric?) a sufficient condition is that a single subject's
%97.5% CI is smaller than +-150mm/s

%% simulate estimation of PSE in two subjects with different true PSEs:
truePSE=[150,300];
sigma=100;
measurePoint=200; %Arbitrary and asymmetric
p1=1./(1+exp((truePSE(1)-measurePoint)/sigma));
p2=1./(1+exp((truePSE(2)-measurePoint)/sigma));
N=15;
Nreps=5000;
clear PSEestim
for i=1:Nreps
   x(1)=binornd(N,p1,1,1);
   x(2)=binornd(N,p2,1,1);
   for k=1:2
        responses=zeros(N,1);
        responses(1:x(k))=1;
        unos=ones(N,1);
        tbl=table(responses,unos);
        glm1=fitglm(tbl,'responses~unos-1','Distribution','binomial');
        PSEestim(i,k)=-sigma*glm1.Coefficients.Estimate+measurePoint;
        if k==1
            PSEci(i,:)=-sigma*glm1.coefCI + measurePoint;
        end
   end
end
figure; histogram(PSEestim(:,1),-100:20:500); hold on; histogram(PSEestim(:,2),-100:20:500)

PSEdiff=PSEestim(:,2)-PSEestim(:,1); %Most of the time this should result in a positive value
disp('Probability that we mis-order two subjects whose true PSE difference is 150mm/s:')
2*sum(PSEdiff<0)/numel(PSEdiff) %Twice the percent of inverted ordering (2x because we do not know the order a priori)
disp('Subject 1: true PSE=150mm/s, proability that it is estimated at least below 250mm/s:')
sum(PSEestim(:,1)<250)/Nreps
disp('For each subject, probability that its PSE is estimated more than 50mm/s away from its true value:')
sum(abs(PSEestim-truePSE)>50)/Nreps
disp('For each subject, probability that its PSE is estimated more than 100mm/s away from its true value:')
sum(abs(PSEestim-truePSE)>100)/Nreps
disp('For each subject, probability that its PSE is estimated more than 150mm/s away from its true value:')
sum(abs(PSEestim-truePSE)>150)/Nreps
disp('For subject 1, probability that the returned 95% CI includes 300mm/s:')
sum(PSEci(:,1)>300)/Nreps
disp('For subject 1, probability that the returned 95% CI includes 0mm/s:')
sum(PSEci(:,2)<0)/Nreps