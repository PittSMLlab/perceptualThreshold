%This script tests a principled way to scale the DDM variance with stimulus
%size, to see if our observations are consistent with the scaling

stimDiff=[1:5:500]; %Differences in stimuli
meanStim=1050; %Same as in experiment
s1=meanStim+stimDiff/2;
s2=meanStim-stimDiff/2;
stimScale=sqrt(s1.*s2);
stimFrac=sqrt(s1./s2); %Accuracy should depend only on this

drifts=@(gamma,a) a*(s1.^gamma - s2.^gamma);
noises=@(gamma,b) b*(s1.^gamma+s2.^gamma); %Variance
drifts=@(gamma,a) (s1-s2)/75;
noises=@(gamma,b) 1 + .1*drifts(gamma,b).^1;

accuracy=@(gamma,a,b) 1./(1+exp(-drifts(gamma,a)./noises(gamma,b)));
rt=@(gamma,a,b) (2*accuracy(gamma,a,b)-1)./drifts(gamma,a);%mean rt

figure; 
for k=[-10,10] %=1 is equivalent to a constant variance case
    a=((1050+250).^k + (1050-250).^k)/(((1050+250).^k - (1050-250).^k));
    b=.16; %Guarantees good snr for all models at 500 mm/s diff
    subplot(2,2,1); hold on;
    r=rt(k,a,b);
plot(accuracy(k,a,b),r/r(1),'LineWidth',2,'DisplayName',['\gamma=' num2str(k)])
subplot(2,2,2); hold on;
plot(stimDiff,r/r(1),'LineWidth',2,'DisplayName',['\gamma=' num2str(k)])
subplot(2,2,3); hold on;
plot(stimDiff,accuracy(k,a,b),'LineWidth',2,'DisplayName',['\gamma=' num2str(k)])
subplot(2,2,4); hold on;
plot(stimFrac,drifts(k,a)./noises(k,b),'LineWidth',2,'DisplayName',['\gamma=' num2str(k)])
plot(stimFrac,10*tanh(2*log(stimFrac)),'LineWidth',2,'DisplayName',['\gamma=' num2str(k)])
end
legend

%%
figure; hold on;
drifts=@(gamma,a) a*(s1.^gamma - s2.^gamma);
noises=@(gamma,b) b*(s1.^gamma+s2.^gamma); %Variance
drifts1=@(gamma,a) (s1-s2)/75;
noises1=@(gamma,b) 1 + .1*drifts1(gamma,b).^1;

g=5;
a=1/500;
b=.4*1e-5;
plot(stimDiff,noises1(g,b))
n=noises(g,b);
plot(stimDiff,n/n(1))