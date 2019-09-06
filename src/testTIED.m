%Testing power-scaling of stimuli as derived in the TIED paper

d=10.^[-10:.5:10]; %This represents the ratio of stimulus intensity s_1/s_2
a=5; %Some constant, does not matter much, except that it be large enough
figure;
hold on;
for gamma=[.1,.3,.5,1]
p=1./(1+exp(-a*tanh(d.^(gamma/2))));
RT=(p-.5)./(sinh(d.^(gamma/2)).*(1-((1-d)./(1+d)).^2).^(gamma/2));
plot(p,gamma*RT/(a),'DisplayName',['\gamma=' num2str(gamma)])
end
legend
axis([.5 1 0 1])

%% Do random-walk simulations
dV=[0,10,25,50,75,100,125,150,200,250,300,350];
mV=1050;
v1=mV+dV/2;
v2=mV-dV/2;
gamma=1;
f=5e-7;
ca=f*(1200.^gamma+900.^gamma)/(1200.^gamma-900.^gamma); %normalized so that the drift rate for dV=300 is the same for any gamma
driftRates=ca*(v2.^gamma-v1.^gamma);
cs=f/5;
diffVariances=cs*(v2.^gamma+v1.^gamma);
d=v1./v2;
a=ca/cs;
expectedP=1./(1+exp(-a*tanh((gamma)*log(d))));
expectedRT=-(2*expectedP-1)./driftRates;

%Sim:
barriers=ones(1,3e4); %Sets barriers and determines simulation time
Nsim=1e3;
[endTime,choice]=simulateEZ(barriers,driftRates,Nsim,diffVariances);

acc=sum(choice==-1)./sum(~isnan(choice));
mRT=nanmean(endTime,1);
nonResp=mean(isnan(choice));

%% PLot:
figure; 
subplot(2,2,1)
plot(acc,mRT,'ko')
hold on
plot(expectedP,expectedRT)

subplot(2,2,2)
plot(dV,acc,'ko')
hold on
plot(dV,expectedP)

subplot(2,2,3)
plot(dV,mRT,'ko')
hold on
plot(dV,expectedRT)