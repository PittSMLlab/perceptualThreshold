function [pSize,driftRate,noise,delay,bias,t73,alpha]=fitEZ_mine(dataTable,method)
%Fits an EZ-difussion model, as in Wagenmakers et al. 2007
%Modified from wagenmaker's approach to simultaneously fit the same noise
%and delay across all trials, allowing different stimuli only to have different
%drift rates. All parameters returned normalized to the decision threshold
%(i.e. as if a=1).
%The relation between stimulus size and drift rate is modeled as v=(x-b)^alpha;
%Requires that the table contains a list of trials and the following fields
%for each trial: pertSize (type/size) (taken as categorical, just to fit a
%different drift rate for each), reaction time, response (binary). For non
%response trials, reaction time and response should be NaN.
%For plotting purposes, accuracy is defined as % of 1 responses when the
%sign of pert size was >0
%OUTPUT:
%pSize: unique stimulus sizes
%driftRate: (normalized) drift rate associated with each stimulus size
%noise: (normalized) noise level of the random walk [diffusion component]
%delay: (mean?) non-decision time 
%bias: value of stimulus size that corresponds to 50/50 decisions (PSE)
%t73: change in stimulus size needed to go from PSE to 75/25 [actually,
%to exp(1)/(1+exp(1)) ~ 73%]
%alpha: non-linear exponent relating stimulus size to drift rate

if nargin<2 || isempty(method)
    method='RT';
end
%Params init guesses and bounds:
%Params [Td, noise, bias, scale, alpha]
x0=[.5 1 0 100 1];
lb=[-1 1e-5 -Inf 1e-5 0];
ub=[Inf Inf Inf Inf Inf];
noBias=true;
if noBias
    ub(3)=0;
    lb(3)=0;
    dataTable.response(sign(dataTable.pertSize)<0)=1-dataTable.response(sign(dataTable.pertSize)<0); %Flipping responses for negative probes
    dataTable.response(dataTable.pertSize==0)=.5; %Removing null probes, setting arbitrarily to .5 (responses make no sense in this folded set)
    dataTable.pertSize=abs(dataTable.pertSize); %If no bias, everything is symmetric, and can use the extra power of merging both signs in estimating MRT, accuracy, and VRT
    %This matters for fitting only if using the VRT procedure (for MRT and
    %acc, there is no averaging/pooling of data needed).
    %It improves visualization in all cases.
end

noAlpha=true;
if noAlpha
    ub(5)=1;
    lb(5)=1;
end

%Get relevant data:
[~,idx]=sort(dataTable.pertSize); 
dataTable=dataTable(idx,:);%Sort by pert size
G=findgroups(dataTable.pertSize); 
pSize=splitapply(@nanmean,dataTable.pertSize,G); 
MRT=splitapply(@nanmean,dataTable.reactionTime,G);
Acc=splitapply(@nanmean,dataTable.response,G); %This is called accuracy, but is actually the probability of getting one particular response
VRT=splitapply(@nanvar,dataTable.reactionTime,G);

%ALt: average nothing:
%nr=isnan(dataTable.response);
%pSize=dataTable.pertSize(~nr);
%MRT=dataTable.reactionTime(~nr);
%Acc=dataTable.response(~nr);
%VRT=nan(size(MRT));
switch method
    case 'RT' %Fits all parameters to the best mean RT fit
        %Mean rt should satisfy = Td+.5*z/y * (1-exp(y))/(1+exp(y));
        %y=noise^2\(pSize-b)^alpha


        params=lsqnonlin(@(x) x(1)+rtFactor(pSize,x(3),x(4),x(5),x(2)) -MRT ,x0,lb,ub);
        delay=params(1);
        bias=params(3);
        alpha=params(5);
        noise=params(2);
        scale=params(4);
        t73=scale *(noise^2)^(1/alpha);
    case 'VRT' %Fits all parameters except delay to best VRT fit
        
    case 'acc' %Fits difficulty rates to accuracy, delay and noise to mean RT
        %Fit glm to get MLE psychometric, this returns 
        params=fmincon(@(x) fitglm(),x0,[],[],[],[],lb,ub); %Get bias, t73, alpha
        %Get delay and noise:
        
end
driftRate=nl(pSize,bias,scale,alpha);

% %Diagnostics plots:
% difficulty=dif(pSize,bias,scale,alpha,noise);
% accuracy=accFactor(pSize,bias,scale,alpha,noise);
% mix=0; %Another parameter that could be fit: probability of reporting the wrong decision
% accuracy=(1-mix)*accuracy + mix*(1-accuracy);
% figure
% subplot(3,2,1)
% scatter(pSize,MRT,50,.5*ones(1,3),'filled')
% hold on
% plot(pSize,delay+rtFactor(pSize,bias,scale,alpha,noise),'k');
% grid on
% xlabel('Stimulus size')
% ylabel('Mean RT (s)')
% 
% subplot(3,2,2)
% scatter(pSize,Acc,50,.5*ones(1,3),'filled')
% hold on
% plot(pSize,accuracy,'k');
% grid on
% xlabel('Stimulus size')
% ylabel('Accuracy')
% 
% subplot(3,2,3)
% scatter(MRT,Acc,50,.5*ones(1,3),'filled')
% hold on
% plot(delay+rtFactor(pSize,bias,scale,alpha,noise),accuracy,'k');
% grid on
% xlabel('MRT (s)')
% ylabel('Accuracy')
% 
% subplot(3,2,4)
% scatter(pSize,VRT,50,.5*ones(1,3),'filled')
% hold on
% plot(pSize,vrtFactor(pSize,bias,scale,alpha,noise),'k');
% grid on
% xlabel('Stimulus size')
% ylabel('Var RT (s)')
% 
% subplot(3,2,5)
% bar([0,1],[mean(dataTable.reactionTime(dataTable.response==0)) mean(dataTable.reactionTime(dataTable.response==1))]) 

end