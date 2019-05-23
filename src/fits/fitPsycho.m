function [params,predictedY,logL] = fitPsycho(x,y,mode,fixedBias)
%Model:
%p(y_i) = 1- 1/(1+exp((x_i-mu)/sigma))
%Fitted through max likelihood (MLE), min squared error (MSE), or min L1
%distance (MAE). 


if nargin<3 || isempty(mode)
    mode='MLE';
end
if nargin<4
    fixedBias=[];
end
%Check: both x and y are column vectors:
x=reshape(x,length(x),1);
y=reshape(y,length(y),1);

%options = optimoptions('fminunc','SpecifyObjectiveGradient',false,'TolX',1e-11,'TolFun',1e-11,'MaxFunctionEvaluations',1e4); %The gradient
options = optimoptions('fmincon','Display','off','SpecifyObjectiveGradient',false,'TolX',1e-11,'TolFun',1e-11,'MaxFunEvals',1e4); %The gradient
  
u0=[nanmean(x),nanstd(x)];
lb=[-Inf 0];
ub=[Inf Inf];
switch mode
    case 'MLEsat'
        u0=[u0 .9];
        params = fmincon(@(u) MLE(u,x,y),u0,[],[],[],[],[lb 0],[ub .2],[],options);
    otherwise %MLE, MSE, or MAE
        eval(['fun=@(u) ' mode '(u,x,y);']); 
        if isempty(fixedBias)
            params = fmincon(fun,u0,[],[],[],[],lb,ub,[],options);
        else
           f2=@(u) fun([fixedBias u]);
           params = fmincon(f2,u0(2),[],[],[],[],lb,ub,[],options);
           params=[fixedBias params];
        end
end
predictedY=psycho(params,x);
logL=sum(y.*log0(predictedY) + (1-y).*log0(1-predictedY));

end

function [f] = MSE(u,x,y) %Defining loss function for MSE estimation
    [ff,gg]=psycho(u,x);
    f=sum((y-ff).^2);
    g = sum(bsxfun(@times,2*(y-1./(1+exp((x-u(1))/u(2)))),gg),1);
end

function [f] = MAE(u,x,y) %Defining loss function for MAE estimation
    [ff,gg]=psycho(u,x);
    f=sum(abs(y-ff));
    g = sum(bsxfun(@times,-sign(y-1./(1+exp((x-u(1))/u(2)))),gg),1);
end

function [f] = MLE(u,x,y) %Defining loss function for MLE
    [ff,~]=psycho(u,x);
    %Linear comb:
    %f=-sum(log(y.*(ff) + (1-y).*(1-ff))); %log-likelihood actually
    %g=-1*sum(bsxfun(@times,(2*y -1)./(y.*(ff) + (1-y).*(1-ff)),gg),1);
    %Product comb: (proper way)
    %f=-sum(log(ff.^y .* (1-ff).^(1-y))); 
    f=-sum(y.*log0(ff) + (1-y).*log0(1-ff));
    g=[sum(y-ff)/u(2) sum((y-ff).*(x-u(1)))/u(2)^2]; %Doxy
    
end

function l=log0(x)
   if x~=0
       l=log(x);
   else
       l=-Inf;
   end
end