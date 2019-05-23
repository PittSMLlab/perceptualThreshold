function [params, predictedY, Likelihood] = fitGenPsycho(x,y,method,fixedBias)
%Fits the bets fit function of the form (generalized psychometric):
%y=1/1+exp(z), where z=(x^alpha-u)/s;
%The strategy is a greedy search over alpha: for each possible value, we
%evaluate the likelihood under the most likely values of u,s.
%Uses the fitPsycho function to do the fitting given alpha.
%Output: params=[u,s,alpha], 


if nargin<3 || isempty(method)
    method='MLE';
end
if nargin<4
    fixedBias=[];
end
%Check: both x and y are column vectors:
x=reshape(x,length(x),1);
y=reshape(y,length(y),1);
%Discard missing obs:
missingObs=isnan(x) | isnan(y);
x=x(~missingObs);
y=y(~missingObs);

%Define options:
options = optimoptions('fminunc','SpecifyObjectiveGradient',false,'TolX',1e-11,'TolFun',1e-11,'MaxFunctionEvaluations',1e4); %The gradient

%Optimize:
alpha = fmincon(@(alpha) -getThirdOutput(@(v) fitPsycho(nonlin(x,v),y,method,fixedBias),alpha),.9,[],[],[],[],0,[],[],options);
[params,predictedY,Likelihood]=fitPsycho(nonlin(x,alpha),y,method);
params=[params alpha];
end

function third=getThirdOutput(fun, varargin)
[out{1:3}] = fun(varargin{:});
third=out{3};
end
function y2=nonlin(y,alpha)
y2=sign(y).*abs(y).^alpha;
end