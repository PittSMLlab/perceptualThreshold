function [p,peval,L] = fitPsycho(x,y,mode)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Model:
%p(y_i) = 1- 1/(1+exp((x_i-mu)/sigma))

%Check: both x and y are column vectors:
x=reshape(x,length(x),1);
y=reshape(y,length(y),1);

options = optimoptions('fminunc','SpecifyObjectiveGradient',false,'TolX',1e-11,'TolFun',1e-11); %The gradient
u0=[nanmean(x),nanstd(x)];
switch mode
    case 'MLE'
        %This mode assumes the variables 'y' are BINARY and fits the model
        %which results in the highest likelihood of the observed data
        p = fminunc(@(u) MLE(u,x,y),u0,options);
        ff=psycho(p,x);
        L=prod(y.*(ff) + (1-y).*(1-ff));
    case 'MSE'
        %This mode accepts continuous variables distributed in the 0-1 range 
        %(not just binary) and assumes they are measurements of 'p'
        %directly, fitting the best curve
        
        p = fminunc(@(u) MSE(u,x,y),u0,options);
    case 'MAE'
        %L1 norm distance minimization
        p = fminunc(@(u) MAE(u,x,y),u0,options);
        %L1 norm minimzation just returns a
%very tight psychom function, as it is minimized by choosing a function
%that is exactly 1 at speeds where p>.5 and 0 at speeds where p<.5, not
%quite sure what would happen if the empirical data did not show a
%monotonic increase with p
    case 'MLEsat'
        u0=[u0 .9];
        options = optimoptions('fmincon','SpecifyObjectiveGradient',false,'TolX',1e-11,'TolFun',1e-11); %The gradient
        p = fmincon(@(u) MLE(u,x,y),u0,[],[],[],[],[-Inf -Inf 0],[Inf Inf 1],[],options);
        ff=psycho(p,x);
        L=prod(y.*(ff) + (1-y).*(1-ff));
end
peval=psycho(p,x);

function [f,g] = MSE(u,x,y) %Defining loss function for MSE estimation
    [ff,gg]=psycho(u,x);
    f=sum((y-ff).^2);
    g = sum(bsxfun(@times,2*(y-1./(1+exp((x-u(1))/u(2)))),gg),1);
end

function [f,g] = MAE(u,x,y) %Defining loss function for MAE estimation
    [ff,gg]=psycho(u,x);
    f=sum(abs(y-ff));
    g = sum(bsxfun(@times,-sign(y-1./(1+exp((x-u(1))/u(2)))),gg),1);
end

function [f,g] = MLE(u,x,y) %Defining loss function for MLE
    [ff,gg]=psycho(u,x);
    %Linear comb:
    %f=-sum(log(y.*(ff) + (1-y).*(1-ff))); %log-likelihood actually
    %g=-1*sum(bsxfun(@times,(2*y -1)./(y.*(ff) + (1-y).*(1-ff)),gg),1);
    %Product comb: (proper way)
    f=-sum(log(ff.^y .* (1-ff).^(1-y))); 
    g=[sum(y-ff)/u(2) sum((y-ff).*(x-u(1)))/u(2)^2]; %Doxy
    
end

% function [f,g] =MLEsat(u,x,y)
%     [ff,gg]=psycho(u,x);
%     th=u(3);
%     ff(ff>(1-th))=1-th;
%     ff(ff<th)=th;
%     f=-sum(log(ff.^y .* (1-ff).^(1-y)));
%     g=zeros(length(f),2); %Doxy
% end

end

