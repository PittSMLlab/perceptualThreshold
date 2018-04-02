function [f,g]=psycho(u,x)
x=reshape(x,length(x),1);
if length(u)<3
    u(3)=0;
elseif u(3)>1 || u(3)<0
    error('')
end
    f =u(3) + (1-2*u(3))./(1+exp(-(x-u(1))/u(2)));
    
    g = bsxfun(@times, f.*(1-f) , [-ones(size(x))./u(2)  -(x-u(1))/u(2)^2]);
end
