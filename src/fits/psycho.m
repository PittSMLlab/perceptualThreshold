function [f,g]=psycho(u,x)
x=reshape(x,length(x),1);
if length(u)<3 %No saturation parameter
    u(3)=0;
elseif u(3)>1 || u(3)<0
    error('')
end
if length(u) < 4 %No alpha value
    u(4)=1; 
end

    f =u(3) + (1-2*u(3))./(1+exp(-nonlin((x-u(1))/u(2),u(4))));
    
    g = bsxfun(@times, f.*(1-f) , [-ones(size(x))./u(2)  -(x-u(1))/u(2)^2]);
    %To do: g with alpha ~=1
end

function x2=nonlin(x,alpha)
x2=sign(x).*abs(x).^alpha;
end