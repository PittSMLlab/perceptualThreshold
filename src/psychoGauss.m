function [f,g]=psychoGauss(u,x)
x=reshape(x,length(x),1);
if length(u)<3
    u(3)=0;
elseif u(3)>1 || u(3)<0
    error('')
end
    f =u(3) + (1-2*u(3)).* (.5 + erf((x-u(1))/(sqrt(2)*u(2)))/2);
    
    g = zeros(size(f));
end
