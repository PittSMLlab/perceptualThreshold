function a=accFactor(x,b,t73,alpha,noise)
d=nl(x,b,t73,alpha);
a=1./(1+exp(-d));
end