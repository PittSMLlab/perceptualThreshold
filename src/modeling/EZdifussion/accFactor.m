function a=accFactor(x,b,s,alpha,noise)
d=dif(x,b,s,alpha,noise);
a=1./(1+exp(-d));
end