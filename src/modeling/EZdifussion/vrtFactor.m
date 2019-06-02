function v=vrtFactor(x,b,t73,alpha,noise)
d=nl(x,b,t73,alpha);
z=1/noise^2;
a=accFactor(x,b,t73,alpha);
v=.5*(z.^2./d.^3).*(-2*d.*exp(-d)-exp(-2*d)+1).*a.^2;
end