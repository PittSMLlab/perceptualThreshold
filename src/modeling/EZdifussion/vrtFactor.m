function v=vrtFactor(x,b,s,alpha,noise)
d=dif(x,b,s,alpha,noise);
z=1/noise^2;
a=accFactor(x,b,s,alpha,noise);
v=.5*(z.^2./d.^3).*(-2*d.*exp(-d)-exp(-2*d)+1).*a.^2;
end