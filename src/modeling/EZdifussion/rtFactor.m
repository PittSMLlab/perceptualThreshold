function y=rtFactor(x,b,t73,alpha,noise)
%d=dif(x,b,s,alpha,noise);
d=nl(x,b,t73,alpha);
y=(2./(1+exp(-d)) -1)./d;
y(d==0)=.5;
y=.5*y./noise.^2;
end
