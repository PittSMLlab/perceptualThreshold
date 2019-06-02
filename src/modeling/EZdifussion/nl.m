function y=nl(x,b,t73,alpha)
u=(x-b)/t73;
y=sign(u).*abs(u).^alpha;
end