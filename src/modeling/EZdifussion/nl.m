function y=nl(x,b,s,alpha)
u=(x-b)/s;
y=sign(u).*abs(u).^alpha;
end