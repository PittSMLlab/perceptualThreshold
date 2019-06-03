function [signed,unsigned]=probeColorMap(N)
if nargin<1
    N=23; %Different probe sizes
end
gamma=1;
ex2=[0.2314    0.2980    0.7529];
ex1=[0.7255    0.0863    0.1608];
x=[0:N-1]/(N-1);
signed=[bsxfun(@plus,ex1.^(1/gamma),bsxfun(@times,.8-ex1.^(1/gamma),x'));bsxfun(@plus,ex2.^(1/gamma),bsxfun(@times,.8-ex2.^(1/gamma),fliplr(x)'))].^gamma;
ex0=[.2 .2 .2];
unsigned=flipud([bsxfun(@plus,ex0.^(1/gamma),bsxfun(@times,.8-ex0.^(1/gamma),x'))]);
end

