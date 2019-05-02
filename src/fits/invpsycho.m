function [x]=invpsycho(u,y)
%Inverse psychometric function, such that x=invpsycho(u,psycho(u,x)) for
%any x. Notice that because the psychometric function has a bounded image,
%the composition of functions in the opposite order works only for some
%values of the input argument
y=reshape(y,length(y),1);
if length(u)<3
    u(3)=0;
elseif u(3)>1 || u(3)<0
    error('')
end
if any(y>(1-u(3))) || any(y<=0)
    warning('invPsycho is not defined for arguments outside [0 u(3)]')
    y(y<=0)=.01;
    y(y>=(1-u(3)))=.99*(1-u(3));
end
    
    x=u(1) -u(2)*log(((1-y)-u(3))./(y-u(3)));

end
