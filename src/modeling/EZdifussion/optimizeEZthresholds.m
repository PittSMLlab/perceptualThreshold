%% Optimizing the RT/accuracy curve

%% Desired target:
%[1 2 3 6]->[1 .9 .8 .6] => 1-acc=.08*(rt-1) => acc= 1.08-.08*rt
cost=@(y) sum((y(:,2)-(1.08-.08*y(:,1))).^2);
fun=@(x) getAccRT(x(1),x(2),x(3),3);

%% optimize:
X0=[1 .5 1 3];
X = fminunc(@(x) cost(fun(x)),X0);