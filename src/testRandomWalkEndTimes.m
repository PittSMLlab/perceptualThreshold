%% Simulate random walk threshold times
N=1e4;
s=.08;
th=1;
M=1e4;
P=10;
correct=false(M,P);
t=nan(M,P);
for j=1:P
    m=(j-1)/400;
    for k=1:M %Number of sims
        x=zeros(N,1);
        for i=2:N
            x(i)=x(i-1)+m+s*randn;
            if abs(x(i))>th
                break
            end
        end
        t(k,j)=i;
        correct(k,j)=x(i)>th;
    end
end

%%
figure
subplot(2,1,1)
histogram(t,'Normalization','pdf')
hold on
c=(th/s)^2;
mu=0;
x=[.01:.01:N];
plot(x,1.8*sqrt(c/(2*pi))*exp(-c./(2*(x-mu)))./(x-mu).^(3/2))

subplot(2,1,2)
plot([0:P-1]/400,mean(correct))
