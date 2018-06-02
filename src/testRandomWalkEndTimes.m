%% load empirical data for comparison
%loadAllDataIntoTable
%% Simulate random walk threshold times
N=30; %Simulation steps
s=.62;
th=1;
M=1e5;
drifts=[0 10 25 50 75 100 125 150 200 250 300 350];
P=length(drifts);
correct=false(M,P);
t=nan(M,P);
for j=1:P
    m=drifts(j)/5e2+drifts(j)^2/5e4;
    for k=1:M %Number of sims
        x=zeros(N,1);
        for i=2:N
            x(i)=x(i-1)+m+s*randn;
            if abs(x(i))>th
                t(k,j)=i;
                correct(k,j)=x(i)>th;
                break
            end
        end
    end
end



%%
figure
subplot(2,2,1)
for j=1:P
histogram(t(:,j),'Normalization','probability','EdgeColor','none','FaceAlpha',.2)
hold on
end
title('Simulated')
axis([0 15 0 .5])
subplot(2,2,2)
for j=1:P
        dataT=superSuperT(superSuperT.pertSize==drifts(j) | superSuperT.pertSize==-drifts(j),:);
        histogram(dataT.reactionStride,'Normalization','probability','EdgeColor','none','FaceAlpha',.2)
        hold on
        empAcc(j,:)=sign(dataT.initialResponse)==-sign(dataT.pertSize);
        empMeanRS(j)=nanmean(dataT.reactionStride);
        empStdRS(j)=nanstd(dataT.reactionStride);
end
title('Empirical')
axis([0 15 0 .5])
%c=(th/s)^2;
%mu=0;
%x=[.01:.01:N];
%plot(x,1.8*sqrt(c/(2*pi))*exp(-c./(2*(x-mu)))./(x-mu).^(3/2))

subplot(2,2,4)
title('End time')
plot(drifts,nanmean(t,1),'LineWidth',2,'DisplayName','Mean')
hold on
plot(drifts,empMeanRS,'o','MarkerSize',5,'LineWidth',2,'DisplayName','Empirical mean')
%plot([0:P-1]/600,1e-3./([0:P-1]/600).^2,'LineWidth',1,'DisplayName','Mean (theory)')
plot(drifts,nanstd(t,1),'LineWidth',2,'DisplayName','Std')
plot(drifts,empStdRS,'o','MarkerSize',5,'LineWidth',2,'DisplayName','Empirical std')
xlabel('Perturbation (mm/s)')
ylabel('Reaction stride (s)')
legend

subplot(2,2,3)
plot(drifts,mean(correct),'LineWidth',2,'DisplayName','Model')
hold on
plot(drifts,mean(empAcc,2),'o','MarkerSize',5,'LineWidth',2,'DisplayName','Model')
grid on
title('Accuracy')
ylabel('Correct decisions')