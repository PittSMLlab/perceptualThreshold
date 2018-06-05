%% Simulate random walk threshold times
N=30; %Simulation steps
th=1; %Threshold, arbitrary scale
M=2e4; %Number of simulations for each parameter pair
drifts=[0:.1:1];
noises=[0:.1:2];
Q=length(noises);
P=length(drifts);
correctRate=zeros(P,Q);
endTime=nan(M,P,Q);
for l=1:Q
    s=noises(l);
    for j=1:P
        m=drifts(j);
        for k=1:M %Number of sims
            x=zeros(N,1);
            for i=2:N
                x(i)=x(i-1)+m+s*randn;
                if abs(x(i))>th%/(1-5*i/N))
                    endTime(k,j,l)=i-1;
                    correctRate(j,l)=correctRate(j,l)+(x(i)>0)/M;
                    break
                end
            end
        end
    end
end


%%
figure
%subplot(1,4,1)
%for j=1:P
%histogram(t(:,j),'Normalization','probability','EdgeColor','none','FaceAlpha',.2)
%hold on
%end
%title('Simulated')
%axis([0 15 0 .5])


subplot(2,4,1)
cc=get(gca,'ColorOrder');
cc=(cc(1,:).^.3).*[0:.1:1]';
set(gca,'ColorOrder',cc);
title('Mean end time')
for l=1:2:Q
    hold on
    plot(drifts,nanmean(endTime(:,:,l),1),'LineWidth',2,'DisplayName',['\sigma=' num2str(noises(l))])
end
xlabel('Drift')
ylabel('Reaction step')
grid on
legend
subplot(2,4,2)
title('Std end time')
set(gca,'ColorOrder',cc);
for l=1:2:Q
    hold on
plot(drifts,nanstd(endTime(:,:,l),1),'LineWidth',2,'DisplayName',['\sigma=' num2str(noises(l))])
end
xlabel('Drift')
ylabel('Reaction step')
legend
grid on

subplot(2,4,3)
set(gca,'ColorOrder',cc);
    hold on
plot(drifts,correctRate(:,1:2:end),'LineWidth',2,'DisplayName',['\sigma=' num2str(noises(l))])
hold on
grid on
title('Accuracy')
ylabel('Correct decisions')

subplot(2,4,7)
%set(gca,'ColorOrder',cc);
    hold on
    for ratio=[6:-1:1]
        xind=2:min(floor(20/ratio)+1,11);
        yind=(xind-1)*ratio+1;
        plot(drifts(xind),correctRate(sub2ind(size(correctRate),xind,yind)),'LineWidth',2,'DisplayName',['\sigma/\mu=' num2str(ratio)])
    end
%     plot(drifts([2:4]),[correctRate(2,7) correctRate(3,13) correctRate(4,19) ],'LineWidth',2,'DisplayName',['\sigma/\mu=6'])
%    plot(drifts([2:5]),[correctRate(2,6) correctRate(3,11) correctRate(4,16) correctRate(5,21)],'LineWidth',2,'DisplayName',['\sigma/\mu=5'])
%    plot(drifts([2:6]),[correctRate(2,5) correctRate(3,9) correctRate(4,13) correctRate(5,17) correctRate(6,21)],'LineWidth',2,'DisplayName',['\sigma/\mu=4'])
%    plot(drifts([2:7]),[correctRate(2,4) correctRate(3,7) correctRate(4,10) correctRate(5,13) correctRate(6,16) correctRate(7,19)],'LineWidth',2,'DisplayName',['\sigma/\mu=3'])
%    plot(drifts([2:11]),[correctRate(2,3) correctRate(3,5) correctRate(4,7) correctRate(5,9) correctRate(6,11) correctRate(7,13) correctRate(8,15) correctRate(9,17) correctRate(10,19) correctRate(11,21)],'LineWidth',2,'DisplayName',['\sigma/\mu=2'])
%plot(drifts,diag(correctRate),'LineWidth',2,'DisplayName',['\sigma/\mu=1'])
%plot(drifts([3:2:11]),[correctRate(3,2) correctRate(5,3) correctRate(7,4) correctRate(9,5) correctRate(11,7)],'LineWidth',2,'DisplayName',['\sigma/\mu=.5'])
hold on
grid on
title('Accuracy')
ylabel('Correct decisions')
legend

subplot(2,4,4)
set(gca,'ColorOrder',cc);
hold on
for l=1:2:Q
    hold on
    plot(nanmean(endTime(:,:,l),1),correctRate(:,l),'LineWidth',2,'DisplayName',['\sigma=' num2str(noises(l))])
end
title('Acc vs RT')
ylabel('Accuracy (%)')
xlabel('Mean RT (steps)')
grid on