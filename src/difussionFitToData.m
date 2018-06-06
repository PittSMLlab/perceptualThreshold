%% diffusion decision model 
%This script shows the power of hand-tuned diffusion model to explain
%perceptual data
%% load empirical data for comparison
loadAllDataIntoTable
pSize=[0 10 25 50 75 100 125 150 200 250 300 350];
P=length(pSize);
for j=1:P
        dataT=superSuperT(superSuperT.pertSize==pSize(j) | superSuperT.pertSize==-pSize(j),:);
        empAcc(j,:)=sign(dataT.initialResponse)==-sign(dataT.pertSize);
        empMeanRS(j)=nanmean(dataT.reactionStride);
        empStdRS(j)=nanstd(dataT.reactionStride);
        empMeanRS(j)=nanmean(dataT.reactionTime);
        empStdRS(j)=nanstd(dataT.reactionTime);
end
superSuperT.correctResponse=sign(superSuperT.initialResponse)==-sign(superSuperT.pertSize);
superSuperT.noResponse=isnan(superSuperT.reactionTime);
superSuperT.pertSize=abs(superSuperT.pertSize);
[pSize,driftRate,threshold,delay]=fitEZ(superSuperT);
a=median(threshold);
f=nanmedian(driftRate./abs(pSize));

%Modified model:
%No se cómo hacerlo, pero la ideas es abandonar la noción de que a tiene q
%ser constante, y hacer un fit simultaneo de mean(RT), std(RT) y accuracy
%y=-threshold.*driftRate;
%acc=1./(1+exp(y));
%x=(1-exp(2*y)+2*y.*exp(y))./(1+exp(y)).^2;
%[pp]=polyfit(empMeanRS(2:end)',mean(empAcc(2:end,:),2)-1),1);
%% Define diffusion model and simulate:
N=3e4; %Simulation steps
th=a/3; %Threshold, arbitrary scale
M=2e3; %Number of simulations for each parameter pair
drifts=[0,.05,.1,.2,.4,.8,1];
noises=1;
Q=length(noises); 
P=length(drifts);
%correctRate=zeros(P,Q);
correctResponse=false(M,P,Q);
endTime=nan(M,P,Q);
clockStep=.001;
for l=1:Q
    s=sqrt(noises(l));
    for j=1:P
        m=drifts(j);
        for k=1:M %Number of sims
            x=zeros(N,1);
            for i=2:N
                x(i)=x(i-1)+m*clockStep+s*randn*sqrt(clockStep);
                if abs(x(i))>=th*(1+(.4*i*clockStep)).^.5
                    endTime(k,j,l)=(i)*clockStep;
                    correctResponse(k,j,l)=x(i)>0;
                    %correctRate(j,l)=correctRate(j,l)+(x(i)>0)/M;
                    break
                end
            end
        end
    end
end
correctRate=reshape(mean(correctResponse),P,Q);
correctedRate=reshape(sum(correctResponse)./sum(~isnan(endTime)),P,Q); %Excluding non-responses


%%
figure
%subplot(1,4,1)
%for j=1:P
%histogram(t(:,j),'Normalization','probability','EdgeColor','none','FaceAlpha',.2)
%hold on
%end
%title('Simulated')
%axis([0 15 0 .5])


subplot(2,5,1)
cc=get(gca,'ColorOrder');
cc=(cc(1,:).^.3).*[1:length(noises)]'/length(noises);
set(gca,'ColorOrder',cc(1:end,:));
title('Mean end time')
for l=1:Q
    hold on
    plot(drifts,nanmean(endTime(:,:,l),1)+mean(delay),'LineWidth',2,'DisplayName',['\sigma^2=' num2str(noises(l))])
end
xlabel('Drift')
ylabel('Reaction step')
plot(pSize*f,empMeanRS,'o','MarkerSize',5,'LineWidth',2,'DisplayName','Empirical std')
grid on
legend

subplot(2,5,6)
set(gca,'ColorOrder',cc(1:end,:));
title('Mean end time | correct')
correctEndTime=endTime;
correctEndTime(~correctResponse)=nan;
ppp=[];
for l=1:Q
    hold on
    %To Do: plot avg. reaction time of correct responses ONLY
    ppp(end+1)=plot(drifts,nanmean(correctEndTime(:,:,l),1),'LineWidth',2,'DisplayName',['\sigma^2=' num2str(noises(l))]);
end
for l=1:2:Q
    y=2*th*(drifts+1e-9)/noises(l);
    plot(drifts,th./(drifts+1e-9) .*(1-exp(-y))./(1+exp(-y)),'LineWidth',1,'DisplayName',['\sigma^2=' num2str(noises(l))],'Color',get(ppp(ceil(l/2)),'Color'))
end
xlabel('Drift')
ylabel('Reaction step')
grid on
legend(ppp)

subplot(2,5,2)
title('Std RT')
set(gca,'ColorOrder',cc);
for l=1:Q
    hold on
plot(drifts,nanstd(endTime(:,:,l),1),'LineWidth',2,'DisplayName',['\sigma^2=' num2str(noises(l))])
end
xlabel('Drift')
ylabel('Reaction step')
hold on
plot(pSize*f,empStdRS,'o','MarkerSize',5,'LineWidth',2,'DisplayName','Empirical std')
%legend
grid on

subplot(2,5,7)
title('Std end time | correct')
set(gca,'ColorOrder',cc);
ppp=[];
for l=1:Q
    hold on
ppp(end+1)=plot(drifts,nanstd(correctEndTime(:,:,l),1),'LineWidth',2,'DisplayName',['\sigma^2=' num2str(noises(l))]);
end
for l=1:Q
    hold on
    y=-2*th*(drifts+1e-3)/noises(l);
    vrt=(th*noises(l)./(drifts+1e-3).^3).*(2*y.*exp(y)-exp(2*y)+1)./(exp(y)+1).^2;
plot(drifts,sqrt(vrt),'LineWidth',1,'DisplayName',['\sigma^2=' num2str(noises(l))],'Color',get(ppp(l),'Color'))
end
xlabel('Drift')
ylabel('Reaction step')
%legend(ppp)
grid on

subplot(2,5,3)
set(gca,'ColorOrder',cc);
hold on
plot(drifts,squeeze(mean(isnan(endTime(:,:,1:end)))),'LineWidth',2,'DisplayName',['\sigma^2=' num2str(noises(l))])
grid on
ylabel('% non-response')
xlabel('Drift')
title('Non-response rate')


subplot(2,5,4)
set(gca,'ColorOrder',[cc;cc;get(gca,'ColorOrder');]);
    hold on
plot(drifts,correctedRate(:,1:end),'LineWidth',2,'DisplayName',['\sigma^2=' num2str(noises(l))])
%plot(drifts,correctRate(:,1:end),'LineWidth',2,'DisplayName',['\sigma^2=' num2str(noises(l))])
plot(drifts,1./(1+exp(-2*th*drifts./(noises(1:end))'))','LineWidth',1,'DisplayName',['\sigma^2=' num2str(noises(l))])
%plot(drifts,correctRate(:,1:2:end),'LineWidth',2,'DisplayName',['\sigma^2=' num2str(noises(l))])
hold on
grid on
title('Accuracy')
ylabel('Correct decisions')
plot(pSize*f,mean(empAcc,2),'o','MarkerSize',5,'LineWidth',2,'DisplayName','Empirical')

ppp=[];
    hold on
    for ratio=[10,5,4,2,1,.5]
        yind=nan(size(drifts));
        for i=1:length(drifts)
            aux=find(drifts(i)*ratio==noises);
            if ~isempty(aux)
                yind(i)=aux;
            end
        end
%         if ~all(isnan(yind))
%             xind=find(~isnan(yind));
%             yind=yind(xind);
%         ppp(end+1)=plot(drifts(xind),correctedRate(sub2ind(size(correctedRate),xind,yind)),'o','LineWidth',2,'DisplayName',['\sigma^2/\mu=' num2str(ratio)]);
%         plot([0 1],1./(1+exp(-2*th/ratio))*[1 1],'Color',get(ppp(end),'Color'))
%         end
    end

hold on
grid on
legend(ppp,'Location','SouthEast')

subplot(2,5,5)
set(gca,'ColorOrder',cc);
hold on
for l=1:Q
    hold on
    plot(nanmean(endTime(:,:,l),1),correctedRate(:,l),'LineWidth',2,'DisplayName',['\sigma^2=' num2str(noises(l))])
end
plot(empMeanRS,mean(empAcc,2),'o','MarkerSize',5,'LineWidth',2,'DisplayName','Model')
title('Acc vs RT')
ylabel('Accuracy (%)')
xlabel('Mean RT (steps)')
grid on