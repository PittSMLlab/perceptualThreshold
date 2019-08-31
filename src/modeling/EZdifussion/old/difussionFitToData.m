%% diffusion decision model 
%This script shows the power of hand-tuned diffusion model to explain
%perceptual data
%% load empirical data for comparison
loadAllDataIntoTable
pSize=[0 10 25 50 75 100 125 150 200 250 300 350];
P=length(pSize);
superSuperT.correctResponse=sign(superSuperT.initialResponse)==-sign(superSuperT.pertSize);
superSuperT.noResponse=isnan(superSuperT.reactionTime);
superSuperT.pertSize=abs(superSuperT.pertSize);
empMeanRS=nan(P,1);
empStdRS=nan(P,1);
empCorrectMeanRS=nan(P,1);
empCorrectStdRS=nan(P,1);
empIncorrectMeanRS=nan(P,1);
empIncorrectStdRS=nan(P,1);
for j=1:P
        dataT=superSuperT(superSuperT.pertSize==pSize(j) | superSuperT.pertSize==-pSize(j),:);
        empAcc(j,:)=dataT.correctResponse;
        %empMeanRS(j)=nanmean(dataT.reactionStride);
        %empStdRS(j)=nanstd(dataT.reactionStride);
        empMeanRS(j)=nanmean(dataT.reactionTime);
        empStdRS(j)=nanstd(dataT.reactionTime);
        empCorrectMeanRS(j)=nanmean(dataT.reactionTime(dataT.correctResponse));
        empCorrectStdRS(j)=nanstd(dataT.reactionTime(dataT.correctResponse));
        empIncorrectMeanRS(j)=nanmean(dataT.reactionTime(~dataT.correctResponse));
        empIncorrectStdRS(j)=nanstd(dataT.reactionTime(~dataT.correctResponse));
                empMedianRS(j)=nanmedian(dataT.reactionTime);
        empCorrectMedianRS(j)=nanmedian(dataT.reactionTime(dataT.correctResponse));
        empIncorrectMedianRS(j)=nanmedian(dataT.reactionTime(~dataT.correctResponse));
end

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
th=a/2; %Threshold, arbitrary scale
M=2e3; %Number of simulations for each parameter pair
drifts=[0,.05,.1,.2,.4,.8,1];
noises=1;
Q=length(noises); 
P=length(drifts);
%correctRate=zeros(P,Q);
correctResponse=false(M,P,Q);
endTime=nan(M,P,Q);
clockStep=.001;
thresholdCurve=th*ones(1,N);
%thresholdCurve=1.1*th*(1+(2*[1:N]/N).^.5);
for l=1:Q
    s=sqrt(noises(l));
    for j=1:P
        m=drifts(j);
        for k=1:M %Number of sims
            x=zeros(N,1);
            for i=2:N
                x(i)=x(i-1)+m*clockStep+s*randn*sqrt(clockStep);
                if abs(x(i))>=thresholdCurve(i)
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

%% Plot model results
figure;
subplot(4,1,2)
plot([1:N]*clockStep,thresholdCurve,'k')
hold on
plot([1:N]*clockStep,-thresholdCurve,'k')
for l=1:2:P
    subplot(4,1,2)
    hold on
ppp=plot(endTime(~isnan(endTime(:,l,1)),l,1),-thresholdCurve(round(endTime(~isnan(endTime(:,l,1)),l,1)/clockStep))'.*(-1).^(correctResponse(~isnan(endTime(:,l,1)),l,1)==1),'o');
subplot(4,1,1)
hold on
histogram(endTime(correctResponse(:,l,1),l,1),[0:.33:30],'FaceColor',ppp.Color,'FaceAlpha',.4,'EdgeColor','none','Normalization','probability')
axis([0 30 0 .15])
subplot(4,1,3)
hold on
histogram(endTime(~correctResponse(:,l,1),l,1),[0:.33:30],'FaceColor',ppp.Color,'FaceAlpha',.5,'EdgeColor','none','Normalization','probability')
axis([0 30 0 .15])
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
plot(pSize*f,empMeanRS,'o','MarkerSize',5,'LineWidth',2,'DisplayName','Empirical mean')
plot(pSize*f,empMedianRS,'d','MarkerSize',5,'LineWidth',2,'DisplayName','Empirical median')
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
    plot(drifts,nanmedian(correctEndTime(:,:,l),1),'LineWidth',2,'DisplayName',['\sigma^2=' num2str(noises(l))],'Color',get(ppp(end),'Color'));
end
for l=1:2:Q
    y=-2*th*(drifts+1e-9)/noises(l);
    mrt=th./(drifts+1e-9) .*(1-exp(y))./(1+exp(y));
    plot(drifts,mrt,'LineWidth',1,'DisplayName',['\sigma^2=' num2str(noises(l))],'Color',get(ppp(ceil(l/2)),'Color'))
end
xlabel('Drift')
ylabel('Reaction step')
grid on
legend(ppp)
plot(pSize*f,empCorrectMeanRS,'o','MarkerSize',5,'LineWidth',2,'DisplayName','Empirical mean correct')
plot(pSize*f,empIncorrectMeanRS,'ko','MarkerSize',5,'LineWidth',2,'DisplayName','Empirical mean incorrect')
plot(pSize*f,empCorrectMedianRS,'d','MarkerSize',5,'LineWidth',2,'DisplayName','Empirical median correct')
plot(pSize*f,empIncorrectMedianRS,'kd','MarkerSize',5,'LineWidth',2,'DisplayName','Empirical median incorrect')

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
plot(pSize*f,empCorrectStdRS,'.','MarkerSize',5,'LineWidth',2,'DisplayName','Empirical std correct')
plot(pSize*f,empIncorrectStdRS,'k.','MarkerSize',5,'LineWidth',2,'DisplayName','Empirical std incorrect')

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
    ppp=plot(nanmean(endTime(:,:,l),1),correctedRate(:,l),'LineWidth',2,'DisplayName',['\sigma^2=' num2str(noises(l))]);
    plot(th./drifts' .* (2*(1./(1+exp(-2*th*drifts')))-1),correctedRate(:,l),'Color',ppp.Color)
end
plot(empMeanRS,mean(empAcc,2),'o','MarkerSize',5,'LineWidth',2,'DisplayName','Model')
plot(empMedianRS,mean(empAcc,2),'d','MarkerSize',5,'LineWidth',2,'DisplayName','Model')
title('Acc vs RT')
ylabel('Accuracy (%)')
xlabel('Mean/median RT (steps)')
grid on

subplot(2,5,10)
set(gca,'ColorOrder',cc);
hold on
for l=1:Q
    hold on
    ppp=plot(nanstd(endTime(:,:,l),1)./nanmean(endTime(:,:,l),1),correctedRate(:,l),'LineWidth',2,'DisplayName',['\sigma^2=' num2str(noises(l))]);
    acc=[.501:.001:.999];
    y=log(1./acc-1);
    vrt_z2=(-2).*(2*y.*exp(y)-exp(2*y)+1)./((exp(y)+1).^2 .*y);
    mrt_z=(1-exp(y))./(1+exp(y));
    plot(sqrt(vrt_z2(2:end))./mrt_z(2:end),acc(2:end),'Color',ppp.Color) %This curve is interesting because it links CV(RT) to accuracy INDEPENDENTLY of parameter fits
end
plot(empStdRS./empMeanRS,mean(empAcc,2),'o','MarkerSize',5,'LineWidth',2,'DisplayName','Model')
title('Acc vs CV(RT)')
ylabel('Accuracy (%)')
xlabel('CV RT')
grid on

subplot(2,5,8)
title('RT histogram | correct')
hold on
for j=1:3:P
        dataT=superSuperT(superSuperT.pertSize==pSize(j) | superSuperT.pertSize==-pSize(j),:);
   histogram(dataT.reactionTime,[0:.5:15],'EdgeColor','none')
end
subplot(2,5,9)
title('RT histogram | incorrect')
