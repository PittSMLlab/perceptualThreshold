function datlogAnalysis(datlog)


%% Parse datlog to get HS, profiles, sent speeds
vRsent=datlog.TreadmillCommands.sent(:,1);
vLsent=datlog.TreadmillCommands.sent(:,2);
vSentT=datlog.TreadmillCommands.sent(:,4);
vRread=datlog.TreadmillCommands.read(:,1);
vLread=datlog.TreadmillCommands.read(:,2);
vReadT=datlog.TreadmillCommands.read(:,4);

vRload=datlog.speedprofile.velR;
vLload=datlog.speedprofile.velL;

RTOt=datlog.stepdata.RTOdata(:,4);
LTOt=datlog.stepdata.LTOdata(:,4);
RHSt=datlog.stepdata.RHSdata(:,4);
LHSt=datlog.stepdata.LHSdata(:,4);

try
    audiostart=datlog.audioCues.start;
    audiostop=datlog.audioCues.stop;
catch
    audiostart=[];
    audiostop=[];
end

if length(RHSt)>length(LHSt) %This means we started counting events at an RTO, and therefore RTOs mark new strides for the GUI and datlog
    pTOt=RTOt; %Primary event
    sTOt=LTOt; %Secondary event
    %vDload=[vRload(2:end)'-vLload(1:end-1)' 0; vRload(2:end)'-vLload(2:end)' 0];
    %vDload=vDload(:);
else %LTOs mark new strides
    pTOt=LTOt;
    sTOt=RTOt;
    %vDload=[vRload(1:end-1)'-vLload(2:end)' 0; vRload(2:end)'-vLload(2:end)' 0];
    %vDload=vDload(:);
end
    %aTOt=[pTOt'; sTOt'];
    %aTOt=aTOt(:);

vR=interp1(vReadT,vRread,pTOt,'nearest');
vL=interp1(vReadT,vLread,pTOt,'nearest');
vR=interp1(vSentT,vRsent,pTOt,'previous'); %Last sent speed BEFORE the event
vL=interp1(vSentT,vLsent,pTOt,'previous'); %Last sent speed BEFORE the event
vD=vR-vL;

inds=find([isnan(vRload(2:end)) & ~isnan(vRload(1:end-1))])-1;
pSize=vRload(inds)-vLload(inds);

Lpress=strcmp(datlog.addLog.keypress(:,1),'leftarrow') | strcmp(datlog.addLog.keypress(:,1),'numpad4');
LpressT=(cell2mat(datlog.addLog.keypress(Lpress,2))-datlog.framenumbers.data(1,2))*86400;
LpressT=(cell2mat(datlog.addLog.keypress(Lpress,2))); %New ver
Rpress=strcmp(datlog.addLog.keypress(:,1),'rightarrow') | strcmp(datlog.addLog.keypress(:,1),'numpad6');
RpressT=(cell2mat(datlog.addLog.keypress(Rpress,2))-datlog.framenumbers.data(1,2))*86400 ;
RpressT=(cell2mat(datlog.addLog.keypress(Rpress,2))) ;


vD_atLpress=interp1(vReadT,vRread-vLread,LpressT,'previous');
vD_atRpress=interp1(vReadT,vRread-vLread,RpressT,'previous');
%% Plot
pp=unique(pSize);
cmap=parula(length(pp)); %parula, jet, redbluecmap
 cmap=cmap*.8;
 
figure('Units','Normalized','OuterPosition',[0 0 1 1])
subplot(3,1,1)
hold on
%plot(vReadT,vRread-vLread)
%plot(pTOt,vD,'k') %Actual speeds
plot(vReadT,vRread-vLread,'k')
yy=interp1(vSentT,vRsent-vLsent,[0:.01:vSentT(end)],'previous');
plot([0:.01:vSentT(end)],yy,'g') %Actual speeds
set(gca,'XTick',sort([audiostart; audiostop],'ascend'))
aa=axis;
axis([0 aa(2:4)])
grid on
plot(pTOt(inds(inds<length(pTOt))),pSize(inds<length(pTOt)),'kx')
p2(1)=plot(RpressT,vD_atRpress,'o','Color',cmap(1,:).^.6,'MarkerFaceColor',cmap(1,:).^.6,'MarkerEdgeColor','none','MarkerSize',4);
p2(2)=plot(LpressT,vD_atLpress,'o','Color',cmap(end,:),'MarkerFaceColor',cmap(end,:),'MarkerEdgeColor','none','MarkerSize',4);


title('Trial time-course')
legend('Speeds read', 'Sent commands','Speed schedule','Trial begin','<- press','-> press')

subplot(3,4,5)
hold on
patch([0 23 23 0],[-400 -400 400 400],.7*ones(1,3),'FaceAlpha',.5,'EdgeColor','none')


for i=1:length(pp)
    data=vD(bsxfun(@plus,inds(pSize==pp(i)),[-5:30])');
    p1(i)=plot([-5:30],mean(data,2),'LineWidth',2,'Visible','off','Color',cmap(i,:));
    plot([-5:30],data','Color',p1(i).Color)  
end
axis([-5 30 -400 400])
title('Individual trials')
xlabel('Strides')
ylabel('Speed difference (mm/s)')
grid on

subplot(3,4,6)
hold on
patch([0 23 23 0],[-400 -400 400 400],.7*ones(1,3),'FaceAlpha',.5,'EdgeColor','none')
axis([-5 30 -400 400])
clear data2 data
for i=1:length(pp)
    data{i}=vD(bsxfun(@plus,inds(pSize==pp(i)),[-5:30])');
    %if size(data{i},2)>1
    %    aD=mean(data{i},1);
    %else
    %    aD=data{i}; %No reps
    %end
    plot([-5:30],mean(data{i},2),'LineWidth',2,'Color',p1(i).Color);
    if pp(i)~=0
    data2(:,:,i)=sign(-1*diff(data{i},1,2))==sign(data{i}(:,1:end-1));
    end
end
title('Average perfomance by perturbation size')
xlabel('Strides')
ylabel('Speed difference (mm/s)')
grid on

subplot(3,4,7)
hold on
bins=[-425:50:425];
vL1=hist(vD_atLpress,bins);
vR1=hist(vD_atRpress,bins);
v2=hist(vD,bins);
%histogram(vD_atLpress,[-351:20:351])
%histogram(vD_atRpress,[-351:20:351])

plot(bins,vL1./v2,'Color',p2(2).Color,'LineWidth',2)
plot(bins,vR1./v2,'Color',p2(1).Color,'LineWidth',2)
plot(bins(vL1==0 & vR1==0),log10(v2(vL1==0 & vR1==0)))
plot(bins,log10(v2),'k')
xlabel('Speed diff (mm/s)')
ylabel('Keypresses per stride')
legend('Leftarrow','Rightarrow','No response strides (log10)','Strides at speed (log10)')
title('Avg. keypress/stride vs. speed diff')

subplot(3,4,8)
histogram(vD(inds(pSize~=0)+22),[-150 -100 -75 -50 -30 -15 -8 0 8 15 30 50 75 100 150],'Normalization','pdf')
title('Histogram of final speed diff. excluding null trials')

subplot(3,4,9) %Response time: Time to first keypress
allPressT=sort([LpressT; RpressT]);
respTime=[];
pertSize=[];
hold on
for i=1:length(pp)
    indsAux=inds(pSize==pp(i));
    it=nan(size(indsAux));
    for j=1:length(indsAux)
        relEvent=pTOt(indsAux(j));
        aux=find(allPressT > (relEvent-1),1,'first');
        if ~isempty(aux) && (allPressT(aux)-relEvent)<30
            it(j)=allPressT(aux)-relEvent;
        end
    end
    respTime=[respTime; it];
    pertSize=[pertSize; pp(i)*ones(size(it))];
    plot(pp(i),it,'o','MarkerFaceColor',p1(i).Color,'MarkerEdgeColor','none','MarkerSize',4)
end

title('Reaction time [until first keypress]')
ylabel('Log-Time (s)')
xlabel('Perturbation speed (mm/s)')
xx=abs(pertSize(~isnan(respTime)));
yy=log(respTime(~isnan(respTime)));
tt=[xx ones(size(xx))]\yy;
if imag(tt)~=0
    error('Regression had imaginary part')
end
plot([-350:350],exp(tt(2))*exp(abs([-350:350])*tt(1)),'r');
%tt=[xx.^2 ones(size(xx))]\yy;
%plot([-350:350],exp(tt(2))*exp(abs([-350:350]).^2 *tt(1)));
%plot(xx,exp(yy),'.')
set(gca,'YScale','log','YTick',[.01 .1 1 10 100])
grid on
axis([-360 360 .01 30])
%set(gca,'Color',.8*ones(1,3))

subplot(3,4,10) %First keypress is in the right direction (accuracy)
allPressT=sort([LpressT; RpressT]);
auxLPT=[LpressT; 1e12];
auxRPT=[RpressT; 1e12];
hold on
clear it
for i=1:length(pp)
    indsAux=inds(pSize==pp(i));
    it{i,1}=zeros(size(indsAux));
    for j=1:length(indsAux)
        switch sign(pp(i))
            case 1
                aux1=find(auxRPT> pTOt(indsAux(j)),1,'first');
                aux2=find(auxLPT > pTOt(indsAux(j)),1,'first');
                aux=auxRPT(aux1)>auxLPT(aux2); %Pressed the right (->) button first
            case -1
                aux1=find(auxLPT > pTOt(indsAux(j)),1,'first');
                aux2=find(auxRPT > pTOt(indsAux(j)),1,'first');
                aux=auxLPT(aux1)<auxRPT(aux2); %Pressed the right (->) button first
            case 0
                aux1=find(allPressT > pTOt(indsAux(j)),1,'first');
                aux=isempty(aux1) || (allPressT(aux1)-pTOt(indsAux(j)))>30; %No presses
        end
        it{i}(j)=aux;
    end
    if sign(pp(i))==-1
        ph(i)=plot(pp(i),sum(~it{i})/length(it{i}),'o','MarkerFaceColor',p1(i).Color,'MarkerEdgeColor','none','MarkerSize',4);
    else
        ph(i)=plot(pp(i),sum(it{i})/length(it{i}),'o','MarkerFaceColor',p1(i).Color,'MarkerEdgeColor','none','MarkerSize',4);
    end
end
title('Accuracy of first keypress')
ylabel('%')
xlabel('Perturbation speed (mm/s)')
axis([-400 400 0 1])
grid on

subplot(3,4,11)
hold on
%Also: plot best-fit psychometric curve based on probability of subjects
%moving in the right direction on the first stride
auxY=cell2mat(it(pp~=0));
auxX=[];
clear ph
for i=1:length(pp)
    if pp(i)~=0
        ph(i)=plot(pp(i),sum(it{i})/length(it{i}),'o','MarkerFaceColor',p1(i).Color,'MarkerEdgeColor','none','MarkerSize',5);
        auxX=[auxX; pp(i)*ones(size(it{i}))];
    end
end
pp0=plot(auxX+5*(randn(size(auxY))), auxY,'k.'); %Indiv trials with some noise to see multiple equal responses
[p,~] = fitPsycho(auxX,auxY,'MLE');
[ff,gg]=psycho(p,auxX);
f=-sum(log(auxY.*(ff) + (1-auxY).*(1-ff))) %log-likelihood quantif.
xx=[-400:400];
pp1=plot(xx,psycho(p,xx),'k');
%[p,~] = fitPsycho(auxX,auxY,'MSE');
%[ff,gg]=psycho(p,auxX);
%f=-sum(log(auxY.*(ff) + (1-auxY).*(1-ff)))
%pp2=plot(xx,psycho(p,xx),'r');

%[p,peval] = fitPsycho(pp,cellfun(@nanmean,it),'MSE'); %MSE on the average
%data is the same as MSE on individual, binary, datapoints. MSE = bias^2 +
%std^2 at any given speed, and since the std is fixed, MSE is minimized by
%minimizing the bias. The bias, at any given speed, is by definition
%proportional to the distances to each individual point. Perhaps the result
%is different if we have a different number of samples at each datapoint.
%[ff,gg]=psycho(p,auxX);
%f=-sum(log(auxY.*(ff) + (1-auxY).*(1-ff)))
%pp3=plot(xx,psycho(p,xx));
%[p,~] = fitPsycho(auxX,auxY,'MAE'); %L1 norm minimzation just returns a
%very tight psychom function, as it is minimized by choosing a function
%that is exactly 1 at speeds where p>.5 and 0 at speeds where p<.5, not
%quite sure what would happen if the empirical data did not show a
%monotonic increase with p
%[ff,gg]=psycho(p,auxX);
%f=-sum(log(auxY.*(ff) + (1-auxY).*(1-ff)))
%pp4=plot(xx,psycho(p,xx));
legend([ph(1) pp0 pp1],{'Avg. data','Indiv. trials','MLE fit'},'Location','SouthEast')
title('Psychometric fit to initial ''<-'' responses')
xlabel('Speed diff (mm/s)')
ylabel('Prob. of pressing ''<-'' as first key')
set(gca,'YTick',[0 .25 .5 .75 1],'XTick',[-350 -250 -150 -75 0 75  150 250 350])
grid on

subplot(3,4,12)
%Also: plot steady-state dependence with initial speed
hold on
patch([-350 0 0], [-350 -350 0],[.7 0 0],'EdgeColor', 'none','FaceAlpha',.6)
patch([350 0 0 350], [0 0 -350 -350],.7*ones(1,3),'EdgeColor', 'none','FaceAlpha',.5)
patch([350 0 0], [350 350 0],[.7 0 0],'EdgeColor', 'none','FaceAlpha',.6)
patch([-350 0 0 -350], [0 0 350 350],.7*ones(1,3),'EdgeColor', 'none','FaceAlpha',.5)
text(100, -200, 'Overshoot','Color','k')
text(-250, -250, 'Wrong correction','Color',[.6 0 0])
allData=[];
for i=1:length(pp)
plot(pp(i), data{i}(29,:),'o','MarkerFaceColor',p1(i).Color,'MarkerEdgeColor','none','MarkerSize',4)
allData=[allData data{i}(29,:)];
end
pv=prctile(allData,[10,90]);
plot(350*[-1 1],pv(1)*[1 1],'k')
text(350,pv(1),[num2str(10) '%'])
plot(350*[-1 1],pv(2)*[1 1],'k')
text(350,pv(2),[num2str(90) '%'])
title('Final speed vs. initial speed')
xlabel('Initial speed (mm/s)')
ylabel('Final speed (mm/s) [+24 strides]') 
axis tight
grid on

end

