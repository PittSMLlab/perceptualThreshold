function fh=datlogAnalysis(datlog,goodOnly)
if nargin<2 || isempty(goodOnly)
    goodOnly=0;
end

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

trialStrides=isnan(vRload);
inds=find(trialStrides(2:end) & ~trialStrides(1:end-1)); %Last automated stride control: start of trial!
%inds=find([isnan(vRload(2:end)) & ~isnan(vRload(1:end-1))]);
pDuration=find(~trialStrides(inds(1)+1:end),1,'first');

Lpress=strcmp(datlog.addLog.keypress(:,1),'leftarrow') | strcmp(datlog.addLog.keypress(:,1),'numpad4') | strcmp(datlog.addLog.keypress(:,1),'pageup');
%LpressT=(cell2mat(datlog.addLog.keypress(Lpress,2))-datlog.framenumbers.data(1,2))*86400;
LpressT=(cell2mat(datlog.addLog.keypress(Lpress,2))); %New ver
Rpress=strcmp(datlog.addLog.keypress(:,1),'rightarrow') | strcmp(datlog.addLog.keypress(:,1),'numpad6') | strcmp(datlog.addLog.keypress(:,1),'pagedown');
%RpressT=(cell2mat(datlog.addLog.keypress(Rpress,2))-datlog.framenumbers.data(1,2))*86400 ;
RpressT=(cell2mat(datlog.addLog.keypress(Rpress,2))) ;

%% Find reaction times & accuracy of first keypress
[allPressT,sortIdxs]=sort([LpressT; RpressT]);
pressedKeys=[-1*ones(size(LpressT)); ones(size(RpressT))]; %R=1, L=-1
pressedKeys=pressedKeys(sortIdxs);
isItRpress=pressedKeys==1;

%Define some variables for each trial:
reactionTime=nan(size(inds));
reactionStride=nan(size(inds));
pertSize=vD(inds);
pertSign=sign(pertSize);
reactionSign=nan(size(inds)); %Positive if vR>vL
accurateReaction=nan(size(inds));
%goodTrial=zeros(size(inds));
pressTrial=nan(size(allPressT));
startCue=pTOt(inds);
endCue=pTOt(inds+pDuration-1);
audioStartCue=datlog.audioCues.start;
audioStopCue=datlog.audioCues.stop;
includedStrides=nan(length(inds),pDuration);

for i=1:length(inds)
        includedStrides(i,:)=inds(i) + [0:pDuration-1];
        relEvent=startCue(i); %Last TO where belt-speeds were under automated control
        relEvent2=endCue(i);
        aux=find(allPressT > (relEvent-1),1,'first');
        aux2=find(pTOt > allPressT(aux),1,'first');
        if ~isempty(aux) && (aux2-inds(i))<pDuration
            reactionTime(i)=allPressT(aux)-relEvent;
            reactionStride(i)=aux2-inds(i);
            reactionSign(i) = pressedKeys(aux);
            if (pertSign(i) == -1*reactionSign(i)) %|| pertSign(i)==0 %Correct choice!
                accurateReaction(i)=true;
            else
                accurateReaction(i)=false;
            end
        end
        aux=find((allPressT > relEvent) & (allPressT < relEvent2));
        pressTrial(aux)=i;
end

%% Filter trials & presses to consider:
if goodOnly==0
    mask=ones(size(inds)); %-> Override: Every trial is 'good'!
elseif goodOnly==-1
    mask=1-accurateReaction;
elseif goodOnly==1
    mask=accurateReaction;
end

    %Mask everything in trials:
    inds=inds(mask==1); %Keep only 'good' trials
    reactionTime=reactionTime(mask==1);
    reactionStride=reactionStride(mask==1);
    pertSize=vD(inds);
    pertSign=pertSign(mask==1);
    reactionSign=reactionSign(mask==1); %Positive if vR>vL
    accurateReaction=accurateReaction(mask==1);
    startCue=startCue(mask==1);
    endCue=endCue(mask==1);
    audioStartCue=audioStartCue(mask==1);
    audioStopCue=audioStopCue(mask==1);
    includedStrides=includedStrides(mask==1,:);
    
    %Fake trial number for presses outside a valid trial & expand goodTrial vector:
    pressTrial(isnan(pressTrial))=length(mask)+1; 
    mask(end+1)=0;
    
    %Filter presses to only those that happend during goodTrials
    LpressT=allPressT(pressedKeys==-1 & mask(pressTrial));
    RpressT=allPressT(pressedKeys==1 & mask(pressTrial));
    
    %Filter trialStrides to only those that happened during good trials:
    trialStrides=false(size(trialStrides));
    aux=bsxfun(@plus,inds,[0:pDuration-1]);
    trialStrides(aux(:))=true;
    trialStrides(~isnan(vRload))=false; %This shouldn't do anything, but just in case
    
    % Compute some derived quantities:
pSize=vRload(inds)-vLload(inds);
vD_atLpress=interp1(vReadT,vRread-vLread,LpressT,'previous');
vD_atRpress=interp1(vReadT,vRread-vLread,RpressT,'previous');


%% Do some response counting: (this can be improved using the previously defined variables instead, no need to do another for-loop)
accumRresponses=nan(size(pTOt));
accumLresponses=nan(size(pTOt));
firstResponse=nan(size(pTOt));
responseTime=nan(size(pTOt));
for i=1:length(pTOt)
    if trialStrides(i)
        aux1=find(allPressT> pTOt(i-1) & allPressT<= pTOt(i) ,1,'first');
        accumRresponses(i) = sum(RpressT> pTOt(i-1) & RpressT<= pTOt(i));
        accumLresponses(i) = sum(LpressT> pTOt(i-1) & LpressT<= pTOt(i));
        if ~isempty(aux1) %There was a press in that stride
            aux=1-isItRpress(aux1);
            responseTime(i)=allPressT(aux1)-pTOt(i-1); %Time in secs from prev pTOt, only relevant if this is the first stride in trial
        else
            aux=.5;
        end
        firstResponse(i)=aux; %L or R press
    end
end

%% Plot
pp=unique(pSize);
cmap=parula(length(pp)); %parula, jet, redbluecmap
cmap=cmap*.8;
 
fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
subplot(3,1,1)
hold on
%plot(vReadT,vRread-vLread,'k')
yy=interp1(vSentT,vRsent-vLsent,[0:.01:vSentT(end)],'previous');
plot([0:.01:vSentT(end)],yy,'g') %Actual speeds
set(gca,'XTick',sort([audiostart; audiostop],'ascend'))
aa=axis;
axis([0 aa(2:4)])
axis tight
grid on
plot(pTOt(inds(inds<length(pTOt))-1),pSize(inds<length(pTOt)),'kx')
p2(1)=plot(RpressT,vD_atRpress,'o','Color',cmap(1,:).^.6,'MarkerFaceColor',cmap(1,:).^.6,'MarkerEdgeColor','none','MarkerSize',4);
p2(2)=plot(LpressT,vD_atLpress,'o','Color',cmap(end,:),'MarkerFaceColor',cmap(end,:),'MarkerEdgeColor','none','MarkerSize',4);


title('Trial time-course')
legend('Sent commands','Trial begin','-> press','<- press')

subplot(3,4,5)
hold on
patch([0 pDuration-1 pDuration-1 0],[-400 -400 400 400],.7*ones(1,3),'FaceAlpha',.5,'EdgeColor','none')
rangeOfInterest=[-5:pDuration+5];
for i=1:length(pp)
    data{i}=vD(bsxfun(@plus,inds(pSize==pp(i)),rangeOfInterest)');
    p1(i)=plot(rangeOfInterest,mean(data{i},2),'LineWidth',2,'Visible','off','Color',cmap(i,:));
    plot(rangeOfInterest,data{i}','Color',p1(i).Color)  
end
axis([-5 30 -400 400])
title('Individual trials')
xlabel('Strides')
ylabel('Speed difference (mm/s)')
grid on


subplot(3,4,7)
bins=[-425:50:425];
bins=[-500 -351 -180 -110 -75 -50 -30 -15 -5 5 15 30 50 75 110 180 351 500];
binsC=bins(1:end-1)+.5*diff(bins);
%vL1=hist(vD_atLpress,bins);
%vR1=hist(vD_atRpress,bins);
aux=discretize(vD,bins);
auxL=nan(size(binsC));
auxR=nan(size(binsC));
v2=nan(size(binsC));
NR=nan(size(binsC));
for i=1:length(binsC)
    v2(i)= sum(trialStrides & aux==i); %Number of valid strides in bin
    if v2(i)>10 %Cut-off criteria for binning: with very few strides this becomes too noisy
    NR(i)= sum(accumLresponses(trialStrides & aux==i)==0 & accumRresponses(trialStrides & aux==i)==0); %Number of no-response strides in bin
    auxL(i)=sum(accumLresponses(trialStrides & aux==i));
    auxR(i)=sum(accumRresponses(trialStrides & aux==i));
    end
end

h1=gca;
h2=axes;
set(h2,'Position',get(h1,'Position'),'box','off','YAxisLocation','Right')
hold on
histogram(vD(trialStrides),bins)
%histogram(vD(trialStrides & accumRresponses==0 & accumLresponses==0),bins)

ylabel('Stride count')
legend({'Total strides'},'Location','NorthWest')
axis tight
aa=axis;
axis([-400 400 0 aa(4)])

axes(h1)
hold on
set(h1,'Color','none')
plot(binsC,auxL./v2,'Color',p2(2).Color,'LineWidth',2)
plot(binsC,auxR./v2,'Color',p2(1).Color,'LineWidth',2)
plot(binsC,NR./v2,'Color','k','LineWidth',2)
axis tight
aa=axis;
axis([-400 400 0 3.5])
xlabel('Speed diff (mm/s)')
ylabel('Keypresses per stride')
legend('<-','->','NR prob.')
title('Avg. keypress/stride vs. speed diff')
linkaxes([h1,h2],'x')

subplot(3,4,6)
% histogram(vD(inds+pDuration-1),bins)
% title('Histogram of final speed diff.')
% axis([-400 400 0 10])
% xlabel('Speed diff (mm/s)')
title('Keypress rates vs. stride #')

xx=[-3:pDuration+3];
auxL=accumLresponses(bsxfun(@plus,inds,xx));
auxL(isnan(auxL))=0;
auxL=mean(auxL,1);
plot(xx,auxL,'Color',p2(2).Color,'LineWidth',2)
hold on
auxR=accumRresponses(bsxfun(@plus,inds,xx));
auxR(isnan(auxR))=0;
auxR=mean(auxR,1);

plot(xx,auxR,'Color',p2(1).Color,'LineWidth',2)
aux=accumRresponses(bsxfun(@plus,inds,xx)) - accumLresponses(bsxfun(@plus,inds,xx));
aux(isnan(aux))=0;
aux=mean(abs(aux),1);
plot(xx,aux,'Color','k','LineWidth',1)
aux=accumRresponses(bsxfun(@plus,inds,xx)) - accumLresponses(bsxfun(@plus,inds,xx));
aux(isnan(aux))=0;
aux=mean((aux),1);
plot(xx,aux,'Color','g','LineWidth',1)
aux=accumRresponses(bsxfun(@plus,inds,xx)) + accumLresponses(bsxfun(@plus,inds,xx));
aux(isnan(aux))=0;
aux=nanmean(aux,1);
plot(xx,aux,'Color','r','LineWidth',1)
legend('<-','->','Abs. diff','Diff','All')
xlabel('Stride count from cue')
ylabel('Keypresses per stride')
grid on
ppA=patch([0 pDuration-1 pDuration-1 0],[-400 -400 400 400],.7*ones(1,3),'FaceAlpha',.5,'EdgeColor','none');
uistack(ppA,'bottom')
axis([-5 30 -.5 3.5])



subplot(3,4,8)
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
    plot(pp(i), data{i}(29,:),'o','MarkerFaceColor','none','MarkerEdgeColor',p1(i).Color,'MarkerSize',4)
    if size(data{i},2)>1
        plot(pp(i), mean(data{i}(29,:)),'o','MarkerFaceColor',p1(i).Color,'MarkerEdgeColor','none','MarkerSize',8)
    end
    allData=[allData data{i}(29,:)];
end
pv=prctile(allData,[16,84,50]);
plot(350*[-1 1],pv(1)*[1 1],'k')
text(350,pv(1),[num2str(16) '% =' num2str(pv(1),2)])
plot(350*[-1 1],pv(2)*[1 1],'k')
text(350,pv(2),[num2str(84) '% =' num2str(pv(2),2)])
plot(350*[-1 1],pv(3)*[1 1],'--k')
text(350,pv(3),[num2str(50) '% =' num2str(pv(3),2)])
title('Final speed vs. initial speed')
xlabel('Initial speed (mm/s)')
ylabel('Final speed (mm/s) [+24 strides]') 
%axis tight
axis([-350 350 -350 350])
grid on


subplot(3,4,9) %Response time: Time to first keypress
hold on
for i=1:length(pp)
    it=reactionTime(pertSize==pp(i));
    itA=it;
    itA(isnan(it))=30;
    plot(pp(i),itA,'o','MarkerFaceColor','None','MarkerEdgeColor',p1(i).Color,'MarkerSize',4)
    if length(it)>1
        plot(pp(i),nanmean(it),'o','MarkerFaceColor',p1(i).Color,'MarkerEdgeColor','None','MarkerSize',8)
        %plot(pp(i),exp(nanmean(log(it))),'o','MarkerFaceColor',p1(i).Color,'MarkerEdgeColor','None','MarkerSize',8)
    end
end

title('Reaction time [until first keypress]')
ylabel('Log-Time (s)')
xlabel('Perturbation speed (mm/s)')
xx=(pertSize(~isnan(reactionTime)));
yy=(reactionTime(~isnan(reactionTime)));
%Centered linear regression on log-space:
tt=[abs(xx) ones(size(xx))]\log(yy);
if imag(tt)~=0
    error('Regression had imaginary part')
end
%plot([-350:350],exp(tt(2)+abs([-350:350])*tt(1)),'r');
%Centered decaying exponential regression: (same model as before, but different
%weighting of errors) 
tt=fminunc(@(u) sum((yy-exp(-(abs(xx)*u(1) +u(2)))).^2),[1,0]);
plot([-350:350],exp(-(tt(2)+abs([-350:350])*tt(1))),'r');
%Uncentered regression:
%tt2=fminunc(@(u) sum((yy-exp(-(abs(xx-u(1))*u(2) +u(3)))).^2),[0 tt]);
%plot([-350:350],exp(-(tt2(3)+abs([-350:350]-tt2(1))*tt2(2))),'r--');

%inverse regression:
%tt=fminunc(@(u) sum((yy-1./((abs(xx)*u(1) +u(2)))).^2),[1,1]);
%plot([-350:350],1./((tt(2)+abs([-350:350])*tt(1))),'k');

set(gca,'YScale','log','YTick',[.1 1 10 30],'YTickLabel',{'.1','1','10','NR'})
grid on
axis([-360 360 .1 30])


auxLPT=[LpressT; 1e12];
auxRPT=[RpressT; 1e12];
% hold on
clear it
for i=1:length(pp)
    indsAux=inds(pSize==pp(i));
    it{i,1}=zeros(size(indsAux));
    for j=1:length(indsAux)
        aux1=find(auxRPT> pTOt(indsAux(j)),1,'first');
        aux2=find(auxLPT > pTOt(indsAux(j)),1,'first');
        if ~isempty(aux1) && (auxRPT(aux1)-pTOt(indsAux(j)))<30
            auxR=auxRPT(aux1);
        else
            auxR=1e6;
        end
        if ~isempty(aux2) && (auxLPT(aux2)-pTOt(indsAux(j)))<30
            auxL=auxLPT(aux2);
        else
            auxL=1e6;
        end
        if auxR>auxL; %Pressed the right (->) button first
            aux=1;
        elseif auxR==auxL %Tie OR no response
            aux=.5;
        else
            aux=0;
        end
        it{i}(j)=aux;
    end
end

subplot(3,4,10) %Reaction time vs accuracy (when you do respond)
hold on
auxX=[];
auxY=[];
for i=1:length(pp)
    if pp(i)~=0 && ~all(isnan(reactionTime(pertSize==pp(i))))%Excluding null perturbations from regression, as accuracy means something different for those.
        auxX=[auxX; nanmean(reactionTime(pertSize==pp(i)))];
        auxY=[auxY; nanmean(accurateReaction(pertSize==pp(i)))];
    
        plot(nanmean(reactionTime(pertSize==pp(i))),nanmean(accurateReaction(pertSize==pp(i))),'o','MarkerFaceColor',p1(i).Color,'MarkerEdgeColor','none','MarkerSize',5)
    end
end
p=polyfit(auxX,auxY,1);
r=corrcoef(auxX,auxY);
[r,pp]=corrcoef(auxX,auxY);
plot(1:10,[1:10]*p(1)+p(2),'r')
text(6,.9,['r=' num2str(r(1,2))])
hold off
xlabel('Avg. reaction time (s)')
ylabel('Avg. accuracy (%)')
title('Accuracy vs. reaction time')
grid on

subplot(3,4,11) %Psychometric curve of first keypress per trial being right-ward
hold on
auxY=[];
auxX=[];
clear ph
for i=1:length(pp)
    %if pp(i)~=0
        ph(i)=plot(pp(i),nanmean(reactionSign(~isnan(reactionTime) & pertSize==pp(i))==-1),'o','MarkerFaceColor',p1(i).Color,'MarkerEdgeColor','none','MarkerSize',5);
        auxY(i)=sum(isnan(reactionTime(pertSize==pp(i))))/sum(pertSize==pp(i)); %No reactions / all trials
        auxYY(i)=sum(accurateReaction(pertSize==pp(i))==1)/sum(~isnan(accurateReaction(pertSize==pp(i)))); %accurate / all reactions
        %auxX=[auxX; pp(i)*ones(size(it{i}))];
    %end
end
try
%[p,~] = fitPsycho(pertSize(~isnan(reactionTime)),reactionSign(~isnan(reactionTime))==-1,'MSE'); %Fitting only actual responses
[p2,~] = fitPsychoGauss(pertSize(~isnan(reactionTime)),reactionSign(~isnan(reactionTime))==-1,'MLE'); %Fitting only actual responses
catch
    warning('Could not fit psycho!')
    p=[0,1];
end
%[p,~] = fitPsycho(auxX,auxY,'MLEsat');
xx=[-400:400];
%pp1=plot(xx,psycho(p,xx),'k');
pp4=plot(xx,psychoGauss(p2,xx),'k');
[aa,ii]=sort(pertSize(~isnan(reactionTime)));
aux=reactionSign(~isnan(reactionTime));
MM=100;
aux=[ones(MM,1); aux(ii);-1*ones(MM,1)]; %Sorted by perturbation size
aa=[ones(MM,1)*-500 ; aa; ones(MM,1)*500];
leftCumProb=cumsum(aux==1)./cumsum(abs(aux)); %Leftward choices  over all choices, for more negative perturbations than the current one
rightCumProb=(sum(aux==-1)-cumsum(aux==-1))./(sum(abs(aux))-cumsum(abs(aux))) ; %Rightward over all, for more positive perturbations
k1=find(leftCumProb>=rightCumProb,1,'last');
k2=find(leftCumProb<=rightCumProb,1,'first');
%pp2=plot(pp,auxYY,'rx');
pp3=plot(pp,auxY,'b');
legend([ph(1) pp4 pp3],{'Empiric % of ''<-''','MLE of ''<-''','% of NR'},'Location','SouthEast')
title('Psychometric fit to initial keypress')
xlabel('Speed diff (mm/s)')
ylabel('Prob.')
set(gca,'YTick',[0:.1:1],'XTick',[-350 -250 -150 -75 0 75  150 250 350])
grid on
%text(-300, .75,['\mu=' num2str(p(1))])
%text(-300, .65,['\sigma=' num2str(p(2))])
text(-300, .95,['\mu=' num2str(p2(1))])
text(-300, .85,['\sigma=' num2str(p2(2))])
try
text(-300, .75,['th=' num2str(p2(3))])
end
text(-300,.65,['m=' num2str(aa(k1))])
plot(aa(k1)*[1 1],[0 1],'k--')
plot(aa(k2)*[1 1],[0 1],'k--')

subplot(3,4,12) %Psychometric curve of first keypress per STRIDE being right-ward
hold on
trialStrides2=zeros(size(trialStrides));
N=1;
aux=bsxfun(@plus,inds,[0:N:pDuration-1]);
trialStrides2(aux(:))=1;
trialStrides2=trialStrides2 & trialStrides;
auxY=firstResponse(trialStrides2(:));
auxX=vD(trialStrides2(:));
clear ph
for i=1:length(inds)
    ph(i)=plot(vD(inds(i)),firstResponse(inds(i)+1),'o','MarkerFaceColor',p1((vD(inds(i))==pp)).Color,'MarkerEdgeColor','none','MarkerSize',5);
end
pp0=plot(auxX+5*(randn(size(auxY))), auxY,'.','Color',.7*ones(1,3)); %Indiv trials with some noise to see multiple equal responses
try
[p,~] = fitPsycho(auxX(auxY~=.5),auxY(auxY~=.5),'MLE');
catch
    warning('Could not fit psycho!')
    p=[0,1e6];
end
%[p,~] = fitPsycho(auxX(auxY~=.5),auxY(auxY~=.5),'MLEsat');
[ff,gg]=psycho(p,auxX);
xx=[-400:400];
pp1=plot(xx,psycho(p,xx),'k');
legend([ph(1) pp0 pp1],{'First stride','All strides','MLE fit'},'Location','SouthEast')
title('Psychometric fit to first keypress in all strides')
xlabel('Speed diff (mm/s)')
ylabel('Prob. of pressing ''<-'' as first key')
set(gca,'YTick',[0:.1:1],'XTick',[-350 -250 -150 -75 0 75  150 250 350])
grid on
text(-300, .75,['\mu=' num2str(p(1))])
text(-300, .65,['\sigma=' num2str(p(2))])
try
text(-300, .55,['th=' num2str(p(3))])
end




end

