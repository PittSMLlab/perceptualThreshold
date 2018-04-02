function fh=accuracyPlots(trialData,goodOnly)
if nargin<2 || isempty(goodOnly)
    goodOnly=0;
end
%%
B=findgroups(trialData.pertSize); %pertSize>0 means vR>vL
pp=unique(trialData.pertSize);
cmap=parula(length(pp)); %parula, jet, redbluecmap
cmap=cmap*.8;
 
fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
Q=6;
%% First plot: overall accuracy bars
subplot(2,Q,1)
nullTrials=trialData.pertSize==0;
correctResponses=trialData.initialResponse==-sign(trialData.pertSize) & ~nullTrials; %Negative response means LEFT IS SLOW (RIGHT IS FAST) choice
nonResponse=isnan(trialData.initialResponse) & ~nullTrials;
incorrectResponses=trialData.initialResponse==sign(trialData.pertSize) & ~nullTrials;
d=[sum(correctResponses) sum(incorrectResponses)  sum(nonResponse)]/sum(~nullTrials);
dN=[sum(nonResponse) sum(nullTrials)]/numel(nullTrials);
d2=[sum(correctResponses & (trialData.pertSize<0))/sum(trialData.pertSize<0)];
d3=sum(correctResponses & (trialData.pertSize>0))/sum(trialData.pertSize>0);
bar([d,dN])
set(gca,'XTickLabel',{'OK','BAD','NR','NR all','Null'})
ylabel('% of total responses')

subplot(2,Q,1+Q)
b=bar([d2,d(:,1),d3]');
set(gca,'XTickLabel',{'vL>vR','ALL','vR>vL'})
ylabel('% CORRECT')
aa=axis;
axis([aa(1:2) .5 1])

%% Second plot: % <- choices as function of perturbation size
subplot(2,Q,2)
S=splitapply(@(x) (sum(x==-1)+.5*sum(isnan(x)))/length(x),trialData.initialResponse,B); %Counting LEFT IS SLOW choices plus HALF of no response
scatter(pp,S,50,pp,'filled')
grid on
ylabel('% ''<-'' (left is slow) responses') 
xlabel('vL>vR     PERTURBATION (mm/s)      vL<vR')
%TODO: add std

subplot(2,Q,2+Q)
B2=findgroups(abs(trialData.pertSize)); %pertSize>0 means vR>vL
pp2=unique(abs(trialData.pertSize));
S2=splitapply(@(x) (sum(x))/length(x),correctResponses+.5*nonResponse,B2); %Counting LEFT IS SLOW choices + half of NR
S2(1)=NaN;
scatter(pp2,S2,80,.4*ones(1,3),'filled')
grid on
ylabel('% CORRECT') 
xlabel('ABS SPEED PERTURBATION (mm/s)')
axis([0 400 .5 1])
hold on 
%Adding vR>vL and vL<vR overlapping
scatter(pp(pp>0),S(pp>0),20,cmap(end,:),'filled')
scatter(abs(pp(pp<0)),1-S(pp<0),20,cmap(1,:),'filled')
%TODO: add STD
%% Third plot: reaction times
rt=trialData.reactionTime;
%rt=trialData.reactionStride;
%fun=@nanmean;
fun=@nanmedian;
subplot(2,Q,4)
%scatter(trialData.pertSize,rt,10,zeros(1,3))
%hold on
%scatter(pp,S3,50,pp)
S4=splitapply(fun,rt,B);
scatter(pp,S4,70,pp,'filled')
grid on
ylabel('Reaction log-time (s))') 
xlabel('vL>vR      PERTURBATION         vL<vR')
legend({'Indiv. trials','Median'},'Location','NorthWest')
set(gca,'YScale','log')

subplot(2,Q,4+Q)
S=splitapply(fun,rt,B); 
S2=splitapply(fun,rt,B2); 
scatter(pp2,S2,80,.4*ones(1,3),'filled')
grid on
ylabel('Reaction log-time (s)') 
xlabel('ABS SPEED PERTURBATION (mm/s)')
hold on 
%Adding vR>vL and vL<vR overlapping
scatter(pp(pp>0),S(pp>0),20,cmap(end,:),'filled')
scatter(abs(pp(pp<0)),S(pp<0),20,cmap(1,:),'filled')
%TODO: add STD
set(gca,'YScale','log')

%% Fourth plot: steady-state as function of perturbation size
subplot(2,Q,3)
%scatter(trialData.pertSize,trialData.lastSpeedDiff,10,zeros(1,3))
%hold on
S4=splitapply(@median,trialData.lastSpeedDiff,B);
scatter(pp,S4,70,pp,'filled')
grid on
ylabel('Final speed (mm/s)') 
xlabel('vL>vR      PERTURBATION         vL<vR')
legend({'Indiv. trials','Median'},'Location','NorthWest')
axis([-360 360 -150 150])

subplot(2,Q,3+Q)
S=splitapply(@nanmedian,trialData.lastSpeedDiff,B); 
S2=splitapply(@nanmedian,trialData.lastSpeedDiff .* sign(trialData.pertSize),B2); 
scatter(pp2,S2,80,.4*ones(1,3),'filled')
hold on
grid on
ylabel('Final speed (mm/s)') 
xlabel('ABS SPEED PERTURBATION (mm/s)')
scatter(pp(pp>0),S(pp>0),20,cmap(end,:),'filled')
scatter(abs(pp(pp<0)),-S(pp<0),20,cmap(1,:),'filled')

%% Fifth: Reaction time vs. accuracy
subplot(2,Q,5)
M=min(20,round(size(trialData,1)/10));
S=splitapply(@nanmedian,trialData.reactionTime,B);
T=splitapply(@(x) (sum(x))/length(x),correctResponses,B);
T(12)=NaN;
[x,idx]=sort(trialData.reactionTime,'ascend');
y=correctResponses(idx);
y=conv(y,ones(M,1)/M,'same');
y(1:floor(M/2))=NaN;
y(end-ceil(M/2):end)=NaN;
plot(x,y,'k')
hold on
scatter(S,T,70,pp,'filled')
grid on
set(gca,'XScale','log')
xlabel('Reaction log-time (s)')
ylabel('% CORRECT')

subplot(2,Q,5+Q)
S=splitapply(@nanmedian,trialData.reactionTime,B2);
T=splitapply(@(x) (sum(x))/length(x),correctResponses,B2);
T(1)=NaN;
rt=trialData.reactionTime(trialData.pertSize>0);
[x,idx]=sort(rt,'ascend');
y=correctResponses(trialData.pertSize>0);
y=y(idx);
y=conv(y,ones(M,1)/M,'same');
y(1:floor(M/2))=NaN;
y(end-ceil(M/2):end)=NaN;
plot(x,y,'Color',cmap(end,:))
hold on
rt=trialData.reactionTime(trialData.pertSize<0);
[x,idx]=sort(rt,'ascend');
y=correctResponses(trialData.pertSize<0);
y=y(idx);
y=conv(y,ones(M,1)/M,'same');
y(1:floor(M/2))=NaN;
y(end-ceil(M/2):end)=NaN;
plot(x,y,'Color',cmap(1,:))
scatter(S,T,70,.4*ones(1,3),'filled')
grid on
set(gca,'XScale','log')
xlabel('Reaction log-time (s)')
ylabel('% CORRECT')

%% Sixth plot: clicking rates
subplot(2,Q,6)
fun=@nanmean;
fun=@nanmedian;
S=splitapply(fun, trialData.Lclicks+trialData.Rclicks,B);
scatter(pp,S,70,.4*ones(1,3),'filled')
hold on
S2=splitapply(fun, trialData.Lclicks-trialData.Rclicks,B);
scatter(pp,abs(S2),70,pp,'filled')
S3=splitapply(fun, (trialData.Rclicks),B);
scatter(pp,S3,20,cmap(1,:),'filled')
S4=splitapply(fun, (trialData.Lclicks),B);
scatter(pp,S4,20,cmap(end,:),'filled')
grid on

subplot(2,Q,Q+6)
%TODO: plot correct presses per trial vs. incorrect presses per trial, as
%function of abs(pertSize)

%% Seventh plot: reaction time vs accuracy

% %% Plot
% 
% 
% 
% subplot(3,4,5)
% % steady-state dependence with initial speed
% hold on
% patch([-350 0 0], [-350 -350 0],[.7 0 0],'EdgeColor', 'none','FaceAlpha',.6)
% patch([350 0 0 350], [0 0 -350 -350],.7*ones(1,3),'EdgeColor', 'none','FaceAlpha',.5)
% patch([350 0 0], [350 350 0],[.7 0 0],'EdgeColor', 'none','FaceAlpha',.6)
% patch([-350 0 0 -350], [0 0 350 350],.7*ones(1,3),'EdgeColor', 'none','FaceAlpha',.5)
% text(100, -200, 'Overshoot','Color','k')
% text(-250, -250, 'Wrong correction','Color',[.6 0 0])
% allData=[];
% allV=[];
% for i=1:length(pp)
%     plot(pp(i), data{i}(pDuration+4,:),'o','MarkerFaceColor','none','MarkerEdgeColor',p1(i).Color,'MarkerSize',4)
%     if size(data{i},2)>1
%         plot(pp(i), median(data{i}(pDuration+4,:)),'o','MarkerFaceColor',p1(i).Color,'MarkerEdgeColor','none','MarkerSize',8)
%     end
%     allData=[allData data{i}(pDuration+4,:)];
%     allV=[allV pp(i)*ones(size(data{i}(pDuration+4,:)))];
% end
% pv=prctile(allData,[16,84,50]);
% plot(350*[-1 1],pv(1)*[1 1],'k')
% text(350,pv(1),[num2str(16) '% =' num2str(pv(1),3)])
% plot(350*[-1 1],pv(2)*[1 1],'k')
% text(350,pv(2),[num2str(84) '% =' num2str(pv(2),3)])
% plot(350*[-1 1],pv(3)*[1 1],'--k')
% text(350,pv(3),[num2str(50) '% =' num2str(pv(3),3)])
% papa=plot([-350 350],[-350 350],'k','LineWidth',2);
% uistack(papa,'bottom')
% text(30,300,['No change line ->'])
% title('Final speed vs. initial speed')
% xlabel('Perturbation speed (mm/s)')
% ylabel('Final speed (mm/s) [+24 strides]') 
% %axis tight
% axis([-350 350 -350 350])
% grid on
% ppP=polyfit(allV,allData,1);
% plot([-350 350],[-350 350]*ppP(1)+ppP(2),'r')
% 
% 
% 
% subplot(3,4,10) %Response time: Time to first keypress
% hold on
% for i=1:length(pp)
%     it=reactionTime(pertSize==pp(i));
%     itA=it;
%     itA(isnan(it))=30;
%     plot(pp(i),itA,'o','MarkerFaceColor','None','MarkerEdgeColor',p1(i).Color,'MarkerSize',4)
%     if length(it)>1
%         plot(pp(i),nanmedian(it),'o','MarkerFaceColor',p1(i).Color,'MarkerEdgeColor','None','MarkerSize',8)
%         %plot(pp(i),exp(nanmedian(log(it))),'o','MarkerFaceColor',p1(i).Color,'MarkerEdgeColor','None','MarkerSize',8)
%     end
% end
% 
% title('Reaction time [until first keypress]')
% ylabel('Log-Time (s)')
% xlabel('Perturbation speed (mm/s)')
% xx=(pertSize(~isnan(reactionTime)));
% yy=(reactionTime(~isnan(reactionTime)));
% %Centered linear regression on log-space:
% tt=[abs(xx) ones(size(xx))]\log(yy);
% if imag(tt)~=0
%     error('Regression had imaginary part')
% end
% %plot([-350:350],exp(tt(2)+abs([-350:350])*tt(1)),'r');
% %Centered decaying exponential regression: (same model as before, but different
% %weighting of errors) 
% tt=fminunc(@(u) sum((yy-exp(-(abs(xx)*u(1) +u(2)))-u(3)).^2),[.01,2,1]);
% plot([-350:350],tt(3)+exp(-(tt(2)+abs([-350:350])*tt(1))),'r');
% %Uncentered regression:
% %tt2=fminunc(@(u) sum((yy-exp(-(abs(xx-u(1))*u(2) +u(3)))).^2),[0 tt]);
% %plot([-350:350],exp(-(tt2(3)+abs([-350:350]-tt2(1))*tt2(2))),'r--');
% 
% %inverse regression:
% %tt=fminunc(@(u) sum((yy-(u(3) + 1./(abs(xx)*u(1) +u(2)+.01))).^2),[.01,.01,1]);
% %plot([-350:350],tt(3)+1./((tt(2)+abs([-350:350])*tt(1)+.01)),'k');
% 
% set(gca,'YScale','log','YTick',[.1 1 10 30],'YTickLabel',{'.1','1','10','NR'})
% grid on
% axis([-360 360 .1 30])
% text(-150,.3,['t~' num2str(tt(3),2) '+' num2str(exp(-tt(2)),2) '*e^{-' num2str(tt(1),2) '|\Delta V|}'])
% 
% 
% auxLPT=[LpressT; 1e12];
% auxRPT=[RpressT; 1e12];
% % hold on
% clear it
% for i=1:length(pp)
%     indsAux=inds(pSize==pp(i));
%     it{i,1}=zeros(size(indsAux));
%     for j=1:length(indsAux)
%         aux1=find(auxRPT> pTOt(indsAux(j)),1,'first');
%         aux2=find(auxLPT > pTOt(indsAux(j)),1,'first');
%         if ~isempty(aux1) && (auxRPT(aux1)-pTOt(indsAux(j)))<30
%             auxR=auxRPT(aux1);
%         else
%             auxR=1e6;
%         end
%         if ~isempty(aux2) && (auxLPT(aux2)-pTOt(indsAux(j)))<30
%             auxL=auxLPT(aux2);
%         else
%             auxL=1e6;
%         end
%         if auxR>auxL %Pressed the right (->) button first
%             aux=1;
%         elseif auxR==auxL %Tie OR no response
%             aux=.5;
%         else
%             aux=0;
%         end
%         it{i}(j)=aux;
%     end
% end
% 
% subplot(3,4,11) %Reaction time vs accuracy (when you do respond)
% hold on
% auxX=[];
% auxY=[];
% allAuxX=[];
% allAuxY=[];
% allAuxZ=[];
% XX=[0:.1:10];
% YY=[0:.01:1]';
% VV(:,1)=reshape(repmat(XX,length(YY),1),[length(XX)*length(YY),1]);
% VV(:,2)=reshape(repmat(YY,1,length(XX)),[length(XX)*length(YY),1]);
% VVV=[VV.^2 2*VV(:,1).*VV(:,2)];
% for i=1:length(pp)
%     if pp(i)~=0 && ~all(isnan(reactionTime(pertSize==pp(i))))%Excluding null perturbations from regression, as accuracy means something different for those.
%         relX=reactionTime(pertSize==pp(i) & ~isnan(reactionTime));
%         relY=accurateReaction(pertSize==pp(i) & ~isnan(reactionTime));
%         auxX=[auxX; nanmean(relX)];
%         allAuxX=[allAuxX; relX(~isnan(relX))];
%         auxY=[auxY; nanmean(relY)];
%         allAuxY=[allAuxY; relY(~isnan(relX))];
%         allAuxZ=[allAuxZ; pp(i)*ones(size(relX(~isnan(relX))))]; %Storing perturbation sizes
%         plot(nanmean(relX),nanmean(relY),'o','MarkerFaceColor',p1(i).Color,'MarkerEdgeColor','none','MarkerSize',5)
%         if length(relX)>1
%         CC=cov(relX,relY);
%         if CC(2,2)==0
%             CC(2,2)=.001;
%         end
%         %VVV=bsxfun(@minus,VV,[nanmean(relX) nanmean(relY)]); %Subtracting mean
%         %VVV=[VVV.^2 2*VVV(:,1).*VVV(:,2)];
%         if CC(2,1)>0 %Only plotting covariances that are positive
%             %CCC=pinv(CC);
%             %CCC=[CCC(1,1),CCC(2,2),CCC(2,1)];
%             %Z=reshape(CCC*(VVV'),[length(YY),length(XX)]); 
%             %contour(XX,YY,Z,[1 1]*1e-1,'Color',p1(i).Color)
%         end
%         end
%     end
% end
% %save(['timeVsAccuracy' num2str(now) ]) 
% NN=25;
% %NN=51;
% [~,inds2]=sort(allAuxX);
% %plot(allAuxX,allAuxY,'o','MarkerFaceColor',.7*[1,1,1],'MarkerEdgeColor','none','MarkerSize',4)
% edges=[0:.5:5,6,7,8,11,30];
% [binAssignment]=discretize(allAuxX,edges);
% auxH=nan(length(edges)-1,1);
% auxC=nan(length(edges)-1,1);
% for i=[1:length(edges)-1]
% auxH(i)=mean(allAuxY(binAssignment==i));
% auxC(i)=sum(binAssignment==i);
% end
% xx=edges(1:end-1)+.5*diff(edges);
% b=bar(edges,[auxH;0],'histc');
% set(b,'FaceAlpha',.3,'EdgeColor','None')
% text(xx-.4,.1*ones(size(auxC))+.05*(-1).^[1:length(auxC)]',num2str(auxC),'FontSize',6)
% try 
% ppp=plot(allAuxX(inds2),conv([ones((NN-1)/2,1);allAuxY(inds2);.5*ones((NN-1)/2,1)],ones(NN,1)/NN,'valid'),'k','LineWidth',1)
% uistack(ppp,'bottom')
% catch
%     warning('Could not smooth time vs. accuracy data')
% end
% uistack(b,'bottom')
% %r=corrcoef(auxX,auxY);
% %[r,pv]=corr(auxX,auxY); %Linear corr
% %text(8,.9,['r=' num2str(r) ' ,p=' num2str(pv,2)],'FontWeight','bold','Color','r')
% %[r,pv]=corr(auxX,auxY,'type','Spearman'); %Spearman rank corr
% %text(8,.8,['sp=' num2str(r) ' ,p=' num2str(pv,2)],'FontWeight','bold','Color','r')
% %p=polyfit(auxX,auxY,1);
% %plot(1:10,[1:10]*p(1)+p(2),'r','LineWidth',2)
% [r,pv]=corr(allAuxX,allAuxY); %Linear corr
% text(8,.75,['r=' num2str(r) ' ,p=' num2str(pv,2)],'Color','r')
% [r,pv]=corr(allAuxX,allAuxY,'type','Spearman'); %Spearman rank corr
% text(8,.7,['sp=' num2str(r) ' ,p=' num2str(pv,2)],'Color','r')
% p=polyfit(allAuxX,allAuxY,1);
% ppp1=plot(1:10,[1:10]*p(1)+p(2),'r','LineWidth',2);
% [r,pv]=corr(log(allAuxX),allAuxY); %Linear corr
% text(8,.4,['r=' num2str(r) ' ,p=' num2str(pv,2)],'Color','k')
% [r,pv]=corr(log(allAuxX),allAuxY,'type','Spearman'); %Spearman rank corr
% text(8,.3,['sp=' num2str(r) ' ,p=' num2str(pv,2)],'Color','k')
% p=polyfit(log(allAuxX),allAuxY,1);
% ppp2=plot(1:20,log([1:20])*p(1)+p(2),'k','LineWidth',2);
% p=fminunc(@(p) sum((allAuxY-p(2)*exp(allAuxX*p(1))-p(3)).^2),[-1 1 .5]); %Need to replace by the MLE estimator instead of MSE
% ppp2=plot(1:20,exp([1:20]*p(1))*p(2)+p(3),'b','LineWidth',2);
% 
% hold off
% xlabel('Avg. reaction time (s)')
% ylabel('Avg. accuracy (%)')
% title('Accuracy vs. reaction time')
% legend([ppp,ppp1,ppp2],{['Running avg. binw=' num2str(NN)],'Linear Regression','Exp. Regression'})
% axis([0 30 0 1])
% grid on
% 
% subplot(3,4,9) %Psychometric curve of first keypress per trial being right-ward
% hold on
% auxY=[];
% auxX=[];
% clear ph
% for i=1:length(pp)
%     %if pp(i)~=0
%         ph(i)=plot(pp(i),nanmean(reactionSign(~isnan(reactionTime) & pertSize==pp(i))==-1),'o','MarkerFaceColor',p1(i).Color,'MarkerEdgeColor','none','MarkerSize',5);
%         auxY(i)=sum(isnan(reactionTime(pertSize==pp(i))))/sum(pertSize==pp(i)); %No reactions / all trials
%         auxYY(i)=sum(accurateReaction(pertSize==pp(i))==1)/sum(~isnan(accurateReaction(pertSize==pp(i)))); %accurate / all reactions
%         %auxX=[auxX; pp(i)*ones(size(it{i}))];
%     %end
% end
% try
%     if psychoLapse==0
%  %Fitting only actual responses
% %[p2,~] = fitPsychoGauss(pertSize(~isnan(reactionTime)),reactionSign(~isna[p,~] = fitPsycho(pertSize(~isnan(reactionTime)),reactionSign(~isnan(reactionTime))==-1,'MLE');n(reactionTime))==-1,'MLE'); %Fitting only actual responses
%     else
%         [p,~] = fitPsycho(pertSize(~isnan(reactionTime)),reactionSign(~isnan(reactionTime))==-1,'MLEsat'); %Fitting only actual responses
%     end
% catch
%     warning('Could not fit psycho!')
%     p=[0,1];
% end
% xx=[-400:400];
% pp1=plot(xx,psycho(p,xx),'k');
% %pp4=plot(xx,psychoGauss(p2,xx),'k');
% [aa,ii]=sort(pertSize(~isnan(reactionTime)));
% aux=reactionSign(~isnan(reactionTime));
% MM=100;
% aux=[ones(MM,1); aux(ii);-1*ones(MM,1)]; %Sorted by perturbation size
% aa=[ones(MM,1)*-500 ; aa; ones(MM,1)*500];
% leftCumProb=cumsum(aux==1)./cumsum(abs(aux)); %Leftward choices  over all choices, for more negative perturbations than the current one
% rightCumProb=(sum(aux==-1)-cumsum(aux==-1))./(sum(abs(aux))-cumsum(abs(aux))) ; %Rightward over all, for more positive perturbations
% k1=find(leftCumProb>=rightCumProb,1,'last');
% k2=find(leftCumProb<=rightCumProb,1,'first');
% %pp2=plot(pp,auxYY,'rx');
% pp3=plot(pp,auxY,'b');
% legend([ph(1) pp1 pp3],{'Empiric % of ''<-''','MLE of ''<-''','% of NR'},'Location','SouthEast')
% title('Psychometric fit to initial keypress')
% xlabel('Perturbation speed (mm/s)')
% ylabel('Prob.')
% set(gca,'YTick',[0:.1:1],'XTick',[-350 -250 -150 -75 0 75  150 250 350])
% grid on
% text(-350, .75,['a=' num2str(p(1))])
% text(-350, .65,['b=' num2str(p(2))])
% text(-350, .45,['p=(1+e^{(\Delta V-a)/b})^{-1}'])
% %text(-300, .95,['\mu=' num2str(p2(1))])
% %text(-300, .85,['\sigma=' num2str(p2(2))])
% try
% text(-350, .55,['th=' num2str(p2(3))])
% end
% text(-300,.95,['m=' num2str(.5*(aa(k1)+aa(k2)),2)])
% plot(aa(k1)*[1 1],[0 1],'k--')
% plot(aa(k2)*[1 1],[0 1],'k--')
% 
% subplot(3,4,12) %Psychometric curve of first keypress per STRIDE being right-ward
% hold on
% trialStrides2=zeros(size(trialStrides));
% N=1;
% aux=bsxfun(@plus,inds,[0:N:pDuration-1]);
% trialStrides2(aux(:))=1;
% trialStrides2=trialStrides2 & trialStrides;
% auxY=firstResponse(trialStrides2(:));
% auxX=vD(trialStrides2(:));
% clear ph
% for i=1:length(inds)
%     ph(i)=plot(vD(inds(i)),firstResponse(inds(i)+1),'o','MarkerFaceColor',p1((vD(inds(i))==pp)).Color,'MarkerEdgeColor','none','MarkerSize',5);
% end
% pp0=plot(auxX+5*(randn(size(auxY))), auxY,'.','Color',.7*ones(1,3)); %Indiv trials with some noise to see multiple equal responses
% try
%     if psychoLapse==0
%         [p,~] = fitPsycho(auxX(auxY~=.5),auxY(auxY~=.5),'MLE');
%     else
%         [p,~] = fitPsycho(auxX(auxY~=.5),auxY(auxY~=.5),'MLEsat');
%     end
% 
% catch
%     warning('Could not fit psycho!')
%     p=[0,1e6];
% end
% %[p,~] = fitPsycho(auxX(auxY~=.5),auxY(auxY~=.5),'MLEsat');
% [ff,gg]=psycho(p,auxX);
% xx=[-400:400];
% pp1=plot(xx,psycho(p,xx),'k');
% legend([ph(1) pp0 pp1],{'First stride','All strides','MLE fit'},'Location','SouthEast')
% title('Psychometric fit to first keypress in all strides')
% xlabel('Speed diff (mm/s)')
% ylabel('Prob. of pressing ''<-'' as first key')
% set(gca,'YTick',[0:.1:1],'XTick',[-350 -250 -150 -75 0 75  150 250 350])
% grid on
% text(-300, .75,['a=' num2str(p(1))])
% text(-300, .65,['b=' num2str(p(2))])
% try
% text(-300, .55,['th=' num2str(p(3))])
% end




end

