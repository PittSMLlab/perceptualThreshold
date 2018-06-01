function fh=ssPlots(trialData)
%%
B=findgroups(trialData.pertSize); %pertSize>0 means vR>vL
pp=unique(trialData.pertSize);
cmap=parula(length(pp)); %parula, jet, redbluecmap
cmap=cmap*.8;
 
fh=figure('Units','Normalized','OuterPosition',[0 0 1 1]);
Q=6;

%%
nullTrials=trialData.pertSize==0;
correctResponses=trialData.initialResponse==-sign(trialData.pertSize) & ~nullTrials; %Negative response means LEFT IS SLOW (RIGHT IS FAST) choice
nonResponse=isnan(trialData.initialResponse) & ~nullTrials;
incorrectResponses=trialData.initialResponse==sign(trialData.pertSize) & ~nullTrials;
B2=findgroups(abs(trialData.pertSize)); %pertSize>0 means vR>vL
pp2=unique(abs(trialData.pertSize));
colorOff=3;

%% Fourth plot: steady-state as function of perturbation size
subplot(2,Q,1:2)
%scatter(trialData.pertSize,trialData.lastSpeedDiff,10,zeros(1,3))
%hold on
fun=@nanmedian;
S4=splitapply(fun,trialData.lastSpeedDiff,B);
s1=scatter(trialData.pertSize,trialData.lastSpeedDiff,5,.5*ones(1,3),'filled');
hold on
s2=scatter(pp,S4,70,pp,'filled');
grid on
ylabel('Final speed (mm/s)') 
xlabel('vL>vR      PERTURBATION         vL<vR')
axis([-360 360 -150 150])
hold on

y=monoLS(S4);
%plot(pp,y,'k')
y=splitapply(@nanmean,trialData.lastSpeedDiff,B);
p1=plot(pp,y,'k');
y2=splitapply(@(x) nanstd(x)/sqrt(numel(x)),trialData.lastSpeedDiff,B);
ptc=patch([pp' fliplr(pp')],[y'+y2' fliplr(y'-y2')],.5*ones(1,3),'FaceAlpha',.3,'EdgeColor','none');
uistack(ptc,'bottom')
legend([s1 s2 p1],{'Indiv. trials','Median','Mean \pm ste'},'Location','NorthWest')

% subplot(2,Q,[1:2]+Q)
% S=splitapply(@nanmedian,trialData.lastSpeedDiff,B); 
% S2=splitapply(@nanmedian,trialData.lastSpeedDiff .* sign(trialData.pertSize),B2); 
% scatter(pp2,S2,80,.4*ones(1,3),'filled')
% hold on
% grid on
% ylabel('Final speed (mm/s)') 
% xlabel('ABS SPEED PERTURBATION (mm/s)')
% scatter(pp(pp>0),S(pp>0),20,cmap(end,:),'filled')
% scatter(abs(pp(pp<0)),-S(pp<0),20,cmap(1,:),'filled')

%% Sixth plot: clicking rates
subplot(2,Q,3:4)
fun=@mean;
fun=@median;
S=splitapply(fun, trialData.Lclicks+trialData.Rclicks,B);
s1=scatter(pp,S,70,.4*ones(1,3),'filled');
hold on
S2=splitapply(fun, trialData.Lclicks-trialData.Rclicks,B);
s2=scatter(pp,S2.*sign(pp+1e-5),70,pp,'filled'); %PLotting as positive if net was in the right direction, negative otherwise
S3=splitapply(fun, (trialData.Rclicks),B);
s3=scatter(pp,S3,20,cmap(1,:),'filled');
S4=splitapply(fun, (trialData.Lclicks),B);
s4=scatter(pp,S4,20,cmap(end,:),'filled');
grid on
xlabel('PERTURBATION (mm/s)')
ylabel('(median) Click rate (per trial)')
legend([s1 s2 s3 s4],{'Total','Net', 'R','L'},'Location','North')
axis([-380 380 0 50])
plot([-350,0,350],[350,0,350]/7,'k','DisplayName','Rate required for full correction')
%To Do: add mean number of clicks required for full correction


%TODO: plot correct presses per trial vs. incorrect presses per trial, as
%function of abs(pertSize)

%%
trialData.prevSize=[0;trialData.pertSize(1:end-1)]; 
trialData.prevSize(trialData.pertSize==-250*((-1).^trialData.blockNo))=0; %Assigning NaN to previous perturbation for first trial in each block (+250 in even blocks, -250 in odd)

%Creating binary response variable(s):
trialData.leftResponse=trialData.initialResponse==-1;
trialData.rightResponse=trialData.initialResponse==1;
trialData.noResponse=isnan(trialData.initialResponse);

%
trialData.invAbsPertSize=1./abs(trialData.pertSize);
trialData.absPertSize=abs(trialData.pertSize);
trialData.pertSign=sign(trialData.pertSize);
%
trialData.correctResponse=trialData.initialResponse==-sign(trialData.pertSize) & ~nullTrials; 
%
aux=trialData.subID;
aux(aux==1)=10;
aux(aux==5)=1;
trialData.ID=categorical(aux);
%
trialData.netClicks=trialData.Lclicks-trialData.Rclicks;
trialData.absNetClicks=trialData.netClicks .* sign(trialData.pertSize+1e-4);
trialData.logNetClicks=sign(trialData.netClicks).*log(abs(trialData.netClicks));
trialData.logPertSize=sign(trialData.pertSize).*log(abs(trialData.pertSize)+1e-9);
%
trialData.incorrectResponse=~trialData.correctResponse;
trialData.missingClicks=trialData.pertSize/7 -trialData.netClicks;
%% Modeling net click as a function of pertSize, etc
X=trialData(:,:);
mm=fitlm(X,'netClicks~pertSize*pertSign+blockNo:pertSign+blockNo:pertSize+correctResponse:pertSize')%,'Distribution','poisson','Link','identity','DispersionFlag',true); %Logisitc regression, excluding null responses
text(380,40,removeTags(evalc('mm.disp')),'FontSize',9,'Clipping','off')
mm1=fitlm(X,'netClicks~incorrectResponse:pertSize+pertSize');%,'Distribution','normal','Link','logit','DispersionFlag',true); %Logisitc regression, excluding null responses
text(380,0,removeTags(evalc('mm1.disp')),'FontSize',9,'Clipping','off')
hold on
plot([-350 0 350],[350 0 350]*mm1.Coefficients.Estimate(2),'r','DisplayName','Best LM fit')
%mm2=fitlm(X,'missingClicks~incorrectResponse:pertSize+pertSize');%,'Distribution','normal','Link','logit','DispersionFlag',true); %Logisitc regression, excluding null responses
%text(380,-15,removeTags(evalc('mm2.disp')),'FontSize',9,'Clipping','off')
%mm2=fitlm(X,'logNetClicks~correctResponse:pertSize+pertSize');%,'Distribution','normal','Link','identity','DispersionFlag',true); %Logisitc regression, excluding null responses
%text(380,-15,removeTags(evalc('mm2.disp')),'FontSize',9,'Clipping','off')
mm3=fitlm(X,'netClicks~pertSize+correctResponse:pertSize+ID*pertSize-1')%,'Distribution','poisson','Link','identity','DispersionFlag',true); %Logisitc regression, excluding null responses
text(400,-55,removeTags(evalc('mm3.disp')),'FontSize',9,'Clipping','off')
%mm2=fitlm(X,'absNetClicks~absPertSize*pertSign+blockNo+correctResponse*absPertSize')%,'Distribution','poisson','Link','identity','DispersionFlag',true); %Logisitc regression, excluding null responses
%text(-350,-60,removeTags(evalc('mm2.disp')),'FontSize',9,'Clipping','off')

%% Modeling final steady-state
%X=trialData(trialData.correctResponse,:);
mm=fitlm(X,'lastSpeedDiff~pertSize+incorrectResponse:pertSize-1')%,'Distribution','poisson','Link','identity','DispersionFlag',true); %Logisitc regression, excluding null responses
text(-400,-30,removeTags(evalc('mm.disp')),'FontSize',9,'Clipping','off')

%X=trialData(trialData.correctResponse,:);
mm=fitlm(X,'lastSpeedDiff~pertSize+incorrectResponse:pertSize+pertSize*ID')%,'Distribution','poisson','Link','identity','DispersionFlag',true); %Logisitc regression, excluding null responses
text(-1400,-45,regexprep(removeTags(evalc('mm.disp')),'model:\n','model: \n'),'FontSize',9,'Clipping','off')
%mm1=fitlm(X,'lastSpeedDiff~pertSize');%,'Distribution','normal','Link','logit','DispersionFlag',true); %Logisitc regression, excluding null responses
%text(-1180,-60,removeTags(evalc('mm1.disp')),'FontSize',9,'Clipping','off')
end

