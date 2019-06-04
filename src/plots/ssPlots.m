function fh=ssPlots(trialData)
%%
B=findgroups(trialData.pertSize); %pertSize>0 means vR>vL
pp1=unique(trialData.pertSize);

[cmap,unsigned]=probeColorMap(23);


%%
trialData.noResponse=isnan(trialData.initialResponse);
trialData.correctResponse= trialData.initialResponse==-sign(trialData.pertSize);
trialData.incorrectResponse= trialData.initialResponse==sign(trialData.pertSize);
trialData.prevSize=[0;trialData.pertSize(1:end-1)]; 
trialData.isFirstInBlock=[1;diff(trialData.blockNo)~=0].*sign(trialData.pertSize);
trialData.prevSize(trialData.isFirstInBlock~=0)=0; %Assigning NaN to previous perturbation for first trial in each block
trialData.prevFinalSpeed=[0;trialData.lastSpeedDiff(1:end-1)].*[0;diff(trialData.subID)==0].*[0;diff(trialData.blockNo)==0];
trialData.pertSign=sign(trialData.pertSize);
trialData.blockNo=trialData.blockNo-1; %To make the first block have 0 contribution in linear model
%% %Remove first trial in each block
trialData=trialData(trialData.isFirstInBlock==0,:); %Test, remove the first trial in each block
%% steady-state as function of perturbation size and correct/incorrect
fh=figure('Units','pixels','InnerPosition',[100 100 3*300 1*300]);
for i=1:2
    switch i
        case 1
            X=trialData(trialData.correctResponse==1,:);
            ttl='Correct trials';
        case 2
            X=trialData(trialData.incorrectResponse==1  | trialData.pertSize==0,:);
            ttl='Incorrect trials';
    end
    
    subplot(1,2,i)
    fun=@nanmean;
    Ba=findgroups(X.pertSize);
    Ya=X.lastSpeedDiff; %All responses
    %Ya(~trialData.correctResponse)=nan;
    S4=splitapply(fun,Ya,Ba);
    M4=splitapply(@(x) nanmedian(x), Ya,Ba);
    pp=splitapply(fun,X.pertSize,Ba);
    %E4=splitapply(@(x) nanstd(x)/sqrt(sum(~isnan(x))),Ya,Ba);
    E4=splitapply(@(x) nanstd(x),Ya,Ba);
    %s1=scatter(trialData.pertSize,trialData.lastSpeedDiff,5,.5*ones(1,3),'filled');
    hold on
    s2=scatter(pp,M4,50,pp,'filled','MarkerEdgeColor','w');
    colormap(cmap)
    ptc=errorbar(pp,S4,E4,'k','LineStyle','none','LineWidth',1);
    E4=splitapply(@(x) prctile(x,75),Ya,Ba)-M4;
    E5=M4-splitapply(@(x) prctile(x,25),Ya,Ba);
    %ptc2=errorbar(pp,M4,E4,E5,'k','LineStyle','none');
    uistack(ptc,'bottom')
    %    uistack(ptc2,'bottom')
    grid on
    ylabel('Reported PSE (mm/s)') 
    xlabel('probe size (mm/s)')
    axis([-360 360 -250*1 250*1])
    hold on

    y=splitapply(@nanmean,Ya,Ba);
    %y2=splitapply(@(x) nanstd(x)/sqrt(numel(x)),Ya,Ba); %ste
    %ptc=patch([pp' fliplr(pp')],[y'+y2' fliplr(y'-y2')],.5*ones(1,3),'FaceAlpha',.3,'EdgeColor','none');
    %uistack(ptc,'bottom')
    if i==1
    legend([s2 ptc],{'median','mean \pm std','25-75 percentile'},'Location','SouthEast','Box','off')
    clear pval
    %For correct responses only: test for differences from 0:
    for j=1:length(pp1)
        resp=X.lastSpeedDiff(X.pertSize==pp1(j));
        [~,pval(j)]=ttest(resp);
    end
    [h,pTh]=BenjaminiHochberg(pval,.05,true);
    disp(['P threshold =' num2str(pTh) ', significantly different from 0:'])
    disp(num2str(pp1(h==1)'))
    end
    title(ttl)
end

%% Analysis of CORRECT TRIALS ONLY
disp('---------------Correct trials only!---------------------')
%Test for foreign effects:
X=trialData(trialData.correctResponse==1,:); %Only correct responses for analysis
%Slope model:
frml='lastSpeedDiff~pertSize+prevSize+pertSize:blockNo+pertSize:pertSign-1';
mm=fitlm(X,frml)
mm=mm.step('Upper',frml,'Criterion','SSE','PEnter',0,'PRemove',0.05,'Nsteps',Inf)

%Dead-zone model:
bestRMSE=Inf;
for th=0:5:300 %Line search over possible thresholds
    X.expl= sign(X.pertSize).*max(abs(X.pertSize)-th,0);
    mm0=fitlm(X,'lastSpeedDiff~expl-1');
    RMSE=mm0.RMSE;
    if RMSE<bestRMSE
        bestRMSE=RMSE;
        bestTh=th;
    end
end
disp(['Deadzone model RMSE: ' num2str(bestRMSE)])
disp(['Deadzone threshold: ' num2str(bestTh)])
disp(['Deadzone slope: ' num2str(mm0.Coefficients.Estimate)])
Nsamp=mm.DFE+1; %Cancels out
Ndof=mm.DFE;
F=(Nsamp*(RMSE-bestRMSE).^2)/1 / ((Nsamp* bestRMSE.^2)/(Ndof-1));
disp(['F-stat: ' num2str(F) ', p=' num2str(1-fcdf(F,1,Ndof-1,'upper'))]) %Deadzone model is nested

%% Analysis of all trials:
disp('----------------------All trials-----------------------')
X=trialData; %Both correct and incorrect trials!
% Test that distribution is bimodal:
mm=fitlm(X,'lastSpeedDiff~pertSize+correctResponse:pertSize-1') %Allowing for different correct and incorrect responses
%mm.anova

%Compare two models for correct responses: linear slope (e.g. subjects correct x%), and fixed
%threshold (e.g. subjects correct until they fall within a certain band)
mm=fitlm(X,'lastSpeedDiff~pertSize-1');
disp(['Linear model RMSE:' num2str(mm.RMSE)])
disp(['Linear model slope:' num2str(mm.Coefficients.Estimate)])
L1=nansum(abs(mm.Residuals.Raw));
%Threshold model:
bestRMSE=Inf;
for th=0:5:200 %Line search over possible thresholds
    X.expl= sign(X.pertSize).*min(abs(X.pertSize),th);
    %RMSE=sqrt(nanmean((X.lastSpeedDiff-X.expl).^2)); %Fixing slope to 1: this does not result in RMSE better than the linear model
    mm0=fitlm(X,'lastSpeedDiff~expl-1'); %Slopey threshold: allows for some correction even of below-threshold probes
    RMSE=mm0.RMSE;
    if RMSE<bestRMSE
        bestRMSE=RMSE;
        bestTh=th;
    end
end
disp(['Threshold model RMSE:' num2str(bestRMSE)])
disp(['Threshold value:' num2str(bestTh)])
disp(['Slope: ' num2str(mm0.Coefficients.Estimate)])
Nsamp=mm.DFE+1; %Cancels out
Ndof=mm.DFE;
F=(Nsamp*(RMSE-bestRMSE).^2)/1 / ((Nsamp* bestRMSE.^2)/(Ndof-1));
disp(['F-stat: ' num2str(F) ', p=' num2str(1-fcdf(F,1,Ndof-1,'upper'))]) %Deadzone model is nested

%% Incorrect trials only:
disp('---------------INCorrect trials only!---------------------')
%Test for foreign effects:
X=trialData(trialData.correctResponse~=1,:); %Only correct responses for analysis
%Slope model:
frml='correction~pertSize-1';
X.correction=X.pertSize-X.lastSpeedDiff;
mm=fitlm(X,frml) %Offset added so stats will reveal difference to no-correction line

nc=X.correction(X.pertSize==0);
m=mean(nc);
s=std(nc);
disp(['Null trial reported PSE: ' num2str(m) ' \pm ' num2str(s) ' (mean \pm std)'])
%figure; hold on;
%mm.plotPartialDependence('pertSize')
%scatter(X.pertSize,X.correction,20,'k','filled')
%%
%figure; %Plotting all data by probe size
%hold on
%pp=unique(trialData.pertSize);
%for i=1:length(pp)
%y=trialData.lastSpeedDiff(trialData.pertSize==pp(i) & trialData.correctResponse==1);
%scatter(pp(i)*ones(size(y)),y,10,'k','filled')
%
%y=trialData.lastSpeedDiff(trialData.pertSize==pp(i) & trialData.correctResponse~=1);
%scatter(pp(i)*ones(size(y)),y,10,.5*ones(1,3),'filled')
%end
%% Figure: clicking rates
% figure
% subplot(2,Q,3:4)
% fun=@mean;
% fun=@median;
% S=splitapply(fun, trialData.Lclicks+trialData.Rclicks,B);
% s1=scatter(pp,S,70,.4*ones(1,3),'filled');
% hold on
% S2=splitapply(fun, trialData.Lclicks-trialData.Rclicks,B);
% s2=scatter(pp,S2.*sign(pp+1e-5),70,pp,'filled'); %PLotting as positive if net was in the right direction, negative otherwise
% S3=splitapply(fun, (trialData.Rclicks),B);
% s3=scatter(pp,S3,20,cmap(1,:),'filled');
% S4=splitapply(fun, (trialData.Lclicks),B);
% s4=scatter(pp,S4,20,cmap(end,:),'filled');
% grid on
% xlabel('PERTURBATION (mm/s)')
% ylabel('(median) Click rate (per trial)')
% legend([s1 s2 s3 s4],{'Total=R+L','Net=|L-R|', 'R','L'},'Location','North')
% axis([-380 380 0 50])
% plot([-350,0,350],[350,0,350]/7,'k','DisplayName','Rate required for full correction')
%To Do: add mean number of clicks required for full correction


%TODO: plot correct presses per trial vs. incorrect presses per trial, as
%function of abs(pertSize)

%%
% trialData.prevSize=[0;trialData.pertSize(1:end-1)]; 
% trialData.prevSize(trialData.pertSize==-250*((-1).^trialData.blockNo))=0; %Assigning NaN to previous perturbation for first trial in each block (+250 in even blocks, -250 in odd)
% 
% %Creating binary response variable(s):
% trialData.leftResponse=trialData.initialResponse==-1;
% trialData.rightResponse=trialData.initialResponse==1;
% trialData.noResponse=isnan(trialData.initialResponse);
% 
% %
% trialData.invAbsPertSize=1./abs(trialData.pertSize);
% trialData.absPertSize=abs(trialData.pertSize);
% trialData.pertSign=sign(trialData.pertSize);
% %
% trialData.correctResponse=trialData.initialResponse==-sign(trialData.pertSize) & ~nullTrials; 
% %
% aux=trialData.subID;
% aux(aux==1)=10;
% aux(aux==5)=1;
% trialData.ID=categorical(aux);
% %
% trialData.netClicks=trialData.Lclicks-trialData.Rclicks;
% trialData.earlyNetClicks=trialData.eLclicks-trialData.eRclicks;
% %trialData.absNetClicks=trialData.netClicks .* sign(trialData.pertSize+1e-4);
% %trialData.logNetClicks=sign(trialData.netClicks).*log(abs(trialData.netClicks));
% %trialData.sqrtPertSize=sign(trialData.pertSize).*sqrt(abs(trialData.pertSize)+1e-9);
% %
% trialData.incorrectResponse=~trialData.correctResponse;
% trialData.missingClicks=trialData.pertSize/7 -trialData.netClicks;
% trialData.totalClicks=trialData.Lclicks+trialData.Rclicks;
% trialData.badTrial=trialData.totalClicks>max(10,2.5*abs(trialData.netClicks)); %If on a trial we observe more than 15 clicks (total), and there is less than a 75/25 split between left and right, we label them bad. %We could do this based on early net clicks being of a different sign, but on the same order, as the last clicks
% trialData.badTrial= trialData.totalClicks> 3*abs(trialData.netClicks) & abs(trialData.earlyNetClicks)>5 ; %This marks 2% of trials as bad
% trialData.correctedAmount=trialData.pertSize-trialData.lastSpeedDiff;
%% Modeling net click as a function of pertSize, etc
%X=trialData(~trialData.badTrial,:);
%X2=trialData(~trialData.incorrectResponse,:);
%mm=fitlm(X,'netClicks~pertSize*pertSign+blockNo:pertSign+blockNo:pertSize+correctResponse:pertSize');%,'Distribution','poisson','Link','identity','DispersionFlag',true); %Logisitc regression, excluding null responses
%text(420,40,removeTags(evalc('mm.disp')),'FontSize',9,'Clipping','off')
%mm1=fitlm(X,'netClicks~incorrectResponse:pertSize+pertSize-1');%,'Distribution','normal','Link','logit','DispersionFlag',true); %Logisitc regression, excluding null responses
%text(420,7,removeTags(evalc('mm1.disp')),'FontSize',9,'Clipping','off')
%hold on
%plot([-350 0 350],[350 0 350]*mm1.Coefficients.Estimate(2),'k','DisplayName','Best LM fit','LineWidth',2);
%mm=fitlm(X,'netClicks~earlyNetClicks+pertSize-1')%,'Distribution','poisson','Link','identity','DispersionFlag',true); %Logisitc regression, excluding null responses
%text(420,-17,removeTags(evalc('mm.disp')),'FontSize',9,'Clipping','off')

%mm=fitlm(X,'netClicks~earlyNetClicks-1');%,'Distribution','poisson','Link','identity','DispersionFlag',true); %Logisitc regression, excluding null responses
%text(420,-40,removeTags(evalc('mm.disp')),'FontSize',9,'Clipping','off')

%mm=fitlm(X,'correctedAmount~earlyNetClicks-1')%,'Distribution','poisson','Link','identity','DispersionFlag',true); %Logisitc regression, excluding null responses
%text(420,-60,removeTags(evalc('mm.disp')),'FontSize',9,'Clipping','off')
%mm2=fitlm(X,'missingClicks~incorrectResponse:pertSize+pertSize');%,'Distribution','normal','Link','logit','DispersionFlag',true); %Logisitc regression, excluding null responses
%text(380,-15,removeTags(evalc('mm2.disp')),'FontSize',9,'Clipping','off')
%mm2=fitlm(X,'logNetClicks~correctResponse:pertSize+pertSize');%,'Distribution','normal','Link','identity','DispersionFlag',true); %Logisitc regression, excluding null responses
%text(380,-15,removeTags(evalc('mm2.disp')),'FontSize',9,'Clipping','off')

%mm3=fitlm(X,'netClicks~pertSize+correctResponse:pertSize+ID*pertSize-1');%,'Distribution','poisson','Link','identity','DispersionFlag',true); %Logisitc regression, excluding null responses
%text(-450,-40,removeTags(evalc('mm3.disp')),'FontSize',9,'Clipping','off')

%mm2=fitlm(X,'absNetClicks~absPertSize*pertSign+blockNo+correctResponse*absPertSize')%,'Distribution','poisson','Link','identity','DispersionFlag',true); %Logisitc regression, excluding null responses
%text(-350,-60,removeTags(evalc('mm2.disp')),'FontSize',9,'Clipping','off')

end

