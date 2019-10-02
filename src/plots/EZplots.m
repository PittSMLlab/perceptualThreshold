function [f1,f2] = EZplots(superSuperT)

superSuperT.isFirstInBlock=[1;diff(superSuperT.blockNo)~=0].*sign(superSuperT.pertSize);
superSuperT=superSuperT(superSuperT.isFirstInBlock==0,:); %Test, remove the first trial in each block
[cmap,unsignedMap]=probeColorMap(23);
%% EZ modeling
clear driftRate
superSuperT.response=1-(superSuperT.initialResponse+1)/2;

%Group model:
[pSize,driftRateA,noiseA,delayA,biasA,t73A,alphaA]=fitEZ_mine(superSuperT);
sA=t73A/(noiseA.^2).^(1/alphaA);
aA=accFactor(pSize,biasA,sA,alphaA,noiseA);
tA=delayA+rtFactor(pSize,biasA,sA,alphaA,noiseA);

%All in a single plot:
f1=figure('Units','pixels','InnerPosition',[100 100 2*300 2*300]);
set(gcf,'Colormap',unsignedMap)
sSize=40;

%Individual model fits:
Nsubs=9;
fitMode={'RT','acc','mixed'};
for j=3%1:length(fitMode) %Different fitting modes
    clear noise aAll tAll driftRate delay bias t73 alpha
    aAll=0;
    tAll=0;
    for i=1:Nsubs
        aux=superSuperT(superSuperT.subID==i,:); %Single subj
        [pSize,driftRate(i,:),noise(i,:),delay(i),bias(i),t73(i),alpha(i)]=fitEZ_mine(aux,fitMode{j});
        a=accFactor(pSize,bias(i),t73(i),alpha(i),noise(i,:));
        aAll=aAll+a;
        t=delay(i)+rtFactor(pSize,bias(i),t73(i),alpha(i),noise(i,:)');
        tAll=tAll+t;
        difficulty(i,:)=driftRate(i,:)./noise(i,:).^2; 
        %Add: estimate a 'mixing' param given the true, vs. the estimated
        %accuracies.
    end
    aAll=aAll/Nsubs;
    tAll=tAll/Nsubs;
    
    subplot(2,2,1)
    hold on
    plot(pSize,aAll,'LineWidth',2,'DisplayName',[fitMode{j} ' fitted']);

    subplot(2,2,2)
    hold on
    plot(pSize,tAll,'LineWidth',2,'DisplayName',[fitMode{j} ' fitted'])
 
    subplot(2,2,3)
    hold on
    plot(aAll,tAll,'LineWidth',2,'DisplayName',[fitMode{j} ' fitted'])
    
    subplot(2,2,4)
    hold on
    plot(pSize,mean(difficulty,1),'LineWidth',2,'DisplayName',[fitMode{j} ' fitted'])

end

%Then plot data:
subplot(2,2,1)
hold on
G=findgroups(abs(superSuperT.pertSize));
superSuperT.cr=double(superSuperT.initialResponse==-sign(superSuperT.pertSize));
superSuperT.cr(isnan(superSuperT.initialResponse))=nan;
acc=splitapply(@(x) nanmean(x),superSuperT.cr,G);
eacc=splitapply(@(x) nanstd(x)/sqrt(sum(~isnan(x))),superSuperT.cr,G);
acc(1)=.5;
ss=scatter(pSize,acc,sSize,pSize,'filled','MarkerEdgeColor','w','DisplayName','group data'); %all data
errorbar(pSize,acc,eacc,'k','LineStyle','none')
xlabel('|probe size| (mm/s)')
ylabel('% correct')
legend('Location','SouthEast','Box','off')
set(gca,'YLim',[.5 1.01])
uistack(ss,'top')

 subplot(2,2,2)
hold on
rt=splitapply(@(x) nanmean(x),superSuperT.reactionTime,G);
ert=splitapply(@(x) nanstd(x)/sqrt(sum(~isnan(x))),superSuperT.reactionTime,G);
ss=scatter(pSize,rt,sSize,pSize,'filled','MarkerEdgeColor','w'); %all data
errorbar(pSize,rt,ert,'k','LineStyle','none')
xlabel('|probe size| (mm/s)')
ylabel('mean RT (s)')
set(gca,'YLim',[0 7])
uistack(ss,'top')

subplot(2,2,3)
hold on
ss=scatter(acc,rt,sSize,pSize,'filled','MarkerEdgeColor','w','DisplayName','group data'); %all data
errorbar(acc,rt,ert,'k','LineStyle','none')
errorbar(acc,rt,eacc,'Horizontal','k','LineStyle','none')
xlabel('% correct')
ylabel('mean RT (s)')
uistack(ss,'top')
set(gca,'XLim',[.5 1.01],'YLim',[0 7])

subplot(2,2,4)
hold on
xlabel('|probe size| (mm/s)')
ylabel('difficulty^{-1}')


%% Alt EZ modeling:
Nsubs=9;
f2=figure('Name','Analysis of choice-fitted indiv. models'); hold on; 
for i=1:Nsubs
    aux=superSuperT(superSuperT.subID==i,:); %Single subj
    aux=aux(aux.pertSize~=250,:); %Eliminating 250mm/s pert
    aux.prevFinalSpeed=[0;aux.lastSpeedDiff(1:end-1)].*[0;diff(aux.blockNo)==0];
    aux.correctResponse=aux.initialResponse==-sign(aux.pertSize);
    aux.absPertSize=abs(aux.pertSize);
    aux.absDistToLastEnd=abs(aux.pertSize-aux.prevFinalSpeed);
    aux.leftResponse=aux.initialResponse==-1+.5*isnan(aux.initialResponse);
    aux.prevSize=[0;aux.pertSize(1:end-1)].*[0;diff(aux.blockNo)==0];
       %Get data:
    for j=1%:2
        switch j
            case 1
            X=aux;
            responseVar='leftResponse';
            mdl='pertSize+prevSize';
            mdl='pertSize';
            case 2
            X=aux(~isnan(aux.initialResponse) & aux.pertSize~=0,:); %No null trials, no non-response trials
            responseVar='correctResponse';
            mdl='absPertSize-1';
        end
        mm=fitglm(X,[responseVar '~' mdl],'Distribution','binomial','Link','logit'); %Logisitc regression, excluding null responses
        ci=mm.coefCI;
        
        subplot(2,5,(j-1)*5+1) % probe size vs. choice
        hold on
        %data
        if j==1
            xvar='pertSize';
        else
            xvar='absPertSize';
        end
        set(gca,'ColorOrderIndex',i)
        G=findgroups(aux.(xvar));
        x=splitapply(@(z) mean(z), aux.(xvar),G);
        y=splitapply(@(z) mean(z),aux.(responseVar),G);
        %x=x(1:end-1)+x(2:end);
        %x=.5*x(1:2:end);
        %y=y(1:end-1)+y(2:end);
        %y=.5*y(1:2:end);
        %scatter(x,y,'filled');
        set(gca,'ColorOrderIndex',i)
        mm.plotPartialDependence(xvar)
                xlabel('probe size (mm/s)')
        ylabel('% left choice')
        title('choice vs. probe size')
        
        for k=1:size(ci,1)
            subplot(2,5,(j-1)*5+k+1)
            hold on
            c=mm.Coefficients.Estimate(k);
            bar(i,c)
            errorbar(i,c,c-ci(k,1),ci(k,2)-c,'k')
            if k==1
                title('bias')
            else 
                title('slope')
            end
            xlabel('subject ID')
            set(gca,'XTick',1:9)
            ylabel('parameter value')
        end

        
        subplot(2,5,(j-1)*5+5) %RT vs choice
        hold on
         set(gca,'ColorOrderIndex',i)
        y=splitapply(@(z) nanmean(z),aux.reactionTime,G);
        xx=splitapply(@(x) nanmean(x),mm.predict(aux),G); %Predicted responses
        [xx,idx]=sort(xx);
        %scatter(y(idx),xx,'filled') % (true) RT vs. expected accuracy for that probe size
         set(gca,'ColorOrderIndex',i)
         auxT=(2*xx-1)./log(xx./(1-xx));
         auxT(xx==.5)=.5; %Forcing limit
         pp=polyfit(auxT(~isnan(auxT)),y(idx(~isnan(auxT))),1);
         rt=pp(1)*auxT+pp(2);
        plot(rt,xx) %EZ DDM prediction
        xlabel('expected RT (s)')
        ylabel('% left choice')
        
        subplot(2,5,(j-1)*5+4) %Probe size vs. RT
                hold on
         set(gca,'ColorOrderIndex',i)
         xx=splitapply(@(x) nanmean(x),aux.(xvar),G); %Probe size
         [xx,idx]=sort(xx);
         y=y(idx);
         %scatter(xx,y,'filled')
         set(gca,'ColorOrderIndex',i)
        plot(xx,rt) %EZ DDM prediction
        ylabel('expected RT (s)')
        xlabel('probe size (mm/s)')
    end
end

