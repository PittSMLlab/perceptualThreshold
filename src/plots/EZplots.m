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

%Individual models:
Nsubs=9;
aAll=0;
tAll=0;
aAll2=0;
tAll2=0;
aAll3=0;
tAll3=0;
for i=1:Nsubs
    aux=superSuperT(superSuperT.subID==i,:); %Single subj
    [pSize,driftRate(i,:),noise(i),delay(i),bias(i),t73(i),alpha(i)]=fitEZ_mine(aux,'RT');
    a=accFactor(pSize,bias(i),t73(i),alpha(i),noise(i));
    aAll=aAll+a;
    t=delay(i)+rtFactor(pSize,bias(i),t73(i),alpha(i),noise(i));
    tAll=tAll+t;
    difficulty(i,:)=driftRate(i,:)/noise(i)^2; 
    %Add: estimate a 'mixing' param given the true, vs. the estimated
    %accuracies.
    
    [pSize,driftRate2(i,:),noise2(i),delay2(i),bias2(i),t732(i),alpha2(i)]=fitEZ_mine(aux,'acc');
    a2=accFactor(pSize,bias2(i),t732(i),alpha2(i),noise2(i));
    aAll2=aAll2+a2;
    t2=delay2(i)+rtFactor(pSize,bias2(i),t732(i),alpha2(i),noise2(i));
    tAll2=tAll2+t2;
    difficulty2(i,:)=driftRate2(i,:)/noise2(i)^2; 
    
    [pSize,driftRate3(i,:),noise3(i),delay3(i),bias3(i),t733(i),alpha3(i),mix3(i)]=fitEZ_mine(aux,'RTalt');
    a3=accFactor(pSize,bias(i),t73(i),alpha(i),noise(i));
    a3=(1-mix3(i))*a + mix3(i)*(1-a);
    aAll3=aAll3+a3;
    t3=delay(i)+rtFactor(pSize,bias3(i),t733(i),alpha3(i),noise3(i));
    tAll3=tAll3+t3;
    difficulty3(i,:)=driftRate3(i,:)/noise3(i)^2; 
end
aAll=aAll/Nsubs;
tAll=tAll/Nsubs;
aAll2=aAll2/Nsubs;
tAll2=tAll2/Nsubs;
aAll3=aAll3/Nsubs;
tAll3=tAll3/Nsubs;

%All in a single plot:
f1=figure('Units','pixels','InnerPosition',[100 100 2*300 2*300]);
set(gcf,'Colormap',unsignedMap)
sSize=40;
subplot(2,2,1)
hold on
%set(gca,'Colormap',cmap)
G=findgroups(abs(superSuperT.pertSize));
superSuperT.cr=double(superSuperT.initialResponse==-sign(superSuperT.pertSize));
superSuperT.cr(isnan(superSuperT.initialResponse))=nan;
acc=splitapply(@(x) nanmean(x),superSuperT.cr,G);
eacc=splitapply(@(x) nanstd(x)/sqrt(sum(~isnan(x))),superSuperT.cr,G);
acc(1)=.5;
ss=scatter(pSize,acc,sSize,pSize,'filled','MarkerEdgeColor','w','DisplayName','Group data'); %all data
p1=plot(pSize,aAll2,'k','LineWidth',2,'DisplayName','Acc. fitted');
p2=plot(pSize,aAll,'LineWidth',2,'DisplayName','RT fitted');
%p3=plot(pSize,aAll3,'LineWidth',2,'DisplayName','RTalt fitted');
errorbar(pSize,acc,eacc,'k','LineStyle','none')
xlabel('|Probe size| (mm/s)')
ylabel('% correct')
legend([ss,p1,p2],'Location','SouthEast','Box','off')
set(gca,'YLim',[.5 1.01])
uistack(ss,'top')
subplot(2,2,2)
hold on
rt=splitapply(@(x) nanmean(x),superSuperT.reactionTime,G);
ert=splitapply(@(x) nanstd(x)/sqrt(sum(~isnan(x))),superSuperT.reactionTime,G);
ss=scatter(pSize,rt,sSize,pSize,'filled','MarkerEdgeColor','w'); %all data
plot(pSize,tAll2,'k','LineWidth',2,'DisplayName','Acc. fitted')
plot(pSize,tAll,'LineWidth',2,'DisplayName','RT fitted')
%plot(pSize,tAll3,'LineWidth',2,'DisplayName','RTalt fitted');
errorbar(pSize,rt,ert,'k','LineStyle','none')
xlabel('|Probe size| (mm/s)')
ylabel('Mean RT (s)')
set(gca,'YLim',[0 7])
uistack(ss,'top')
subplot(2,2,3)
hold on
ss=scatter(acc,rt,sSize,pSize,'filled','MarkerEdgeColor','w','DisplayName','Group data'); %all data
plot(aAll2,tAll2,'k','LineWidth',2,'DisplayName','Acc. fitted')
plot(aAll,tAll,'LineWidth',2,'DisplayName','RT fitted')
%plot(aAll3,tAll3,'LineWidth',2,'DisplayName','RTalt fitted');
errorbar(acc,rt,ert,'k','LineStyle','none')
errorbar(acc,rt,eacc,'Horizontal','k','LineStyle','none')
xlabel('% correct')
ylabel('Mean RT (s)')
%axis tight
uistack(ss,'top')
set(gca,'XLim',[.5 1.01],'YLim',[0 7])
subplot(2,2,4)
hold on
plot(pSize,mean(difficulty2,1),'k','LineWidth',2,'DisplayName','Acc fitted')
plot(pSize,mean(difficulty,1),'LineWidth',2,'DisplayName','Acc fitted')
xlabel('|Probe size| (mm/s)')
ylabel('1/difficulty')
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
                xlabel('Probe size (mm/s)')
        ylabel('% left choice')
        title('Choice vs. probe size')
        
        for k=1:size(ci,1)
            subplot(2,5,(j-1)*5+k+1)
            hold on
            c=mm.Coefficients.Estimate(k);
            bar(i,c)
            errorbar(i,c,c-ci(k,1),ci(k,2)-c,'k')
            if k==1
                title('Bias')
            else 
                title('Slope')
            end
            xlabel('Subject ID')
            set(gca,'XTick',1:9)
            ylabel('Parameter value')
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
        xlabel('Expected RT (s)')
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
        ylabel('Expected RT (s)')
        xlabel('Probe size (mm/s)')
    end
end
end

