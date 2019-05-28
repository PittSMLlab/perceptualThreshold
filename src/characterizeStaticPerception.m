%% Generate accuracy plots for each subject
%%
clear all
addpath(genpath('../'))
run loadAllDataIntoTable.m
addpath(genpath('../../monoLS'))
%%
[f1,f2]=accPlots(superSuperT);
saveFig(f1,'../fig/allstatic/',['accuracyAll'],0)
saveFig(f2,'../fig/allstatic/',['accuracySubjectAndBlockEffects'],0)
%%
%All subjects:
[fh,f2]=rtPlots(superSuperT);
%saveFig(fh,'../fig/allstatic/',['rtAll'],0)
saveFig(f2,'../fig/allstatic/',['EZdd'],0)
%w/o subject 2:
%fh=rtPlots(superSuperT(superSuperT.subID~=2,:));
%saveFig(fh,'../../fig/all/',['rtAll'],0)

%%
fh=ssPlots(superSuperT);
%saveFig(fh,'../../fig/all/',['ssAll'],0)

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
a2=0;
t2=0;
for i=1:Nsubs
    aux=superSuperT(superSuperT.subID==i,:); %Single subj
    [pSize,driftRate(i,:),noise(i),delay(i),bias(i),t73(i),alpha(i)]=fitEZ_mine(aux);
    s=t73(i)/(noise(i).^2).^(1/alpha(i));
    a=accFactor(pSize,bias(i),s,alpha(i),noise(i));
    a2=a2+a;
    t=delay(i)+rtFactor(pSize,bias(i),s,alpha(i),noise(i));
    t2=t2+t;
    %Add: estimate a 'mixing' param given the true, vs. the estimated
    %accuracies.
end
a2=a2/Nsubs;
t2=t2/Nsubs;
%All in a single plot:
figure
subplot(2,2,1)
hold on
G=findgroups(abs(superSuperT.pertSize));
superSuperT.cr=double(superSuperT.initialResponse==-sign(superSuperT.pertSize));
superSuperT.cr(isnan(superSuperT.initialResponse))=nan;
acc=splitapply(@(x) nanmean(x),superSuperT.cr,G);
acc(1)=.5;
scatter(pSize,acc,50,.4*ones(1,3),'filled') %all data
plot(pSize,a,'k')
plot(pSize,a2,'r')
subplot(2,2,2)
hold on
rt=splitapply(@(x) nanmean(x),superSuperT.reactionTime,G);
scatter(pSize,rt,50,.4*ones(1,3),'filled') %all data
plot(pSize,t,'k')
plot(pSize,t2,'r')
subplot(2,2,3)
hold on
scatter(acc,rt,50,.4*ones(1,3),'filled') %all data
plot(a,t,'k')
plot(a2,t2,'r')

%% Alt EZ modeling:
Nsubs=9;
figure; hold on; 
for i=1:Nsubs
    i
    aux=superSuperT(superSuperT.subID==i,:); %Single subj
    aux=aux(aux.pertSize~=250,:); %Eliminating 250mm/s pert
    aux.prevFinalSpeed=[0;aux.lastSpeedDiff(1:end-1)].*[0;diff(aux.blockNo)==0];
    aux.correctResponse=aux.initialResponse==-sign(aux.pertSize);
    aux.absPertSize=abs(aux.pertSize);
    aux.absDistToLastEnd=abs(aux.pertSize-aux.prevFinalSpeed);
    aux.leftResponse=aux.initialResponse==-1+.5*isnan(aux.initialResponse);
    aux.prevSize=[0;aux.pertSize(1:end-1)].*[0;diff(aux.blockNo)==0];
       %Get data:
    for j=1:2
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
        mm=fitglm(X,[responseVar '~' mdl],'Distribution','binomial','Link','logit') %Logisitc regression, excluding null responses
        ci=mm.coefCI;
        subplot(2,5,(j-1)*5+1)
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
        scatter(x,y,'filled');
        set(gca,'ColorOrderIndex',i)
        mm.plotPartialDependence(xvar)
        for k=1:size(ci,1)
            subplot(2,5,(j-1)*5+k+1)
            hold on
            c=mm.Coefficients.Estimate(k);
            bar(i,c)
            errorbar(i,c,c-ci(k,1),ci(k,2)-c,'k')
        end
        subplot(2,5,(j-1)*5+5)
        hold on
         set(gca,'ColorOrderIndex',i)
        y=splitapply(@(z) nanmean(z),aux.reactionTime,G);
        xx=splitapply(@(x) nanmean(x),mm.predict(aux),G); %Predicted responses
        [xx,idx]=sort(xx);
        scatter(y(idx),xx,'filled') % (true) RT vs. expected accuracy for that probe size
         set(gca,'ColorOrderIndex',i)
         auxT=(2*xx-1)./log(xx./(1-xx));
         auxT(xx==.5)=.5; %Forcing limit
         pp=polyfit(auxT(~isnan(auxT)),y(idx(~isnan(auxT))),1);
         rt=pp(1)*auxT+pp(2);
        plot(rt,xx) %EZ DDM prediction
        xlabel('RT')
        ylabel('%')
        subplot(2,5,(j-1)*5+4)
                hold on
         set(gca,'ColorOrderIndex',i)
         xx=splitapply(@(x) nanmean(x),aux.(xvar),G); %Probe size
         [xx,idx]=sort(xx);
         y=y(idx);
         scatter(xx,y,'filled')
         set(gca,'ColorOrderIndex',i)
        plot(xx,rt) %EZ DDM prediction
        ylabel('RT')
        xlabel('Probe size (mm/s)')
    end
end