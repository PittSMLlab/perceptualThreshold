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

%% Do EZ modeling:

Nsubs=10;
figure; hold on; 
for i=1:Nsubs
    aux=superSuperT(superSuperT.subID==i,:); %Single subj
    mm=fitglm(X,'leftResponse~ID*pertSize+pertSize*prevSize+blockNo*pertSize+lastSpeedDiff','Distribution','binomial'); %Logisitc regression, excluding null responses

    %Get data:
    x=aux.pertSize;
    y=aux.initialResponse==-1+.5*isnan(aux.initialResponse);
    acc=aux.initialResponse==-sign(aux.pertSize); %Correct responses
    nr=isnan(aux.initialResponse); %No responses
    nt=x==0; %Null trials
    t=aux.reactionTime;
    for j=1:2
        switch j
            case 1
            x1=x(~nr);
            y=y(~nr);
            fixedBias=[];
            case 2
            x1=abs(x);
            y=acc(~nr & ~nt);
            x1=x1(~nr & ~nt);
            fixedBias=0;
        end
    %[params, predictedY, Likelihood] = fitGenPsycho(x1,y,'MLE',fixedBias); %fit optimal model to signed data
    [params, predictedY, Likelihood] = fitPsycho(x1,y,'MLE',fixedBias); %fit optimal model to signed data
    
    %Plot:
    subplot(2,5,(j-1)*5+1)
    hold on;
    G=findgroups(x1); 
    plot(splitapply(@(z) mean(z),x1,G),splitapply(@(z) mean(z),y,G),'o'); %Mean responses
    plot(sort(x1),sort(predictedY),'Color',.5*ones(1,3));
    text(max(x1)+10,max(predictedY),num2str(i))
    grid on
    subplot(2,5,(j-1)*5+2) %Bias
    hold on
    bar(i,params(1))
    ylabel('Bias')
    subplot(2,5,(j-1)*5+3) %Slope
    hold on
    bar(i,params(2))
    ylabel('Noise')
    subplot(2,5,(j-1)*5+4) %Alpha
    hold on
    bar(i,params(end))
    subplot(2,5,(j-1)*5+5) %rt
    hold on
    G=findgroups(x);
    scatter(splitapply(@(z) mean(z),t,G),splitapply(@(z) mean(z),acc,G))
    xlabel('RT')
    ylabel('Acc')
    end
end
  mm=fitglm(X,'leftResponse~ID*pertSize+pertSize*prevSize+blockNo*pertSize+lastSpeedDiff','Distribution','binomial'); %Logisitc regression, excluding null responses

  
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