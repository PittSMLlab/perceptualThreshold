function [pSize,driftRate,threshold,delay]=fitEZ(dataTable)
%Fits an EZ-difussion model, as in Wagenmakers et al. 2007
%Requires that the table contains a list of trials and the following fields
%for each trial: pertSize (type/size) (taken as categorical, just to fit a
%different drift rate for each), reaction time, correct response (binary),
%no response (binary).


pSize=unique(dataTable.pertSize);
threshold=nan(size(pSize));
driftRate=nan(size(pSize));
delay=nan(size(pSize));
for i=1:length(pSize)
    trials=dataTable(dataTable.pertSize==pSize(i),:);
        RT=trials.reactionTime;
        correct=trials.correctResponse;
        NR=trials.noResponse;
    if pSize(i)~=0
        meanCorrectRT=mean(RT(correct & ~NR));
        varCorrectRT=var(RT(correct & ~NR));
        meanAccuracy=mean(correct(~NR));
        if meanAccuracy==1
            meanAccuracy=1-.5/numel(RT); %Adding half a wrong trial for computation sake
        end
        y=log(1/meanAccuracy -1); %This is minus the product of drift and threshold
        z=varCorrectRT*(1+exp(y))^2/(1-exp(2*y)+2*y*exp(y)); %This is half the threshold over drift^3
        x=-y/z; %This is twice drift^4
        driftRate(i)=sqrt(sqrt(x/2));
        threshold(i)=-y/driftRate(i);
        delay(i)=meanCorrectRT-.5*threshold(i)/driftRate(i)*(2*meanAccuracy-1);
    else
       varCorrectRT=var(RT(~NR)); %tAKING ALL TRIALS
       meanCorrectRT=mean(RT(~NR));
       driftRate(i)=0;
       threshold(i)=sqrt(sqrt(varCorrectRT*24));
       delay(i)=meanCorrectRT-threshold(i)^2/4;
    end
end
    
    


end