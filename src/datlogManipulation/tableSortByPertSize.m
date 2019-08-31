function tbl = tableSortByPertSize(superSuperT)

ps=superSuperT.pertSize; %pertSize
nr=isnan(superSuperT.initialResponse); %No responses
cr=superSuperT.initialResponse==-sign(superSuperT.pertSize) & ~nr; %Correct responses
rt=superSuperT.reactionTime;
%rt=superSuperT.reactionStride;
pList=unique(ps);
for i=1:length(pList)
    relTrials=ps==pList(i);
    noResponseTrials =relTrials & nr;
    correctTrials=relTrials & cr & ~nr;
    incorrectTrials=relTrials & ~cr & ~nr;
   noResponseRate(i)=mean(nr(relTrials));
   correctRate(i)=mean(cr(relTrials) & ~nr(relTrials) & sign(pList(i))~=0);
   auxRT=rt(relTrials);
   meanRT(i)=nanmean(auxRT);
   medianRT(i)=nanmedian(auxRT);
   auxRT(isnan(auxRT))=Inf;
   medianRTprojected(i)=median(auxRT);
   correctMeanRT(i)=mean(rt(correctTrials));
   incorrectMeanRT(i)=mean(rt(incorrectTrials));
   correctMedianRT(i)=median(rt(correctTrials));
   incorrectMedianRT(i)=median(rt(incorrectTrials));
end

vn={'DeltaV','NRrate','accRate','meanRT','medRT','projMedRT','corMeanRT','incMeanRT','corMedRT','incMedRT'};
tbl=table(pList,noResponseRate',correctRate', meanRT',medianRT', medianRTprojected',correctMeanRT',incorrectMeanRT',correctMedianRT',incorrectMedianRT','VariableNames',vn);
end

