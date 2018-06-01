function [trialData,strideData]=datlogSummarize(datlog)
goodOnly=0;
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
        firstResponse(i)=NaN; %L or R press
    end
end

%% Returns a matrix with the following columns:
%trial #, date+time, startTime, endTime, pertSize, reactionTime, reactionStride, initialResponse, steadyState, Lclicks, Rclicks, 
Ntrials=length(audiostart); %Number of trials
Nblocks =Ntrials/24; %Each block is 24 trials long
Nstrides= Ntrials*46 +10*Nblocks; %There should be 46 strides in each trial (24 active, 22 inactive), and 10 extra strides at end of each block
startTime=audiostart;
endTime=audiostop;
date=datenum(datlog.buildtime);
pertSize;
reactionTime;
reactionStride;
initialResponse=reactionSign;
%Reshape into blocks:
list={'accumLresponses','accumRresponses','responseTime','RTOt','RHSt','LTOt','LHSt','vRload','vLload','vR','vL','vD','firstResponse'}; %tasks start at primary toe-off (pTO), which is always RTOt by design
LHSt=cat(1,LHSt,ones(Nblocks,1)); %Adding final LHSt event to each block (Trial ends before that)
for k=1:length(list)
    eval(['aux=' list{k} ';']);
    
    aux=reshape(aux,46*24+10,Nblocks);
    aux=reshape(aux(2:(46*24+1),:),46,24,Nblocks);
    aux=aux(21:45,:,:);
    eval([list{k} '=aux;'])
end
Lclicks=squeeze(nansum(accumLresponses,1));
Rclicks=squeeze(nansum(accumRresponses,1));
Lclicks=Lclicks(:);
Rclicks=Rclicks(:);
speedDiffAtStart=vD(1,:)';
velRprofile=vRload;
velLprofile=vLload;
speedRactual=vR;
speedLactual=vL;
lastSpeedDiff=vD(end,:)';
date=date*ones(size(startTime));
trialData=table(date,startTime,endTime,pertSize,reactionTime,reactionStride,initialResponse,Lclicks,Rclicks,lastSpeedDiff);

%% Returns a tensor of Ntrials x 25 strides x M fields, containing:
%startTime,endTime,Lclicks,Rclicks, initialClick, timeToFirstClick,

strideData.startTime=RTOt(:,:)';
strideData.endTime=LHSt(:,:)'; %Not quite, but first approximation
strideData.Lclicks=accumLresponses(:,:)';
strideData.Rclicks=accumRresponses(:,:)';
strideData.initialClick=firstResponse(:,:)';
strideData.timeToFirstClick=responseTime(:,:)';
strideData.speedDiff=vD(:,:)';
%strideData=cat(3,startTime,endTime,Lclicks,Rclicks,initialClick,timeToFirstClick);


end

