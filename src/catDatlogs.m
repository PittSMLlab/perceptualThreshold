function [datlog] = catDatlogs(datlogNameList)

load(datlogNameList{1}); %First datlog
oTime=datevec(datlog.buildtime);
for i=2:length(datlogNameList)
    %Load:
   new=load(datlogNameList{i});
   %Do the cat-ting:
   %First: find relative buildtime from original datlog to new datlog
   nTime=datevec(new.datlog.buildtime);
   dT=etime(nTime,oTime);
   dT=datlog.framenumbers.data(end,3)+500; %Fake time gap
   datlog.framenumbers.data=[datlog.framenumbers.data; new.datlog.framenumbers.data(:,1:2) new.datlog.framenumbers.data(:,3)+dT];
   datlog.stepdata.RHSdata=[datlog.stepdata.RHSdata; new.datlog.stepdata.RHSdata(:,1:3) new.datlog.stepdata.RHSdata(:,4)+dT];
   datlog.stepdata.LHSdata=[datlog.stepdata.LHSdata; new.datlog.stepdata.LHSdata(:,1:3) new.datlog.stepdata.LHSdata(:,4)+dT];
   datlog.stepdata.RTOdata=[datlog.stepdata.RTOdata; new.datlog.stepdata.RTOdata(:,1:3) new.datlog.stepdata.RTOdata(:,4)+dT];
   datlog.stepdata.LTOdata=[datlog.stepdata.LTOdata; new.datlog.stepdata.LTOdata(:,1:3) new.datlog.stepdata.LTOdata(:,4)+dT];
   datlog.audioCues.stop=[datlog.audioCues.stop; new.datlog.audioCues.stop+dT];
   datlog.audioCues.start=[datlog.audioCues.start; new.datlog.audioCues.start+dT];
   datlog.speedprofile.velL=[datlog.speedprofile.velL; new.datlog.speedprofile.velL];
   datlog.speedprofile.velR=[datlog.speedprofile.velR; new.datlog.speedprofile.velR];
   datlog.TreadmillCommands.read=[datlog.TreadmillCommands.read; new.datlog.TreadmillCommands.read(:,1:end-1) new.datlog.TreadmillCommands.read(:,4)+dT];
   datlog.TreadmillCommands.sent=[datlog.TreadmillCommands.sent; new.datlog.TreadmillCommands.sent(:,1:end-1) new.datlog.TreadmillCommands.sent(:,4)+dT];
   new.datlog.addLog.keypress(:,2)=mat2cell(cell2mat(new.datlog.addLog.keypress(:,2))+dT,ones(size(new.datlog.addLog.keypress,1),1),1); 
   datlog.addLog.keypress=[datlog.addLog.keypress;new.datlog.addLog.keypress];
end

end

