%% main

%% Load
load('../threshProfile_290816.mat');
vD=velR-velL;
Ts=1; %Step-time
RTOt=2*Ts*[0:(length(vD)-1)]'; %Stride times
RHSt=RTOt+Ts/3;
LTOt=RTOt+Ts;
LHSt=RTOt+4*Ts/3;
aux=isnan(velR);


%% Generate pertSchedule
velocityScheduleL=timeseries(velR(~isnan(velR))*1000,LTOt(~isnan(velR)));
velocityScheduleR=timeseries(velL(~isnan(velL))*1000,RTOt(~isnan(velL)));
goCue=timeseries(double(aux),RTOt);
%% Run sim
%set_param('fullLoop', 'StopTime', num2str(length(vD)-1))
mu=0; %Noise bias
sigma=50; %Noise std
a=0; %Pole-position for perception filter
map=1; %Static-map btw perception and action. 1=linear, 2=threshold, 3=quadratic
k=5/200; %Gain for all maps.
th=20; %Threshold for thresholded map

simout=sim('./fullLoop.slx','StopTime',num2str(length(vD)-1));
beltSpeedDiff=simout.get('beltSpeedDiff');
keypresses=simout.get('keypresses');
%% Generate fake datlog
m=1050;

beltSpeedDiff.data=reshape(beltSpeedDiff.data(:),length(beltSpeedDiff.time),1);
sentT=sort([RTOt; LTOt],'ascend');
datlog.TreadmillCommands.sent=[m+beltSpeedDiff.Data/2 m-beltSpeedDiff.Data/2 zeros(size(beltSpeedDiff.Time)) beltSpeedDiff.Time];
datlog.TreadmillCommands.read=datlog.TreadmillCommands.sent;

datlog.speedprofile.velR=velR*1000;
datlog.speedprofile.velL=velL*1000;

datlog.stepdata.LTOdata=[nan(length(LTOt),3) LTOt]; 
datlog.stepdata.LHSdata=[nan(length(RTOt),3) LHSt];
datlog.stepdata.RTOdata=[nan(length(RTOt),3) RTOt];%First event
datlog.stepdata.RHSdata=[nan(length(RTOt),3) RHSt];
datlog.stepdata.RHSdata=datlog.stepdata.RHSdata(1:end-1,:); %REmoving last rHS, to be consistent with the real GUI performance

datlog.audioCues.start=goCue.Time(diff(goCue.Data)==1);
datlog.audioCues.stop=goCue.Time(diff(goCue.Data)==-1);



datlog.addLog.keypress=cell(sum(keypresses.data(:)),2);
counter=0;
for i=1:length(keypresses.time)
   Rpress=keypresses.data(i,1);
   Lpress=keypresses.data(i,2);
   datlog.addLog.keypress(counter+[1:Rpress],1)={'rightarrow'};
   for j=1:Rpress
    datlog.addLog.keypress(counter+j,2)={Ts*(keypresses.time(i)+j/(Rpress+1))};
   end
   counter=counter+Rpress;
   datlog.addLog.keypress(counter+[1:Lpress],1)={'leftarrow'};
   for j=1:Lpress
    datlog.addLog.keypress(counter+j,2)={Ts*(keypresses.time(i)+j/(Lpress+1))};
   end
   counter=counter+Lpress;
end
datlog.framenumbers.data=[0 0 0];

save([num2str(now) 'simDatlog.mat'], 'datlog')
save(['lastSimDatlog.mat'], 'datlog')
%% Run analysis of datlog
datlogAnalysis(datlog)