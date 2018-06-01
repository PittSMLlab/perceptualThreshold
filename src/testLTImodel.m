%% Testing some linear modeling


%% Loading data
speeds=cell(9,4);
baseDir='../data/AA0';
figDir='../fig/AA0';
for sub=1:9
    subDir=[baseDir num2str(sub) '/'];
    for blockNo=1:4
        try
            speeds{sub,blockNo}=csvread([subDir 'Block' num2str(blockNo) '_speedDiff.csv']);
        end
    end
end

%% Plotting some data
fun=@nanmean;
s=reshape(speeds,9*4,1);
matSpeeds=cell2mat(s(cellfun(@(x) ~isempty(x),s)));
pp=unique(matSpeeds(:,1));
cmap=parula(length(pp)); %parula, jet, redbluecmap
cmap=cmap*.8;
figure
m=splitapply(fun,matSpeeds./sign(matSpeeds(:,1)),findgroups(abs(matSpeeds(:,1))));
figure
subplot(3,1,1)
hold on
set(gca,'ColorOrder',cmap)
plot(m')
ylabel('Speed diff')
title('Speed diff. evolution')
xlabel('Stride #')
subplot(3,1,2)
hold on
set(gca,'ColorOrder',cmap)
scaledMatSpeeds=matSpeeds./matSpeeds(:,1);
sm=splitapply(fun,scaledMatSpeeds,findgroups(matSpeeds(:,1)));
plot(sm')
ylabel('Normalized Speed diff (as % of initial)')
title('Normalized')
xlabel('Stride #')
subplot(3,1,3)
hold on
set(gca,'ColorOrder',cmap)
scaledMatSpeeds=matSpeeds./matSpeeds(:,1);
sm=splitapply(fun,scaledMatSpeeds,findgroups(matSpeeds(:,1)));
for i=1:size(sm,1)
    plot(sm(i,find(sm(i,:)==1,1,'last'):end))
end
ylabel('Normalized Speed diff (as % of initial)')
title('Normalized+aligned')
xlabel('Stride #')
%scatter(24*ones(size(m)),m,20,pp,'filled')
%% Simulation of models
N=200;
Dv=zeros(N,1);
action=zeros(N,1);
changeAction=zeros(N,1);
perceptualAction=zeros(N,1);
perceptualSymmetry=zeros(N,1);
pertSize=300;
Dv(2)=pertSize;
for i=2:N
    perceptualSymmetry(i)=.96*perceptualSymmetry(i-1)+.02*Dv(i);
    changeAction(i)=.5*changeAction(i)+.2*(Dv(i)-Dv(i-1));
    perceptualAction(i)=.2*(Dv(i)-perceptualSymmetry(i));
    action(i)=changeAction(i)+perceptualAction(i);
    Dv(i+1)=Dv(i)-action(i); %Closed-loop
    %Dv(i+1)=Dv(i); %Open-loop
end

figure
plot(Dv,'DisplayName','\Delta v','LineWidth',2)
hold on
plot(perceptualSymmetry,'DisplayName','Perception','LineWidth',2)
plot(action,'DisplayName','Action taken','LineWidth',2)
axis([0 30 0 pertSize])