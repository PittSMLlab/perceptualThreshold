%%
addpath(genpath('../../src/'))
%% Load
b1=load('20180705T123618_PD_Baseline1.mat');
b2=load('20180705T124346_PD_Baseline2.mat');
a1=load('20180705T125434_PD_Adaptation.mat');
p1=load('20180705T131443_PD_PostAdaptation.mat');

%% Basic plot
figure; hold on;
t0=0;
for i=1:4
    switch i
        case 1
            aux=b1;
        case 2
            aux=b2;
        case 3
            aux=a1;
        case 4
            aux=p1;
    end
    t=aux.datlog.TreadmillCommands.read(:,4)+t0;
plot(t,aux.datlog.TreadmillCommands.read(:,1:2))
t1=aux.datlog.TreadmillCommands.sent(:,4)+t0;
plot(t1,aux.datlog.TreadmillCommands.sent(:,1:2),'ko')
t0=t(end);
end

%Fast task: need to quantify, subject seems to undershoot significantly
%Adaptation: even 25 strides in there is a huge change in behavior. WHY? a
%very fast process? or something about rapid changes is a more prominent
%signal? (is there a modeling difference between these two)?
%Adaptation: what is the most generous interpretation of the learning?
%60mm/s between first and last tests?
%PostAdaptation: what is the most generous interpretation of aftereffects?

%Is it worth it to start subjects at both +100 and +200? Do subjects
%correct differently for both? Should we also start subjects at +500 (where
%they are) to see what they do?
