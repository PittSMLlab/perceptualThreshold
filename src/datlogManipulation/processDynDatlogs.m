%% load datlog
addpath(genpath('../src/datlogManipulation/'))
dataDir='./';
subList=dir([dataDir 'AB*']);
%% Load fast & slow baseline trials
for j=1:length(subList)
    %Normal (fast) baseline:
    base2trials=dir([dataDir subList(j).name '/*Baseline2.mat']);
    for k=1:length(base2trials)
        load([dataDir subList(j).name '/' base2trials(k).name],'datlog');
        [trialData,strideData]=datlogSummarizeFast(datlog);
        writetable(trialData,[dataDir subList(j).name '/Baseline' num2str(k) '.csv']);
        fF=fields(strideData);
        for i=1:length(fF)
            csvwrite([dataDir subList(j).name '/Baseline' num2str(k) '_' fF{i} '.csv'],strideData.(fF{i}));
        end
    end
    %Slow baseline:
    base1trials=dir([dataDir subList(j).name '/*Baseline1.mat']);
    load([dataDir subList(j).name '/' base1trials(1).name],'datlog');
    [trialData,strideData]=datlogSummarizeFast(datlog);
    writetable(trialData,[dataDir subList(j).name '/BaselineSlow.csv']);
    fF=fields(strideData);
    for i=1:length(fF)
        csvwrite([dataDir subList(j).name '/BaselineSlow_' fF{i} '.csv'],strideData.(fF{i}));
    end
    %Adapt data:
    trials=dir([dataDir subList(j).name '/*_Adaptation.mat']);
    load([dataDir subList(j).name '/' trials(1).name],'datlog');
    [trialData,strideData]=datlogSummarizeFast(datlog);
    writetable(trialData,[dataDir subList(j).name '/Adaptation.csv']);
    fF=fields(strideData);
    for i=1:length(fF)
        csvwrite([dataDir subList(j).name '/Adaptation_' fF{i} '.csv'],strideData.(fF{i}));
    end
    %Postadapt data:
        trials=dir([dataDir subList(j).name '/*_PostAdaptation.mat']);
    load([dataDir subList(j).name '/' trials(1).name],'datlog');
    [trialData,strideData]=datlogSummarizeFast(datlog);
    writetable(trialData,[dataDir subList(j).name '/PostAdaptation.csv']);
    fF=fields(strideData);
    for i=1:length(fF)
        csvwrite([dataDir subList(j).name '/PostAdaptation_' fF{i} '.csv'],strideData.(fF{i}));
    end
end
