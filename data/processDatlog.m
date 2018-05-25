%% load datlog

%load('18_Aug_2016_16_56_27_threshProfile_160816.mat')
%load('16_Aug_2016_18_10_26_threshProfile_160816.mat')
%load('23_Aug_2016_13_22_11_threshProfile_230816.mat')
%load('./sim/lastSimDatlog.mat')
%load('26_Aug_2016_18_01_14_threshProfile_260816.mat')
%profileName=['./AA01/20160830T200954_threshProfile_300816.mat'];
%profileName=['./AA01/20160830T204312_threshProfile_300816opp.mat'];


clear profile

profile{1}{1}=['../data/AA01/20160830T200954_threshProfile_300816.mat'];
profile{1}{2}=['../data/AA01/20160830T204312_threshProfile_300816opp.mat'];

profile{2}{1}=['../data/AA02/20160920T180946_threshProfile_300816.mat'];
profile{2}{2}=['../data/AA02/20160920T184328_threshProfile_300816opp.mat'];
profile{2}{3}=['../data/AA02/20160920T192048_threshProfile_300816.mat'];
profile{2}{4}=['../data/AA02/20160920T195216_threshProfile_300816opp.mat'];

profile{3}{1}=['../data/AA03/20160922T191249_threshProfile_300816.mat'];
profile{3}{2}=['../data/AA03/20160922T194458_threshProfile_300816opp.mat'];
%profile{3}{3}=['../data/AA03/20160922T201415_threshProfile_300816.mat'];

profile{4}{1}=['../data/AA04/20161004T185906_threshProfile_300816.mat'];
profile{4}{2}=['../data/AA04/20161004T192712_threshProfile_300816opp.mat'];
profile{4}{3}=['../data/AA04/20161004T195550_threshProfile_300816.mat'];
profile{4}{4}=['../data/AA04/20161004T202348_threshProfile_300816opp.mat'];

profile{5}{1}=['../data/AA05/20161006T184845_threshProfile_300816.mat'];
profile{5}{2}=['../data/AA05/20161006T191454_threshProfile_300816opp.mat'];
profile{5}{3}=['../data/AA05/20161006T194202_threshProfile_300816.mat'];
profile{5}{4}=['../data/AA05/20161006T200843_threshProfile_300816opp.mat'];

profile{6}{1}=['../data/AA06/20161129T134420_threshProfile_300816.mat'];
profile{6}{2}=['../data/AA06/20161129T141117_threshProfile_300816opp.mat'];
profile{6}{3}=['../data/AA06/20161129T144041_threshProfile_300816.mat'];
profile{6}{4}=['../data/AA06/20161129T150649_threshProfile_300816opp.mat'];

profile{7}{1}=['../data/AA07/20180228T154736_threshProfile_300816.mat'];
profile{7}{2}=['../data/AA07/20180228T160815_threshProfile_300816opp.mat'];

profile{8}{1}=['../data/AA08/20180319T161447_threshProfile_300816.mat'];
profile{8}{2}=['../data/AA08/20180319T164257_threshProfile_300816opp.mat'];
profile{8}{3}=['../data/AA08/20180319T170901_threshProfile_300816.mat'];
profile{8}{4}=['../data/AA08/20180319T173454_threshProfile_300816opp.mat'];

profile{9}{1}=['../data/AA09/20180321T132152_threshProfile_300816.mat'];
profile{9}{2}=['../data/AA09/20180321T134649_threshProfile_300816opp.mat'];
%profile{9}{3}=['../data/AA09/20180321T141011_threshProfile_300816.mat'];

%% Process each block, generate tables, save as csv
for i=1:length(profile)
    for j=1:length(profile{i})
        load(profile{i}{j},'datlog');
        [trialData,strideData]=datlogSummarize(datlog);
        rootDir=profile{i}{j}(1:13);
        writetable(trialData,[rootDir 'Block' num2str(j) '.csv']);
        fF=fields(strideData);
        for k=1:length(fF)
            csvwrite([rootDir 'Block' num2str(j) '_' fF{k} '.csv'],strideData.(fF{k}));
        end
    end
end
%%
% close all
% for sub=7:9
% %% 1: individual blocks:
% for i=1:length(profile{sub})
%     load(profile{sub}{i})
%     fh=datlogAnalysis(datlog,0);
%     saveName=profile{sub}{i};
%     saveFig(fh,'./',saveName)
% end
% %% 2: all blocks together:
% [datlog] = catDatlogs(profile{sub});
% %% 
% psychoLapse=1;
% for goodOnly=[-1:1]
% fh=datlogAnalysis(datlog,goodOnly,psychoLapse);
% saveName=[profile{sub}{1}(3:6) 'All'];
%     if goodOnly==1
%         saveName=[saveName '_good'];
%     elseif goodOnly==-1
%         saveName=[saveName '_bad'];
%     end
%     if psychoLapse==1
%         saveName=[saveName 'Sat'];
%     end
% saveFig(fh,profile{sub}{1}(1:7),saveName)
% end
% end

% %% All subs together:
% allProfile={};
% for sub=1:length(profile)
%    allProfile=[allProfile profile{sub}]; 
% end
% [datlog] = catDatlogs(allProfile);
% %
% psychoLapse=1;
% for goodOnly=0%[-1:1]
% fh=datlogAnalysis(datlog,goodOnly,psychoLapse);
% saveName=['All_upto' datestr(datevec(now),30)];
%     if goodOnly==1
%         saveName=[saveName '_good'];
%     elseif goodOnly==-1
%         saveName=[saveName '_bad'];
%     end
%     if psychoLapse==1
%         saveName=[saveName 'Sat'];
%     end
% saveFig(fh,'./',saveName)
% end
% 
% %% All subs together, discard first block:
% allProfile={};
% for sub=1:length(profile)
%    allProfile=[allProfile profile{sub}(2:end)]; 
% end
% [datlog] = catDatlogs(allProfile);
% %
% psychoLapse=1;
% for goodOnly=0%[-1:1]
% fh=datlogAnalysis(datlog,goodOnly,psychoLapse);
% saveName=['All_noFirstBlock_upto' datestr(datevec(now),30)];
%     if goodOnly==1
%         saveName=[saveName '_good'];
%     elseif goodOnly==-1
%         saveName=[saveName '_bad'];
%     end
% saveFig(fh,'./',saveName)
% end
% 
% %% All subs together, second block only:
% allProfile={};
% for sub=1:length(profile)
%    allProfile=[allProfile profile{sub}(2)]; 
% end
% [datlog] = catDatlogs(allProfile);
% %%
% for goodOnly=[-1:1]
% fh=datlogAnalysis(datlog,goodOnly);
% saveName=['All_onlySecondBlock_upto' datestr(datevec(now),30)];
%     if goodOnly==1
%         saveName=[saveName '_good'];
%     elseif goodOnly==-1
%         saveName=[saveName '_bad'];
%     end
% saveFig(fh,'./',saveName)
% end
% 
% %% All subs together, even block only:
% allProfile={};
% for sub=1:length(profile)
%    allProfile=[allProfile profile{sub}(2:2:end)]; 
% end
% [datlog] = catDatlogs(allProfile);
% %
% psychoLapse=1;
% for goodOnly=0%[-1:1]
% fh=datlogAnalysis(datlog,goodOnly,psychoLapse);
% saveName=['All_evenBlocks_upto' datestr(datevec(now),30)];
%     if goodOnly==1
%         saveName=[saveName '_good'];
%     elseif goodOnly==-1
%         saveName=[saveName '_bad'];
%     end
% saveFig(fh,'./',saveName)
% end
% 
% %% All subs together, odd block only:
% allProfile={};
% for sub=1:length(profile)
%    allProfile=[allProfile profile{sub}(1:2:end)]; 
% end
% [datlog] = catDatlogs(allProfile);
% %
% psychoLapse=1;
% for goodOnly=0%[-1:1]
% fh=datlogAnalysis(datlog,goodOnly,psychoLapse);
% saveName=['All_oddBlocks_upto' datestr(datevec(now),30)];
%     if goodOnly==1
%         saveName=[saveName '_good'];
%     elseif goodOnly==-1
%         saveName=[saveName '_bad'];
%     end
% saveFig(fh,'./',saveName)
% end