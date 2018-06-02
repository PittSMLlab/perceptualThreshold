addpath(genpath('./'));
%% This file assumes the existence of csv files that summarize block results. If this is not the case, run processDatlogs.m
%%
baseDir='../data/AA0';
figDir='../fig/AA0';
    subID=[];
for j=1:9 %All subjects
    exit=false;
    blockNo=0;
    subDir=[baseDir num2str(j) '/'];
    figDir1=[figDir num2str(j) '/'];
   while ~exit 
       blockNo=blockNo+1;
       try
          t=readtable([subDir 'Block' num2str(blockNo) '.csv']);
          if false %blockNo==3 %To drop block 3 from subjects with odd number of blocks
              try %Try loading block 4, if fails (no block 4), dont use block 3 either
                  t4=readtable([subDir 'Block' num2str(blockNo+1) '.csv']);
              catch
                  exit=true;
                  break
              end
          end
       catch %Went trhough all blocks for this sub
           exit=true;
           break
       end
          %fh=accuracyPlots(t);
          %saveFig(fh,subDir,['accBlock' num2str(blockNo)],0)
          taux=table(j.*ones(size(t.date,1),1), blockNo*ones(size(t.date,1),1),'VariableNames',{'subID','blockNo'});
          t=cat(2,t,taux);
          if blockNo==1
              superT=t;
          else
             superT=cat(1,superT,t);
          end

   end
    %fh=accuracyPlots(superT);
    %fh.Name=['Accuracy for subject AA0' num2str(j)];
    %saveFig(fh,figDir1,['accuracyAA0' num2str(j)],0)
    if j==1
        superSuperT=superT;
    else
        superSuperT=cat(1,superSuperT,superT);
    end
end