%% Generate accuracy plots for each subject
addpath(genpath('./'));
%% This file assumes the existence of csv files that summarize block results. If this is not the case, run processDatlogs.m
%%
baseDir='../data/AA0';
figDir='../fig/AA0';
for j=1:9 %All subjects
    exit=false;
    blockNo=0;
    subDir=[baseDir num2str(j) '/'];
    figDir1=[figDir num2str(j) '/'];
   while ~exit 
       blockNo=blockNo+1;
       try
          t=readtable([subDir 'Block' num2str(blockNo) '.csv']); 
          %fh=accuracyPlots(t);
          %saveFig(fh,subDir,['accBlock' num2str(blockNo)],0)
          if blockNo==1
              superT=t;
          else
             superT=cat(1,superT,t);
          end
       catch %Went through ALL blocks for this subject
           exit=true;
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
    fh=accuracyPlots(superSuperT);
    %saveFig(fh,'../fig/all/',['accuracyAll'],0)
%% Block 1 vs. block 2 accuracy, and block 2 vs. block 3 accuracy (when present)
%This (loosely) compares history before perturbation as a factor in performance

%% Block 1 vs block 3, and block 2 vs block 4 accuracy (when present)
%This assesses learning in the task

%% Dominance effects? Separate based on handedness? (footedness is R for all subjects)