function fitEZ(dataTable)
%Fits an EZ-difussion model, as in Wagenmakers et al. 2007
%Requires that the table contains a list of trials and the following fields
%for each trial: perturbationSize (type/size) (taken as categorical, just to fit a
%different drift rate for each), reaction time, correct response (binary).

pp=unique(dataTable.perturbationSize);


end