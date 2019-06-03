function extendedPanelWidth(fh,w)
if nargin<2
    w=.1; %10%
end
ph=findobj(fh,'Type','Axes');
for i=1:length(ph)
   ph(i).Position=ph(i).Position+(ph(i).Position-[.5 0 0 0]).*[w 0 w 0]; 
end
end

