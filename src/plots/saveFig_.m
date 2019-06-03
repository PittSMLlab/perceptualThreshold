function saveFig_(h,dir,fileName,sizeFlag)
%saveFig saves figure h as .fig and .png for further reference
fName='OpenSans';
txt=findobj(h,'Type','Text');
set(txt,'FontName',fName);
ax=findobj(h,'Type','Axes');
set(ax,'FontName',fName);
for i=1:length(ax)
    ax(i).Title.FontWeight='normal';
end

if nargin<4 || isempty(sizeFlag)
set(h,'Units','Normalized','OuterPosition',[0 0 1 1])
end
fullName=[dir fileName];
if ~exist(dir,'dir')
    mkdir(dir)
end

%set(h,'Color','None')
%print(h, '-painters', '-dpng', '-r900', [fullName '.png']);
%savefig(h,[fullName '.fig'],'compact') ;
hgexport(h, [fullName '.png'], hgexport('factorystyle'), 'Format', 'png');
hgexport(h,[fullName '.eps'], hgexport('factorystyle'), 'Format', 'eps');
%hgexport(h,[fullName '.svg'], hgexport('factorystyle'), 'Format', 'svg');
end
