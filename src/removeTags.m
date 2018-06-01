function [stringNoTags] = removeTags(string,removeBlankLinesFlag)
%removeTags takes a string and removes tags of the form <SOME_TAG> (html style)
%from all its ocurrences
if nargin<2 || isempty(removeBlankLinesFlag)
    removeBlankLinesFlag=true;
end
stringNoTags=regexprep(string,'<[^>]*>','');
stringNoTags=regexprep(stringNoTags,'~','='); %Also removing '~'

if removeBlankLinesFlag
    stringNoTags = regexprep(stringNoTags, '\n\n+', '\n');
end
end

