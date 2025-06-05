function SaveStructAsText(datInfo,fullSaveName,varargin)
% ----------------------------------------------------------------------- %
% Save structure as a text file
% fieldName = fieldValue \newline
% 
% ----------------------------------------------------------------------- %
% Alistair Boettiger
% CC BY Aug 8 2017
% ----------------------------------------------------------------------- %

defaults = cell(0,3);
defaults(end+1,:) = {'verbose','boolean',false};
pars = ParseVariableArguments(varargin,defaults,mfilename); %  

fullSaveName = regexprep(fullSaveName,{'.i4d','.i5d'},'.inf');
fid = fopen(fullSaveName,'w+');
names = fields(datInfo);
values = struct2cell(datInfo);
numVals = length(names);
for i=1:numVals
    if ~ischar(values{i}) && ~isstring(values{i})
        value = num2str(values{i});
    else
        value = values{i};
    end
    stringOut = [names{i},' = ',value];
    fprintf(fid,'%s\r\n',stringOut);
end

fclose(fid); 
if pars.verbose
   disp(['wrote ',fullSaveName]); 
end