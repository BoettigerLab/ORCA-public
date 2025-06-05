function structOut = ReadTextToStruct(textFile,varargin)
% ----------------------------------------------------------------------- %
% Read data in text file in which each field is indicated as a separate
% line organized fieldName = fieldValue.
% ----------------------------------------------------------------------- %
% Alistair Boettiger
% CC BY Aug 8 2017
% ----------------------------------------------------------------------- %


%%
defaults = cell(0,3);
defaults(end+1,:) = {'verbose','boolean',false};
pars = ParseVariableArguments(varargin,defaults,mfilename); %  

textFile = regexprep(textFile,{'.i4d','.i5d'},'.inf');
if pars.verbose
   disp(['reading ',textFile]); 
end

% Open file
fid = fopen(textFile);
if fid == -1
    error([textFile ' is not a valid .inf file']);
end

textFields = {'encoding','dataType','description','dataFolder','notes'};

% Read file
while ~feof(fid)
    textIn = fgetl(fid);
    % textIn = text{1}
    fieldParts = strsplit(textIn,' = ');
    fieldName = fieldParts{1};
    fieldValue = fieldParts{2};
    if ~strcmp(fieldName,textFields)
        fieldValue = str2num(fieldValue); %#ok<ST2NM>
    end
    structOut.(fieldName) = fieldValue;
    
end
fclose(fid);

