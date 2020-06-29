function newValues = ChrTracer_ParameterGUI(currentValues,varargin)
%  newValues = ParameterMenu(currentValues,varargin)
% Creat a dialogue box that allows strings, numbers, and logicals in the
% structure "currentValues" to be changed

defaults = cell(0,3);
defaults(end+1,:) = {'title','string','Select Parameters'};
parameters = ParseVariableArguments(varargin,defaults,mfilename);

%% Main Function
newValues = currentValues; % default is no change (in case aborted)

parameterNames = fieldnames(currentValues);
parameterValues = struct2cell(currentValues);

% remove parameters (such as tables) that can't be (readily) displayed in a popup. 
isValid = cellfun(@(x) (isnumeric(x) | islogical(x) | ischar(x) | iscell(x)),parameterValues);
isValid = isValid & cellfun(@(x) size(x,1),parameterValues) < 2;
parameterValues = parameterValues(isValid);
parameterNames = parameterNames(isValid); 

% conver cell array of strings to comma separted list
isCellPar = cellfun(@iscell,parameterValues);
parameterValues(isCellPar) = cellfun(@CellToCSL,parameterValues(isCellPar),'UniformOutput',false) ;

% convert numerics and logicals to strings
isNumericPar = cellfun(@isnumeric,parameterValues);
isLogicalPar = cellfun(@islogical,parameterValues);
isChar = cellfun(@ischar,parameterValues);
isArrayPar = cellfun(@(x) size(x,2),parameterValues)>1  & ~isChar & ~isCellPar;
parameterValues(isNumericPar & ~isArrayPar) = cellfun(@num2str,parameterValues(isNumericPar & ~isArrayPar),'UniformOutput',false);
parameterValues(isLogicalPar & ~isArrayPar) = cellfun(@logical2str,parameterValues(isLogicalPar & ~isArrayPar),'UniformOutput',false);
parameterValues(isArrayPar & (isNumericPar | isLogicalPar)) = cellfun(@(x) ['[',num2str(x),']'], parameterValues(isArrayPar & (isNumericPar | isLogicalPar)), 'UniformOutput',false); % arrays need to be grouped in brackets to be parsed correctly

% could write a little loop that parses structures into component parts. 
isNotString = ~cellfun(@ischar,parameterValues);
parameterValues(isNotString) = [];
parameterNames(isNotString) = []; 
isNumericPar(isNotString) = [];
isLogicalPar(isNotString) = [];

% get user choices (these will be returned as strings)
choices = inputdlg(parameterNames,parameters.title,1,parameterValues);

% parse choices back into the desired formats
if ~isempty(choices)
    for i=1:length(parameterNames)
        if isNumericPar(i) || isArrayPar(i)
            newVal = str2num(choices{i}); %#ok<ST2NM>
        elseif isLogicalPar(i)
            if strcmp(choices{i},'true')
                newVal = true;
            else 
                newVal = false;
            end
        elseif isCellPar(i)
            newVal = strsplit(choices{i},',');
        else
            newVal = choices{i};
        end
        newValues.(parameterNames{i}) = newVal;
    end
end
