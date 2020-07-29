function regData = LoadRegDataCSV(saveFolder,varargin)
% load all *_regData.csv files in folder as a cell array of structures
%
% % see also LoadRegDataCSV (combine?)

% defaults = cell(0,3);
% defaults(end+1,:) = {'expectFOVs','nonnegative',0};  
% pars = ParseVariableArguments(varargin,defaults,mfilename);

regDataFiles = ls([saveFolder,filesep,'*_regData.csv']);
if ~isempty(regDataFiles)
    regDataFiles = cellstr(regDataFiles);
end
numFOVs = length(regDataFiles);
regData = cell(numFOVs,1); 

% main function
for f=1:numFOVs
    % tableName = [saveFolder,'fov',num2str(f,'%03d'),'_regData.csv'];
    tableName = [saveFolder,filesep,regDataFiles{f}];
    tableData = readtable(tableName);
    regData{f} = table2struct(tableData);
end