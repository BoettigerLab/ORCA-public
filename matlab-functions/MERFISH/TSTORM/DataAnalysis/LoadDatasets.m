function datasets = LoadDatasets(dataPath, varargin)
% ------------------------------------------------------------------------
% datasets = LoadDatasets(dataPath)
% This function loads an xls spreadsheet containing information on
% multiFISH experiments to analyze.
%--------------------------------------------------------------------------
% Necessary Inputs
% dataPath/(\\cajal\TSTORMdata\Datasets.xls): The path to the xls/xlsx
% spreadsheet to be loaded.
%--------------------------------------------------------------------------
% Outputs
% datasets/A structure array containing the information necessary to
% analyze the desired data
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Alistair Boettiger, Jeffrey Moffitt
% jeffmoffitt@gmail.com
% October 12, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Define standard fields
% -------------------------------------------------------------------------
standardField = {'dataPath', 'numWords', 'libraryNum', 'experimentNum', ...
    'hybOrder', 'numHybs', 'codebookPath', 'FPKMDataPath', 'correctable', ...
    'listTag', 'imageTag'};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    dataPath = '\\cajal\TSTORMdata\Datasets.xls';
end

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true}; 

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Load Data
% -------------------------------------------------------------------------
try
    tableData = readtable(dataPath); % xls or csv
    datasets = table2struct(tableData);
catch er
    warning(er.getReport);
    error('matlabFunctions:invalidArguments', ['data path ',dataPath,' is invalid']);
end

% -------------------------------------------------------------------------
% Check for the necessary fields
% -------------------------------------------------------------------------
if ~isempty(setdiff(standardField, fields(datasets)))
    warning('matlabFunctions:missingFields', 'Some standard fields were not found');
end

% -------------------------------------------------------------------------
% Remove empty entries
% -------------------------------------------------------------------------
if ismember('dataPath', fields(datasets))
    goodInds = cellfun(@(x)~isempty(x), {datasets.dataPath});
    datasets = datasets(goodInds);
end

% -------------------------------------------------------------------------
% Display Info on Files To Be Analyzed
% -------------------------------------------------------------------------
if parameters.verbose
    display('--------------------------------------------------------------');
    display(['Found ' num2str(length(datasets)) ' data sets to analyze']);
    for i=1:length(datasets)
        display(['   ' datasets(i).dataPath]);
    end
end


