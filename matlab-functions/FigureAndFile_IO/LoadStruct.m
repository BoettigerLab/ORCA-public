function [structName,origStructName] = LoadStruct(structMatFile,varargin)
% load a structure saved by SaveStruct 
% FUNCTION Requires testing
% LoadStruct(structMatFile);

% Approach: convert cells to table and table to structure

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'requires a structure');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
varInput = varargin(2:end);
parameters = ParseVariableArguments(varInput, defaults, mfilename);

%% main function
% structMatFile = '\\cajal\TSTORMdata\AdditionalAnalysis\141101_AnnotatedMosaics\words.mat';

structName = 'savedStructure';
localVars = who;
load(structMatFile);
allVars = who;
[~,idx] = intersect(allVars,localVars);
loadedVars = allVars(idx); 
origStructName = structName; 

T = table(loadedVars{:});
eval([structName ' = table2struct(T)']);


if parameters.verbose
    disp(['loaded ',structMatFile]);
end

    