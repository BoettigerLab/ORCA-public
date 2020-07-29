[numHits, maxLength] = ParseBLAST(blastData, vargin)
% ------------------------------------------------------------------------
% [numHits, maxLength] = ParseBLAST(blastData)
% This function returns a quick parse of a blast data structure.
%--------------------------------------------------------------------------
% Necessary Inputs
% fastaFile: String to a fasta file
%--------------------------------------------------------------------------
% Outputs
% dataBaseOut: Path to the created dataBase 
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1 || ~isstruct(blastData)
    error('matlabFunctions:invalidArguments', 'A valid blast data structure is required.');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Loop
% -------------------------------------------------------------------------
numHits = arrayfun(@(x) length(x.Hits), blastData);
hits = arrayfun(@(x) [x.Hits.HSPs], blastData, 'UniformOutput', false, ');
maxHits = arrayfun(@(x) max(arrayfun(@(y) abs(diff(y.HSPs(1).QueryIndices)), x.Hits),'UniformOutput', false), blastData);
