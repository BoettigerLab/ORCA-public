function rnaCnts = WordsPerCell(words,libGenes,varargin) 
% rnaCnts = WordsPerCell(words,libGenes,varargin) 
% 

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'perfectMatch', 'boolean', false};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'a nx3 data array is required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

%%

if parameters.perfectMatch
    keep = [words.isExactMatch];
else
    keep = [words.isExactMatch] | [words.isCorrectedMatch];
end
words = words(keep);

% Extract counts per cell from stripped words structure
[cellNums,~] = unique([words.cellID],'stable');
numCells = length(cellNums); 
rnaCnts = zeros(140,numCells); 
allGenes = {words.geneName};
for c =1:numCells;
    k = cellNums(c);
    cellGenes = allGenes([words.cellID] == k); 
    [uniqueGenes,geneCnts] = occurrences(cellGenes);
    [~,ia,ib] = intersect(libGenes,uniqueGenes,'stable');
    rnaCnts(ia,c) = geneCnts(ib);
end