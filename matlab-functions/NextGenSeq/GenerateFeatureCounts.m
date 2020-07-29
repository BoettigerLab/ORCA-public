function [countStruct] = GenerateFeatureCounts(varargin)
%--------------------------------------------------------------------------
% [countStruct] = GenerateFeatureCounts(coverageVector, featureStructure, ...)
% This function generates a count structure based on the features provided
% in the feature structure
%--------------------------------------------------------------------------
% Necessary Inputs
% coverageVector/structure: A coverage vector structure containing counts
%   at each base on both strands
%
% featureStruct/structure: A structure containing the indices corresponding
%   to different features
%
%--------------------------------------------------------------------------
% Outputs
% countStruct/structure: Structure containing the number of counts and
% other statistics across each feature
%
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
% verbose/boolean(false): Determine whether or not progress is reported
% 
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% February 9, 2013
%
% Version 1.0
%--------------------------------------------------------------------------
% Version 1.1; 
%   JRM; 2/10/13
%   Added feature indices to each element of the count structure 
%--------------------------------------------------------------------------
% Creative Commons License CC BY
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Hardcoded Variables
%--------------------------------------------------------------------------
quiet = 0;
flags = {'verbose'};
requiredFieldsCV = {'refName', 'fileName', 'top', 'bottom'};
requiredFieldsFeat = {'type', 'isJoin', 'isComplement', 'indices', 'name'};

%--------------------------------------------------------------------------
% Global Variables
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Default Variables
%--------------------------------------------------------------------------
verbose = false;

%--------------------------------------------------------------------------
% Parse initial structure information
%--------------------------------------------------------------------------
% Confirm coverageVector structure
if ~isstruct(varargin{1});
    error(['The first argument of ''' mfilename ''' must be a coverageVector structure']);
end
foundFields = fields(varargin{1});
for i=1:length(requiredFieldsCV)
    if ~strcmp(requiredFieldsCV{i}, foundFields)
        error(['A valid coverageVector structure must contain ''' requiredFieldsCV{i} ''' ' ]);
    end
end
cvStruct = varargin{1};
if length(cvStruct) > 1
    error('The coverageVector structure cannot be an array');
end

% Confirm feature structure
if ~isstruct(varargin{2});
    error(['The first argument of ''' mfilename ''' must be a feature structure']);
end
foundFields = fields(varargin{2});
for i=1:length(requiredFieldsFeat)
    if ~strcmp(requiredFieldsFeat{i}, foundFields)
        error(['A valid feature structure must contain ''' requiredFieldsFeat{i} ''' ' ]);
    end
end
featureStruct = varargin{2};

% Update variable arguments in
varargin = varargin(3:end);

%--------------------------------------------------------------------------
% Parse variable input
%--------------------------------------------------------------------------
if (mod(length(varargin), 2) ~= 0 ),
    error(['Extra parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
end
parameterCount = length(varargin)/2;

for parameterIndex = 1:parameterCount,
    parameterName = varargin{parameterIndex*2 - 1};
    parameterValue = varargin{parameterIndex*2};
    switch parameterName
        case 'verbose'
            verbose = CheckParameter(parameterValue,'boolean','verbose');
        otherwise
            error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
    end
end

%--------------------------------------------------------------------------
% Initialize count structure properties
%--------------------------------------------------------------------------
numEntries = length(featureStruct);
uniqueID = now;

if verbose
    display('---------------------------------------------------------------');
    display(['Generating feature counts for ' cvStruct.fileName ]);
    display(['    with reference ' cvStruct.refName]);
end

%--------------------------------------------------------------------------
% Populate count structure
%--------------------------------------------------------------------------
for i=1:numEntries
    % Coverage vector general properties
    countStruct(i).seqName = cvStruct.fileName;
    countStruct(i).refName = cvStruct.refName;
    
    % Feature specific properties
    countStruct(i).featureType = featureStruct(i).type;
    countStruct(i).featureName = featureStruct(i).name;
    countStruct(i).isComplement = featureStruct(i).isComplement;
    countStruct(i).featureIndices = featureStruct(i).indices;
    
    % Unique identifier for this count structure
    countStruct(i).ID = uniqueID;

    % Measure counts for feature
    indices = featureStruct(i).indices;
    topStrand = double(cvStruct.top(indices(1):indices(2)));
    bottomStrand = double(cvStruct.bottom(indices(1):indices(2)));
    countStruct(i).topCount = sum(topStrand);
    countStruct(i).bottomCount = sum(bottomStrand);
    
    countStruct(i).topMean = mean(topStrand);
    countStruct(i).bottomMean = mean(bottomStrand);
    
    countStruct(i).topSTD = std(topStrand);
    countStruct(i).bottomSTD = std(bottomStrand);
    
    countStruct(i).featureLength = length(topStrand);
    
    if featureStruct(i).isComplement
        countStruct(i).counts = sum(bottomStrand);
        countStruct(i).senseMean = mean(bottomStrand);
        countStruct(i).antisenseMean = mean(topStrand);
    else
        countStruct(i).counts = sum(topStrand);
        countStruct(i).senseMean = mean(topStrand);
        countStruct(i).antisenseMean = mean(bottomStrand);
    end

end

if verbose
    display(['Tabulated count statistics for ' num2str(length(countStruct)) ' features']);
end

