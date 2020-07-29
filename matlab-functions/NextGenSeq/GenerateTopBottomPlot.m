function fig_handle = GenerateTopBottomPlot(varargin)
%--------------------------------------------------------------------------
% fig_handle = GenerateTopBottomPlot(varargin)
% This function generates a top strand bottom strand abundance plot
%--------------------------------------------------------------------------
% Necessary Inputs
% coverageVectorStructure: A structure containing counts for top and bottom
% strand of a set of sequences. Multiple reference sequences can be
% includes as an array of coverage vector structures. 
%
%--------------------------------------------------------------------------
% Outputs
% fig_handles/cell array: A cell array of handles to the various figures
% created by this function
%
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% February 9, 2013
%
% Version 1.0
%--------------------------------------------------------------------------
% Creative Commons License CC BY
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Hardcoded Variables
%--------------------------------------------------------------------------
quiet = 0;
requiredFields = {'refName', 'fileName', 'top', 'bottom'};
%--------------------------------------------------------------------------
% Global Variables
%--------------------------------------------------------------------------
global defaultDataPath;

%--------------------------------------------------------------------------
% Default Variables
%--------------------------------------------------------------------------
verbose = true;
YMode = 'log';
tsColStyle = 'g';
bsColStyle = 'r';
logBase = 2;
continuous = true;
readLength = 50;
countAllBases = false;
titleString = '';
%--------------------------------------------------------------------------
% Confirm cvStruct input
%--------------------------------------------------------------------------
if ~strcmp(class(varargin{1}), 'struct')
    error(['The first argument to ''' mfilename ''' must be a cvStruct']);
end
foundFields = fields(varargin{1});
for i=1:length(requiredFields)
    if ~strcmp(requiredFields{i}, foundFields)
        error(['A valid cvStruct must contain ''' requiredFields{i} ''' ' ]);
    end
end
cvStruct = varargin{1};
varargin = varargin(2:end);

%--------------------------------------------------------------------------
% Parse variable input
%--------------------------------------------------------------------------
presentationPlotArgIn = {};
if nargin > 1
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
            case 'YMode'
                YMode = CheckList(parameterValue, {'log', 'linear'}, 'YMode');
            case 'countAllBases'
                countAllBases = CheckParameter(parameterValue,'boolean','countAllBases');
            case 'titleString'
                titleString = CheckParameter(parameterValue,'string','titleString');
            otherwise
                presentationPlotArgIn(end+1:end+2) = varargin( (parameterIndex*2 - 1):(parameterIndex*2) );
                %error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

