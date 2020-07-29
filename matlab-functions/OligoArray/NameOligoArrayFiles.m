function designStruct = NameOligoArrayFiles(designStruct, filePath, varargin)
%--------------------------------------------------------------------------
% designStruct = NameOligoArrayFiles(designStruct, filePath)
% This function populates the various file path fields in a design
% structure with a name derived from the various elements
%
%--------------------------------------------------------------------------
% Outputs:
%
% designStruct/struct: A structure containing fields with the required
% parameters for OligoArray2.0
%
%--------------------------------------------------------------------------
% Inputs:
% designStruct/struct: A structure containing fields with the required
% parameters for OligoArray2.0 populated. 
% 
% filePath/path: A path to the final locations of the desired files
%--------------------------------------------------------------------------
% Variable Inputs:
% 'baseName'/string('oligos_'): A string containing the base name for all
% files. 
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% November 21, 2013
%
% Version 1.0
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Define default parameters
%--------------------------------------------------------------------------
baseName = 'oligos';

%--------------------------------------------------------------------------
% Parse variable input
%--------------------------------------------------------------------------
if nargin > 2
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;

    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            case 'baseName'
                baseName = CheckParameter(parameterValue,'string','baseName');
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

%--------------------------------------------------------------------------
% Generate Default String
%--------------------------------------------------------------------------
maxOligosString = ['_mN_' num2str(designStruct.maxOligos)];
lengthString = ['_len_' num2str(designStruct.minLength) ...
    '_' num2str(designStruct.maxLength)];
TmString = ['_Tm_' num2str(designStruct.minTm) ...
    '_' num2str(designStruct.maxTm)];
ssString = ['_SS_' num2str(designStruct.secStructTm)];
crossHybString = ['_CH_' num2str(designStruct.crossHybTm)];
GCString = ['_GC_' num2str(designStruct.minGC) ...
    '_' num2str(designStruct.maxGC)];
minSepString = ['_sep_' num2str(designStruct.minSep)];

defaultNameString = [baseName maxOligosString lengthString ...
    TmString GCString ssString crossHybString minSepString];

%--------------------------------------------------------------------------
% Populate Names
%--------------------------------------------------------------------------
designStruct.outputFile = [filePath defaultNameString '.txt'];
designStruct.failedFile = [filePath defaultNameString '_failed.txt'];
designStruct.logFile = [filePath defaultNameString '_log.txt'];
