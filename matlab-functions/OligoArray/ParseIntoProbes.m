function [probeStruct, rawText] = ParseIntoProbes(varargin)
%--------------------------------------------------------------------------
% probeStruct = ParseIntoProbes(fileName, varargin)
% This function reads a flat text file from OligoArray2.0 and loads into a
% matlab structure
%
%--------------------------------------------------------------------------
% Outputs:
%
% probeStruct/structure array: A structure array in which each element
% contains information on a single probe
%
%--------------------------------------------------------------------------
% Inputs:
%
% fileName/string(''): The name and path to a text file from OligoArray
%
%--------------------------------------------------------------------------
% Variable Inputs:
%
% 'verbose'/boolean (true): Display or hide function progress
%  
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% February 6, 2013
%
% Version 1.0
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Hardcoded variables
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Global Variables
%--------------------------------------------------------------------------
global defaultDataPath;

%--------------------------------------------------------------------------
% Define default parameters
%--------------------------------------------------------------------------
verbose = true;
parseName = true;
flags = {'verbose', 'tagName', 'doNotLoad', 'parseName'};
fileName = defaultDataPath;

%--------------------------------------------------------------------------
% Parse file information
%--------------------------------------------------------------------------
if nargin < 1 || ismember(varargin{1}, flags)
    [tempFile, tempPath] = uigetfile([fileName '\*.txt']);
    if tempFile
        params = [];
        display('Canceled load');
    else
        fileName = [tempPath tempFile];
    end
else
    fileName = varargin{1};
    varargin = varargin(2:end);
end

%--------------------------------------------------------------------------
% Parse variable input
%--------------------------------------------------------------------------
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
            case 'parseName'
                parseName = CheckParameter(parameterValue,'boolean','parseName');
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

%--------------------------------------------------------------------------
% Open file and read text
%--------------------------------------------------------------------------
fid = fopen(fileName);
if fid < 0
    error(['Error opening file: ' fileName]);
end

rawText = textscan(fid, '%s', 'delimiter', '\t');
fclose(fid);

rawText = rawText{1};

if verbose
    display(['Opened: ' fileName]);
    display(['Found ' num2str(length(rawText)/9) ' primers']);
end
probeStruct = [];
%--------------------------------------------------------------------------
% Parse each line of text and build structure array
%--------------------------------------------------------------------------
count = 1;
for i=1:9:length(rawText)
    
    % Determine feature name and type
    entryName = rawText{i};
    if parseName
        colon = strfind(entryName, ':');
        dash = strfind(entryName, '-');
    
        probeStruct(count).featureName = entryName( (colon+1):(dash-1) );
        probeStruct(count).featureType = entryName( (dash+1):end );
    else
        probeStruct(count).featureName = rawText{i};
        probeStruct(count).featureType = '';
    end
    
    % Record properties
    probeStruct(count).start = str2num(['int32(' rawText{i+1} ')'] );
    probeStruct(count).length = str2num(['int32(' rawText{i+2} ')'] );
    probeStruct(count).G = str2num(['double(' rawText{i+3} ')'] );
    probeStruct(count).H = str2num(['double(' rawText{i+4} ')'] );
    probeStruct(count).S = str2num(['double(' rawText{i+5} ')'] );
    probeStruct(count).Tm = str2num(['double(' rawText{i+6} ')'] );
    
    props = oligoprop(rawText{i+8});
    probeStruct(count).GC = props.GC;
    
    %Is unique
    probeStruct(count).isUnique = isempty(strfind(rawText{i+7}, ';')  );
    
    %Record probe
    probeStruct(count).seq = rawText{i+8};
    
    %Name the probe
    baseName = probeStruct(count).featureName;
    position = num2str(probeStruct(count).start);
    probeStruct(count).name = [baseName '-' position];
    count = count+1;
    
    if verbose
        if ~mod(floor(i/9), 1000)
            display(['... ' num2str(floor(i/9)) ]);
        end
    end
    
end
if ~isempty(probeStruct)
    if verbose
        display(['Found ' num2str(sum([probeStruct.isUnique])) ' unique probes']);
    end
end
