function seqStruct = ReadAb1File(varargin)
%--------------------------------------------------------------------------
% seqStruct = ReadAb1File(filename, varargin)
% Read a .ab1 AppliedBiosciences sequencing file and convert it to a
% sequence structure containing the sequence read as well as various
% properties of the sequencing run.
%--------------------------------------------------------------------------
% Necessary Inputs
% fileName/string/(''): Name and path of the annotation to load
%
%--------------------------------------------------------------------------
% Outputs
% seqStruct/structure: A structure containing information on the sequencing
% read provided in the .ab1 file
%
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% August 11, 2013
%
% Version 1.0
%--------------------------------------------------------------------------
% Creative Commons License BY NC SA 2
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Hardcoded Variables
%--------------------------------------------------------------------------
quiet = 0;
flags = {'verbose', 'converterExePath'};

%--------------------------------------------------------------------------
% Global Variables
%--------------------------------------------------------------------------
global defaultDataPath;
global ab1converterPath;

%--------------------------------------------------------------------------
% Default Variables
%--------------------------------------------------------------------------
verbose = false;
fileName = defaultDataPath;
converterExePath = ab1converterPath;

%--------------------------------------------------------------------------
% Parse file information
%--------------------------------------------------------------------------
if nargin < 1 || ismember(varargin{1}, flags)
    [tempFile, tempPath] = uigetfile([fileName '\*.ab1']);
    if ~tempFile
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
            case 'converterExePath'
                converterExePath = CheckParameter(parameterValue, 'filePath', 'converterExePath');
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

%--------------------------------------------------------------------------
% Convert File if Necessary
%--------------------------------------------------------------------------
[~, ~, fileType] = fileparts(fileName);

if ~ismember(fileType, {'.ab1', '.txt'})
    error('File must be a .ab1 or .txt file');
end

if strcmp(fileType, '.ab1')
    
    if exist([fileName '.txt'], 'file')
        fileName = [fileName '.txt'];
        if verbose
            display('----------------------------------------------------------------');
            display(['Using existing converted file: ' fileName ]); 
        end
    else % Convert file
        command = [converterExePath ...
        ' -i ' fileName ...
        ' -t '];

        if verbose
            display('----------------------------------------------------------------');
            display(['Converting ' fileName]);
            display(['    Issued: ' command]);
        end

        failure = dos(command);

        if failure
            error('File conversion failed');
        else
            fileName = [fileName '.txt'];
        end
    end
end
    
%--------------------------------------------------------------------------
% Load text file
%--------------------------------------------------------------------------
fid = fopen(fileName);
if fid < 1
    error(['Error opening file: ' fileName]);
end

if verbose
    display('----------------------------------------------------------------');
    display(['Loading ' fileName]);
end

rawText = textscan(fid, '%s', 'delimiter', '\n'); % Load each line
fclose(fid);
    
rawText = rawText{1}; % Flatten cell

%--------------------------------------------------------------------------
% Record file information in structure
%--------------------------------------------------------------------------
seqStruct.fileName = fileName;
[~, name, fileExt] = fileparts(fileName);
seqStruct.name = [name fileExt];

%--------------------------------------------------------------------------
% Headings To Parse
%--------------------------------------------------------------------------
headingsToParse = {'PBAS', 'DATA', 'PCON'};
entriesToRead = {1, 9:12, 1}; % Determines of all lines with the above tag, which should be read
entryType = {'char', 'num', 'num'};
fieldNames = { {'seq'}, ...
    {'G', 'A', 'T', 'C'} ...
    {'Quality'}};
% data9 = G; data10 = A; data11 = T; data12 = C; 

for i=1:length(headingsToParse)
    entryInds = find(cellfun(@(x)~isempty(regexp(x, headingsToParse{i}, 'once')), rawText));
    
    if verbose
        display(['Reading ' num2str(length(entriesToRead{i})) ' copies of ' headingsToParse{i}]);
    end
    
    for j=1:length(entriesToRead{i})
        rawEntry = rawText{entryInds(j)};
        wsInds = find(isspace(rawEntry));
        switch entryType{i}
            case 'char'
                seqStruct.(fieldNames{i}{j}) = rawEntry( (wsInds(end)+1) : end);
            case 'num'
                temp = textscan(rawEntry( (wsInds(5)+1) : end), '%d', Inf);
                seqStruct.(fieldNames{i}{j}) = temp{1};
        end
    end 
end

%--------------------------------------------------------------------------
% Record sequence length
%--------------------------------------------------------------------------
seqStruct.seqLength = length(seqStruct.seq);


