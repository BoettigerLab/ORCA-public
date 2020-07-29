function [gbStruct, rawText] = LoadGenbank(varargin)
%--------------------------------------------------------------------------
% gbStruct = LoadGenbank(fileName, varargin)
% Read a .gb annotation file and create a structure array containing
% information that element
%--------------------------------------------------------------------------
% Necessary Inputs
% fileName/string/(''): Name and path of the annotation to load
%
%--------------------------------------------------------------------------
% Outputs
% gbStruct/structure array: Structure array containing information on each
% entry in the gb file
%
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
% 
% 'exampleInputFlag'/data type/(defalut value): Definition of the flag
%
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% February 3, 2013
%
% Version 1.0
%--------------------------------------------------------------------------
% Creative Commons License CC BY
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Hardcoded Variables
%--------------------------------------------------------------------------
quiet = 0;
flags = {'verbose', 'uniqueEntry', 'newTag', 'prepend'};
%--------------------------------------------------------------------------
% Global Variables
%--------------------------------------------------------------------------
global defaultDataPath;

%--------------------------------------------------------------------------
% Default Variables
%--------------------------------------------------------------------------
verbose = false;
uniqueEntry = true;
fileName = defaultDataPath;
tags = {'ncRNA', 'tmRNA', 'rRNA', 'tRNA', 'gene'};
newTags = {};
prepend = true;
removeJoins = false;

%--------------------------------------------------------------------------
% Parse file information
%--------------------------------------------------------------------------
if nargin < 1 || ismember(varargin{1}, flags)
    [tempFile, tempPath] = uigetfile([fileName '\*.gb']);
    if ~tempFile
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
            case 'uniqueEntry'
                uniqueEntry = CheckParameter(parameterValue,'boolean','uniqueEntry');
            case 'newTag'
                newTags{end+1} = CheckParameter(parameterValue,'string','newTag');
            case 'prepend'
                prepend = CheckParameter(parameterValue,'boolean','prepend');
            case 'removeJoins'
                removeJoins = CheckParameter(parameterValue,'boolean','removeJoins');
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end
%--------------------------------------------------------------------------
% Add newTags
%--------------------------------------------------------------------------
if ~isempty(newTags)
    numNew = length(newTags);
    if prepend
        tags = cat(2, newTags, tags);
    else
        tags = cat(2, tags, newTags);
    end
end
%--------------------------------------------------------------------------
% Open file and read text
%--------------------------------------------------------------------------
fid = fopen(fileName);
if fid < 1
    error(['Error opening file: ' fileName]);
end

rawText = textscan(fid, '%s', 'delimiter', '\n');
fclose(fid);

rawText = rawText{1};
if verbose
    display('----------------------------------------------------------------');
    display(['Loaded ' fileName]);
end

%--------------------------------------------------------------------------
% Find FEATURES Flag
%--------------------------------------------------------------------------
line = 1;
while isempty(strfind(rawText{line}, 'FEATURES'));
    line = line + 1;
end
line = line + 1;

%--------------------------------------------------------------------------
% Parse Features region 
%--------------------------------------------------------------------------
count = 1;
while isempty(strfind(rawText{line}, 'ORIGIN'))
    if strcmp(rawText{line}(1), '/') || ~isempty(strfind(rawText{line}, '"')) || isempty(strfind(rawText{line}, '..'))
        line = line + 1;
    else
        isMatch = false(1, length(tags));
        for i=1:length(tags)
            isMatch(i) = ~isempty(strfind(rawText{line}, tags{i}));
        end
        if any(isMatch)
            % Find feature type
            gbStruct(count).type = strtrim(tags{isMatch});

            % Find feature region and properties
            isComplement = strfind(rawText{line}, 'complement');
            isJoin = strfind(rawText{line}, 'join');

            gbStruct(count).isJoin = ~isempty(isJoin);
            gbStruct(count).isComplement = ~isempty(isComplement);

            openPara = strfind(rawText{line}, '(');
            closedPara = strfind(rawText{line}, ')');
            
            if numel(openPara) > numel(closedPara)
                rawText{line+1} = [rawText{line} rawText{line+1}]; % Handle newline
                line = line + 1;
                openPara = strfind(rawText{line}, '(');
                closedPara = strfind(rawText{line}, ')');
            end
            
            comma = strfind(rawText{line}, ',');
            dots = strfind(rawText{line}, '..');

            tempIndices = [];
            if ~isempty(openPara)
                if isempty(comma)
                    parseStart = [openPara(end) dots+1]+1;
                    parseEnd = [dots closedPara(1)]-1;
                else
                    parseStart = sort([openPara(end) dots+1 comma])+1;
                    parseEnd = sort([dots comma closedPara(1)])-1;
                end
            else % neither complement or join
                parseStart = [15 dots+1]+1;
                parseEnd = [dots length(rawText{line})+1]-1;
            end
            for i=1:length(parseStart)
                indexString = rawText{line};
                indexString = strtrim(indexString([parseStart(i):parseEnd(i)]));
                tempIndices(i) = str2num(['int32(' indexString  ')']);
            end
            gbStruct(count).indices = unique(sort(tempIndices));

            % Find feature name (Always on next line under ' /gene=".." '
            line = line + 1;
            quotes = strfind(rawText{line}, '"');
            indexString = rawText{line};
            gbStruct(count).name = indexString([(quotes(1)+1):(quotes(2)-1)]);
            count = count + 1;
        end
        line = line + 1;
    end
end

if verbose
    display('----------------------------------------------------------------');
    types = unique({gbStruct.type});
    for j=1:length(types)
        display(['Parsed ' num2str(sum(strcmp({gbStruct.type}, types{j}))) ... 
            ' features of type ''' types{j} ''' ']);
    end
end

%--------------------------------------------------------------------------
% Remove duplicates if necessary: precedence set by order of tags
%--------------------------------------------------------------------------
if uniqueEntry
    toKeep = 1:length(gbStruct);
    names = {gbStruct.name};
    uniqueNames = unique(names);
    for i=1:length(uniqueNames)
        ind = find(strcmp(uniqueNames{i}, names));
        if length(ind) > 1
            types = {gbStruct(ind).type};
            toKeep = setdiff(toKeep, ind); % Remove all initially
            for j=1:length(tags)
                ind2 = find(strcmp(types, tags{j}));
                if ~isempty(ind2)
                    toKeep = sort([toKeep ind(ind2(1))]); % Add back the first instance of the highest tag
                    break; % exit loop
                end
            end
        end
    end
    initLen = length(gbStruct);
    
    % Remove degenerate entries
    gbStruct = gbStruct(toKeep);
    
    finalLen = length(gbStruct);

    if verbose
        display(['Removed ' num2str(initLen - finalLen) ' non-unique entries']);
    end
end

%--------------------------------------------------------------------------
% Remove joins
%--------------------------------------------------------------------------
if removeJoins
    isJoin = [gbStruct.isJoin];
    initLen = length(gbStruct);
    
    %Remove joins
    gbStruct = gbStruct(~isJoin);
    finalLen = length(gbStruct);
    
    if verbose
        display(['Removed ' num2str(initLen - finalLen) ' join entries']);
    end
end
    