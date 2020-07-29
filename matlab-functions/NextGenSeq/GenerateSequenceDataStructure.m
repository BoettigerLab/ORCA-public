function [seqStruct, parameters] = GenerateSequenceDataStructure(bioIndex, varargin)
%--------------------------------------------------------------------------
% seqStruct = GenerateSequenceDataStructure(bioIndex, varargin)
% This function generates a sequence data structure from a bioIndex file
%--------------------------------------------------------------------------
% Necessary Inputs
% -SAMfile/bioIndex class, string(''), or seqData structure: Either a 
%   bioindex class for a sam file, the name and path of a sam file, or the 
%   output structure of a read call to a bioIndex class of a sam file. 
%
%--------------------------------------------------------------------------
% Outputs
% -seqStruct/structure: A structure whose fields are flattened versions of
%   the fields in the SAM file
% -seqData/structure: The matlab default representation of the sam file
% 
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
% -verbose/boolean (true): Display function progress?
%
% -field/string: The name of a SAM file field to include. The default is to
%   include all fields. By specifying a specific field, only that field will
%   be included. This input may be called multiple times to include multiple
%   fields. 
%
% -includeReferenceID/boolean(false): This option determines if the output
%   structure includes a field with the ID of the reference sequence in the
%   fasta file to which the sequence was aligned.  In other words, if the
%   ReferenceName includes an initial number, e.g. '001-Reference Name',
%   then this number is included in the ReferenceID field. 

%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% November 11, 2013
% Version 1.0
% Version 2.0: Improved memory usage: jeffmoffitt@gmail.com
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA (2013)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Hardcoded Variables
%--------------------------------------------------------------------------
possibleFields = cell(0,3);
possibleFields(end+1,:) = {'QueryName', @char, 75}; % Field name, type, and dimension
possibleFields(end+1,:) = {'Flag', @uint16, 1};
possibleFields(end+1,:) = {'ReferenceName', @char, 75};
possibleFields(end+1,:) = {'Position', @uint32, 1};
possibleFields(end+1,:) = {'MappingQuality', @uint8, 1};
possibleFields(end+1,:) = {'CigarString', @char, 75};
possibleFields(end+1,:) = {'Sequence', @char, 200};
possibleFields(end+1,:) = {'Quality', @char, 200};
possibleFields(end+1,:) = {'MD', @char, 75};
possibleFields(end+1,:) = {'ReferenceID', @uint16, 1};

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true}; 
defaults(end+1,:) = {'includeReferenceID', 'boolean', true}; % Flag for image mlist
defaults(end+1,:) = {'fieldsToInclude', 'cell', possibleFields(1:(end-1),1)}; % Do not include referenceID by default
defaults(end+1,:) = {'seqChunk', 'nonnegative', []};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1 || ~ismember(class(bioIndex), {'char', 'BioIndexedFile', 'struct'}) % 7=folder
    error('matlabFunctions:invalidArguments', 'bioIndex must be a path, a BioIndexedFile class, or a structure');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

%--------------------------------------------------------------------------
% Include ReferenceID if required
%--------------------------------------------------------------------------
if parameters.includeReferenceID
    parameters.fieldsToInclude{end+1} = 'ReferenceID';
end

% -------------------------------------------------------------------------
% Handle various types of bioIndex
% -------------------------------------------------------------------------
switch class(bioIndex)
    case 'char'
        display('-----------------------------------------------------------------');
        display(['Load bioIndex Class ' bioIndex]);
        tic;
        bioIndex = BioIndexedFile('sam', bioIndex);
        time = toc;
        display(['Finished in ' num2str(time) ' s']);
    case 'BioIndexedFile'
        % Nothing needed here
    case 'struct'
        if ~isfield(bioIndex, 'QueryName')
            error('Supplied structure does not appear to be derived from a bioIndex of a sam file');
        end
    otherwise
        error([mfilename ' does not support objects of class ' class(bioIndex)]);
end

%--------------------------------------------------------------------------
% Get information on provided file
%--------------------------------------------------------------------------
switch class(bioIndex)
    case 'BioIndexedFile'
        inputFile = bioIndex.InputFile;
        numSeq = bioIndex.NumEntries;
        fileFormat = bioIndex.FileFormat;
        if ~strcmp(fileFormat, 'sam');
            error([mfilename ' does not support bioIndex objects from file type ' fileFormat]);
        end
        info = saminfo(inputFile); 
        numRefSeq = length(info.SequenceDictionary);
        refNames = {info.SequenceDictionary.SequenceName, '*'};
        if parameters.verbose
            display('-----------------------------------------------------------------');
            display(['Loaded bioindex for ' inputFile ]);
            display(['Found ' num2str(numSeq) ' sequences']);
            display(['Found ' num2str(numRefSeq) ' reference sequences']);
        end
    case 'struct'
        numSeq = length(bioIndex);
        refNames = unique({bioIndex.ReferenceName});
        numRefSeq = length(refNames);
        if parameters.verbose
            display('-----------------------------------------------------------------');
            display(['Parsed provided structure']);
            display(['Found ' num2str(numSeq) ' sequences']);
            display(['Found ' num2str(numRefSeq) ' reference sequences']);
        end
end

%--------------------------------------------------------------------------
% Determine subLoop indices
%--------------------------------------------------------------------------
if isempty(parameters.seqChunk)
    parameters.seqChunk = numSeq;
end
numChunks = ceil(numSeq/parameters.seqChunk);
indsPerLoop = {};
for i=1:numChunks
    indsPerLoop{i} =( (parameters.seqChunk*(i-1)) + 1 ) : ( ...
        min( (parameters.seqChunk*i), numSeq) );
end

%--------------------------------------------------------------------------
% Preallocate flatDataStructure
%--------------------------------------------------------------------------
for i=1:length(parameters.fieldsToInclude)
    fieldID = find(strcmp(parameters.fieldsToInclude{i}, possibleFields(:,1)));
    convFun = possibleFields{fieldID, 2};
    fieldLen = possibleFields{fieldID, 3};
    if fieldLen > 1 
        seqStruct.(parameters.fieldsToInclude{i}) = repmat(convFun(' '), [numSeq, fieldLen]);
    else
        seqStruct.(parameters.fieldsToInclude{i}) = convFun(zeros(1, numSeq));
    end
end

%--------------------------------------------------------------------------
% Variables for improving speed of referenceID parsing
%--------------------------------------------------------------------------
detectedDashInd = [];

%--------------------------------------------------------------------------
% Read sequencing data and loop over structure
%--------------------------------------------------------------------------
for i=1:length(indsPerLoop)
    %--------------------------------------------------------------------------
    % Read data chunk
    %--------------------------------------------------------------------------
    switch class(bioIndex)
        case 'BioIndexedFile'
            if parameters.verbose
                display(['Reading entries ' num2str(indsPerLoop{i}(1)) '-' num2str(indsPerLoop{i}(end)) ... 
                    ' of ' num2str(numSeq) ' total sequences']);
            end

            tic;
            seqData = read(bioIndex, indsPerLoop{i});
            elapsedTime = toc;

            if parameters.verbose
                display(['Completed in ' num2str(elapsedTime) ' s']);
            end
        case 'struct'
            seqData = bioIndex;
    end
    
    %----------------------------------------------------------------------
    % Parse structure
    %----------------------------------------------------------------------
    for j=1:length(parameters.fieldsToInclude)
        switch parameters.fieldsToInclude{j}
            case {'QueryName', 'ReferenceName', 'CigarString', 'Sequence', 'Quality'}
                fieldID = find(strcmp(parameters.fieldsToInclude{j}, possibleFields(:,1)));
                fieldLen = possibleFields{fieldID, 3};
                tempData = char({seqData.(parameters.fieldsToInclude{j})});
                dim = size(tempData);
                if dim(2) > fieldLen;
                    warning('matlabFuctions:GenerateSequenceDataStructure', ['Allocated storage is insufficient for ' parameters.fieldsToInclude{j}]);
                end
                seqStruct.(parameters.fieldsToInclude{j})(indsPerLoop{i}, 1:min(dim(2), fieldLen)) = tempData(:,1:min(dim(2), fieldLen));
            case {'Flag', 'Position', 'MappingQuality'}
                seqStruct.(parameters.fieldsToInclude{j})(indsPerLoop{i}) = [seqData.(parameters.fieldsToInclude{j})];
            case 'MD'
                for k=1:length(seqData)
                    if isfield(seqData(k).Tags, 'MD');
                        seqStruct.MD(indsPerLoop{i}(k),1:length(seqData(k).Tags.MD)) = seqData(k).Tags.MD;
                    end
                end
            case 'ReferenceID'
                fieldID = find(strcmp(parameters.fieldsToInclude{j}, possibleFields(:,1)));
                fieldConv = possibleFields{fieldID, 2};
                if isempty(detectedDashInd) % find dash once to significantly speed conversion
                    for k=1:length(seqData)
                        dashInds = regexp(seqData(k).ReferenceName, '-', 'once');
                        if ~isempty(dashInds)
                            detectedDashInd = dashInds(1);
                            if parameters.verbose
                                display(['Detected dash at position ' num2str(detectedDashInd)]);
                            end
                            break;
                        end
                        if k==length(seqData)
                            error('matlabFuctions:GenerateSequenceDataStructure', 'No dash ind dectected. Unable to parse referenceID');
                        end
                    end
                end
                seqStruct.ReferenceID(indsPerLoop{i}) = cellfun(@(x) str2numWithZero(x(1:min((detectedDashInd-1), length(x)))), ...
                    {seqData.ReferenceName});
        end
    end
    if parameters.verbose
        display('...Parsed tags')
    end

end

    function temp = str2numWithZero(x)
        temp = str2num(x);
        if isempty(temp)
            temp = 0;
        end
    end
end

    
    
