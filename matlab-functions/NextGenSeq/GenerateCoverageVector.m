function [cvStruct] = GenerateCoverageVector(varargin)
%--------------------------------------------------------------------------
% [cvStruct] = GenerateCoverageVector(varargin)
% This function generates a strand count structure from a bio-indexed SAM
% file. 
%--------------------------------------------------------------------------
% Necessary Inputs
% SAMfile/bioIndex class or string(''): Either a bioindex class for a sam
% file or the name and path of a sam file
%
%--------------------------------------------------------------------------
% Outputs
% cvStruct/structure array: A structure containing counts per bp for both
% strands of a given reference sequence. Each element in the array
% corresponds to a different reference sequence. 
%
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% February 8, 2013
%
% Version 1.0
%--------------------------------------------------------------------------
% Version 1.1
% February 9, 2013; JRM
% Fixed bug in top/bottom strand identification
%--------------------------------------------------------------------------
% Creative Commons License CC BY
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Hardcoded Variables
%--------------------------------------------------------------------------
quiet = 0;
flags = {'verbose', 'parseChunk'};

%--------------------------------------------------------------------------
% Global Variables
%--------------------------------------------------------------------------
global defaultDataPath;

%--------------------------------------------------------------------------
% Default Variables
%--------------------------------------------------------------------------
fileName = defaultDataPath;


%--------------------------------------------------------------------------
% Parse file information
%--------------------------------------------------------------------------
if nargin < 1
    bioIndex = [];
elseif ~strcmp(class(varargin{1}), 'BioIndexedFile')
    if ismember(varargin{1}, flags)
        bioIndex = [];
    else
        bioIndex = varargin{1};
        varargin = varargin(2:end);
    end
else
    bioIndex = varargin{1};
    varargin = varargin(2:end);
end

if isempty(bioIndex)
    [tempFile, tempPath] = uigetfile([fileName '\*.gb']);
    if ~tempFile
        error('Canceled load');
    else
        bioIndex = [tempPath tempFile];
    end
end

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
    otherwise
        error([mfilename ' does not support objects of class ' class(bioIndex)]);
end

%--------------------------------------------------------------------------
% Parse variable input
%--------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'parseChunk','positive',100000};
defaults(end+1,:) = {'offset','integer',0};

parameters = ParseVariableArguments(varargin,defaults,mfilename);
parseChunk = parameters.parseChunk;
verbose = parameters.verbose;

%--------------------------------------------------------------------------
% Get information on file
%--------------------------------------------------------------------------
inputFile = bioIndex.InputFile;
numSeq = bioIndex.NumEntries;
fileFormat = bioIndex.FileFormat;
if ~strcmp(fileFormat, 'sam');
    error([mfilename ' does not support bioIndex objects from file type ' fileFormat]);
end
info = saminfo(inputFile);
numRefSeq = length(info.SequenceDictionary);
refNames = {info.SequenceDictionary.SequenceName, '*'};
if verbose
    display('-----------------------------------------------------------------');
    display(['Loaded bioindex for ' inputFile ]);
    display(['Found ' num2str(numSeq) ' sequences']);
    display(['Found ' num2str(numRefSeq) ' reference sequences']);
end

%--------------------------------------------------------------------------
% Initialize Coverage Vector Structure
%--------------------------------------------------------------------------
numRepeats = floor(numSeq/parseChunk);

% Initialze strand count structure
uniqueID = now;
for i=1:numRefSeq
    cvStruct(i).refName = info.SequenceDictionary(i).SequenceName;
    cvStruct(i).ID = uniqueID;
    cvStruct(i).fileName = info.Filename;
    cvStruct(i).top = uint32(zeros(1,info.SequenceDictionary(i).SequenceLength));
    cvStruct(i).bottom = uint32(zeros(1,info.SequenceDictionary(i).SequenceLength));
end

%--------------------------------------------------------------------------
% Populate coverage vector
%--------------------------------------------------------------------------
for i=1:numRepeats
    % Read data chunk
    seqData = read(bioIndex, ((i-1)*parseChunk + 1+parameters.offset):(i*parseChunk));
    % Parse sequences
    for j=1:length(seqData)
        % Find reference sequence
        refSeqInd = find(strcmp(refNames, seqData(j).ReferenceName));
        if refSeqInd <= numRefSeq % Check for '*' indicator of no alignment
            % Assign to strand
            if bitget(seqData(j).Flag, 5); %Flag=1 indicates top strand
                cvStruct(refSeqInd).top(seqData(j).Position) = ...
                    cvStruct(refSeqInd).top(seqData(j).Position) + 1;
            else
                cvStruct(refSeqInd).bottom(seqData(j).Position) = ...
                    cvStruct(refSeqInd).bottom(seqData(j).Position) + 1;
            end
        end
    end
    if verbose
        display(['Finished ' num2str(i) ' of ' num2str(numRepeats + 1) ]);
    end
end
% Read last chunk
try
    seqData = read(bioIndex, (numRepeats*parseChunk+1+parameters.offset):numSeq);
catch er
    warning(er.getReport);
    seqData = read(bioIndex, (numRepeats*parseChunk+1+parameters.offset):numSeq-1);
end

for j=1:length(seqData)
    % Find reference sequence
    refSeqInd = find(strcmp(refNames, seqData(j).ReferenceName));
    if refSeqInd <= numRefSeq % Check for '*'
        % Assign to strand
        if bitget(seqData(j).Flag, 5);
            cvStruct(refSeqInd).top(seqData(j).Position) = ...
                cvStruct(refSeqInd).top(seqData(j).Position) + 1;
        else
            cvStruct(refSeqInd).bottom(seqData(j).Position) = ...
                cvStruct(refSeqInd).bottom(seqData(j).Position) + 1;
        end
    end
end


    
    
    
