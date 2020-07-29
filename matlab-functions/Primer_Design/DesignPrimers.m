function [primerPairs, leftPrimers, rightPrimers] = DesignPrimers(sequence, varargin)
%--------------------------------------------------------------------------
% primers = DesignPrimers(sequence, varargin)
% This function produces a structure array of primers derived from the
% input sequence and which fit the specified design constraints
%
% Primer design is accomplished with Primer3: http://primer3.sourceforge.net/
%--------------------------------------------------------------------------
% Outputs:
%
% primers/structure array: This array contains a structure for each
% identified primer in sequence. It contains the following fields  
%   -- Header: Name of the primer derived from the sequence name
%   -- Sequence: primer sequence
%   -- isReverse: a boolean that determines the orientation of the primer
%
%--------------------------------------------------------------------------
% Inputs:
% 
% sequence/structure: A structure containing the following fields
%   -- Sequence: The nt sequence of the desired target
%   -- Header: name of the primer 
%--------------------------------------------------------------------------
% Variable Inputs:
% 
% 'verbose'/boolean (true): Display or hide function progress
% 'primer3Path'/path: Path to primer3 folder. 
% 'outputPath'/path: Path to folder where input and output files will be written
% 'primerTask'/string ('pick_sequencing_primers'): Controls function of
%   primer3
% 'targetTm'/array ([0 Inf]): The desired target Tm for the primer 
%--------------------------------------------------------------------------
% Documentation on Primer3 is found at:
% http://sourceforge.net/projects/primer3/files/primer3/2.3.6/primer3_manual.htm/download
%    Untergasser A, Cutcutache I, Koressaar T, Ye J, Faircloth BC, Remm M and Rozen SG.
%    Primer3--new capabilities and interfaces.
%    Nucleic Acids Res. 2012 Aug 1;40(15):e115. 
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% December 27, 2013
%
% Version 1.0
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Hardcoded variables
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Global Variables
%--------------------------------------------------------------------------
global primer3Path;
global primer3ScratchPath;

%--------------------------------------------------------------------------
% Define default parameters
%--------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true};
defaults(end+1,:) = {'veryVerbose', 'boolean', false};
defaults(end+1,:) = {'primer3Path', 'path', primer3Path};
defaults(end+1,:) = {'outputPath', 'path', primer3ScratchPath};
defaults(end+1,:) = {'primerTask', {'generic'}, 'generic'};
defaults(end+1,:) = {'optimalLength', 'positive', 20};
defaults(end+1,:) = {'lengthRange', 'array', [15 21]};
defaults(end+1,:) = {'numPrimers', 'positive', 1};
defaults(end+1,:) = {'productLengthRange', 'array', []};
defaults(end+1,:) = {'sequenceTarget', 'array', []};
defaults(end+1,:) = {'keepIOFiles', 'boolean', false};

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

%--------------------------------------------------------------------------
% Parse required inputs
%--------------------------------------------------------------------------
if ~strcmp(class(sequence), 'struct')
    error('Sequence must be a struct');
end
    
if ~all(ismember(fieldnames(sequence), {'Header', 'Sequence'}))
    error('Library is missing important fields');
end

%--------------------------------------------------------------------------
% Confirm primer3 location
%--------------------------------------------------------------------------
if isempty(parameters.primer3Path)
    error(['No defined location for primer3: download (http://primer3.sourceforge.net/) ' ...
        ' and define path to exe']);
end

%--------------------------------------------------------------------------
% Define File ID and output file
%--------------------------------------------------------------------------
fileID = round((now*1e2 - floor(now*1e2))*10000);
outputFileName = ['primer3_output_' num2str(fileID) '.txt'];
if parameters.verbose
    display('---------------------------------------------------------------');
    display(['I/O File Timestamp: ' num2str(fileID)]);
end

%--------------------------------------------------------------------------
% Create Primer3 Input File
%--------------------------------------------------------------------------
inputFileName = ['primer3_input_' num2str(fileID) '.txt'];
if exist([parameters.outputPath inputFileName]);
    delete([parameters.outputPath inputFileName]);
end

fid = fopen([parameters.outputPath inputFileName], 'w');
% Write sequence name and sequence
fprintf(fid, '%s\n', ['SEQUENCE_ID=' sequence.Header(~isspace(sequence.Header))]);
fprintf(fid, '%s\n', ['SEQUENCE_TEMPLATE=' sequence.Sequence(~isspace(sequence.Sequence))]);

% Write primer task
fprintf(fid, '%s\n', ['PRIMER_TASK=' parameters.primerTask]);
fprintf(fid, '%s\n', ['PRIMER_PICK_LEFT_PRIMER=1']);
fprintf(fid, '%s\n', ['PRIMER_PICK_INTERNAL_OLIGO=0']);
fprintf(fid, '%s\n', ['PRIMER_PICK_RIGHT_PRIMER=1']);

% Write primer targets
fprintf(fid, '%s\n', ['PRIMER_OPT_SIZE=' num2str(parameters.optimalLength)]);
fprintf(fid, '%s\n', ['PRIMER_MIN_SIZE=' num2str(parameters.lengthRange(1))]);
fprintf(fid, '%s\n', ['PRIMER_MAX_SIZE=' num2str(parameters.lengthRange(2))]);

% Write product properties
if ~isempty(parameters.productLengthRange)
    fprintf(fid, '%s\n', ['PRIMER_PRODUCT_SIZE_RANGE=' ...
        num2str(parameters.productLengthRange(1)) '-' num2str(parameters.productLengthRange(2))]);
end

if ~isempty(parameters.sequenceTarget)
    fprintf(fid, '%s\n', ['SEQUENCE_TARGET=' ...
        num2str(parameters.sequenceTarget(1)) ',' num2str(parameters.sequenceTarget(2))]);
end

% Determine Output File Format
fprintf(fid, '%s\n', 'P3_FILE_FLAG=0');
fprintf(fid, '%s\n', 'PRIMER_EXPLAIN_FLAG=0');

% End file entry
fprintf(fid, '%s\n', '=');
fclose(fid);

%--------------------------------------------------------------------------
% Call Primer3
%--------------------------------------------------------------------------
current_dir = pwd; 
cd(parameters.primer3Path);
primer3Command =    ['primer3_core.exe' ...
                     ' -output ' parameters.outputPath outputFileName];
dos([primer3Command ' ' parameters.outputPath inputFileName]);
cd(current_dir);

%--------------------------------------------------------------------------
% Load Primer3 Output
%--------------------------------------------------------------------------
fid = fopen([parameters.outputPath outputFileName]);
rawLines = textscan(fid, '%s');
rawLines = rawLines{1};
fclose(fid);

%--------------------------------------------------------------------------
% Parse Primer3 Output: Find only PRIMER_LEFT, PRIMER_RIGHT, PRIMER_PAIR
%--------------------------------------------------------------------------
leftPrimers = [];
rightPrimers = [];
primerPairs = [];
internal = [];

if parameters.veryVerbose
    display('---------------------------------------------------------------');
    display('Primer3 output: ');
end

for i=1:length(rawLines)
    lineText = rawLines{i};
    if parameters.veryVerbose
        display(['    ' lineText]);
    end
    
    % Split string into flag and value
    equalsInd = regexp(lineText, '=');
    if ~isempty(equalsInd)
        flag = lineText(1:(equalsInd(1)-1));
        value = lineText((equalsInd(1)+1):end);
    else
        flag = [];
        value = [];
    end
    
    % Handle initial flag cases
    switch flag
        case 'PRIMER_LEFT_NUM_RETURNED'
            numLeftPrimers = str2num(value);
        case 'PRIMER_RIGHT_NUM_RETURNED'
            numRightPrimers = str2num(value);
        case 'PRIMER_INTERNAL_NUM_RETURNED'
            numInternalPrimers = str2num(value);
        case 'PRIMER_PAIR_NUM_RETURNED'
            numPairs = str2num(value);
    end
        
    isSpecificPrimer = ~isempty(regexp(flag, '_[0-9]'));
    usInds = regexp(flag, '_');
    if isSpecificPrimer
        if length(usInds) > 2
            primerID = str2num(flag( (usInds(2)+1):(usInds(3)-1)));
            primerID = primerID + 1; % Take first element from 0 to 1
            elementName = flag( (usInds(3)+1):end);
            primerType = flag( (usInds(1)+1):usInds(2)-1);
            
            switch elementName
                case {'PENALTY', 'TM', 'GC_PERCENT', ...
                        'SELF_ANY_TH', 'SELF_END_TH', ...
                        'HAIRPIN_TH', 'END_STABILITY', ...
                        'COMPLY_ANY_TH', 'COMPLY_END_TH', ...
                        'PRODUCT_SIZE'}
                    value = str2num(value);
                otherwise
                    value = value; % Leave as string
            end
            elementName = lower(elementName); % To match typical naming conventions
            switch primerType
                case 'LEFT'
                    leftPrimers(primerID).(elementName) = value;
                case 'RIGHT'
                    rightPrimers(primerID).(elementName) = value;
                case 'PAIR'
                    primerPairs(primerID).(elementName) = value;
            end
        end
    end
end

%--------------------------------------------------------------------------
% Report status
%--------------------------------------------------------------------------
if parameters.verbose || parameters.veryVerbose
    display('---------------------------------------------------------------');
    display(['Found ' num2str(numLeftPrimers) ' left primers']);
    display(['Found ' num2str(numRightPrimers) ' right primers']);
    display(['Found ' num2str(numInternalPrimers) ' internal primers']);
    display(['Found ' num2str(numPairs) ' primer pairs']);
end

%--------------------------------------------------------------------------
% Compile primer properties into primer pairs
%--------------------------------------------------------------------------
for i=1:length(primerPairs)
    primerPairs(i).leftPrimer = leftPrimers(i);
    primerPairs(i).rightPrimer = rightPrimers(i);
end

%--------------------------------------------------------------------------
% Cleanup Scratch Files
%--------------------------------------------------------------------------
if ~parameters.keepIOFiles
    delete([parameters.outputPath outputFileName]);
    delete([parameters.outputPath inputFileName]);
    if parameters.veryVerbose
        display('---------------------------------------------------------------');
        display(['Deleted: ' parameters.outputPath inputFileName]);
        display(['Deleted: ' parameters.outputPath outputFileName]);
    end
end


