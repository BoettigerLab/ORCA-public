function [oligoData, parameters] = FoldOligos(sequences, varargin)
%--------------------------------------------------------------------------
% oligoData = FoldOligos(sequences, varargin)
% This function provides an interface to UNAfold and oligoArrayAux folding
% functions. 
% http://mfold.rna.albany.edu/
%--------------------------------------------------------------------------
% Outputs:
%
% oligoData/structure array: The folding information for each of the
% input sequences. This structure has the following fields:
%   -- Header: Header information for each sequence
%   -- Sequence: Sequence for each element of sequences
%   -- dG: The free energy associated with the most stable folded state
%   -- fold: A bracket notation representing the basepairing in the folded
%      state
%
%--------------------------------------------------------------------------
% Inputs:
% 
% sequence/structure array: A structure containing the following fields
%   -- Sequence: The nt sequence of the desired target
%   -- Header: name of the primer 
% This input can also be a path to a fasta file containing the desired
% sequences. 
%--------------------------------------------------------------------------
% Variable Inputs:
% 
% 'verbose'/boolean (true): Display or hide function progress
% 'exePath'/path: Path to executables 
% 'scratchPath'/path: Path to folder where input and output files will be written
% 'nucleicAcidType'/string ('DNA'): Determines if the folded nucleic acid is RNA or DNA
% 'temperature'/integer (37): The temperature at which folding is accomplished 
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% June 3, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Global Variables
%--------------------------------------------------------------------------
global oligoArrayAuxPath;
global oligoArrayAuxScratchPath;

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true};
defaults(end+1,:) = {'exePath', 'fileDir', oligoArrayAuxPath};
defaults(end+1,:) = {'scratchPath', 'fileDir', oligoArrayAuxScratchPath};
defaults(end+1,:) = {'nucleicAcidType', {'DNA', 'RNA'}, 'RNA'};
defaults(end+1,:) = {'tmin', 'nonnegative', 37};
defaults(end+1,:) = {'tmax', 'nonnegative', 37};
defaults(end+1,:) = {'tinc', 'nonnegative', 1};
defaults(end+1,:) = {'dosWindow', 'boolean', false};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabFunctions:invalidArguments', 'A structure array is required.');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

%--------------------------------------------------------------------------
% Create unique ID
%--------------------------------------------------------------------------
import java.util.UUID;
uid = char(UUID.randomUUID());

%--------------------------------------------------------------------------
% Parse input type
%--------------------------------------------------------------------------
if isstruct(sequences)
    if ~all(ismember({'Header', 'Sequence'}, fieldnames(sequences)))
        error('matlabFunctions:invalidArguments', 'Header and Sequence fields must be provided.');
    end
    sequenceFileName = [parameters.scratchPath 'scratchFasta-' uid '.fasta'];
    fastawrite(sequenceFileName, sequences);
    if parameters.verbose
        display('---------------------------------------------------------');
        display(['Wrote sequences: ' sequenceFileName]);
    end
elseif ischar(sequences)
    if exist(sequences) ~= 2 % File
        error('matlabFunctions:invalidArguments', 'Provided file path is not valid.');
    end
    sequenceFileName = sequences;
    sequences = fastaread(sequenceFileName);
    if parameters.verbose
        display('---------------------------------------------------------');
        display(['Read sequences: ' sequenceFileName]);
    end
end

%--------------------------------------------------------------------------
% Define output file
%--------------------------------------------------------------------------
outputFileName = [parameters.scratchPath 'output-' uid];
if parameters.verbose
    display('---------------------------------------------------------');
    display(['Output sequence: ' outputFileName]);
end

%--------------------------------------------------------------------------
% Construct hybrid-ss-min call string
%--------------------------------------------------------------------------
commandString = 'hybrid-ss-min';
commandString = [commandString ' --output=' outputFileName];

% Nucleic Acid
commandString = [commandString ' --NA=' parameters.nucleicAcidType];

%--------------------------------------------------------------------------
% UNDER CONSTRUCTION: Construct hybrid-ss-min call string
%--------------------------------------------------------------------------
commandString = [commandString ' --tmin=' num2str(parameters.tmin)];
commandString = [commandString ' --tmax=' num2str(parameters.tmax)];
commandString = [commandString ' --tinc=' num2str(parameters.tmax)];
%commandString = [commandString ' --sodium=' parameters.sodium];

commandString = [commandString ' ' sequenceFileName];

%--------------------------------------------------------------------------
% Call hybrid-ss-min
%--------------------------------------------------------------------------
if parameters.dosWindow
    commandString = [commandString ' &'];
end
dos(commandString);

%--------------------------------------------------------------------------
% Load and parse energy file
%--------------------------------------------------------------------------
fid = fopen([outputFileName '.dG'], 'r');
data = textscan(fid, '%s');
data = data{1}(6:end);
energy = str2num(char(data(2:3:end)));
fclose(fid);
%--------------------------------------------------------------------------
% Load and parse Ct file
%--------------------------------------------------------------------------
fid = fopen([outputFileName '.ct'], 'r');
brackets = {};
done = false;
while ~feof(fid)
    % Read header line
    headerLine = fgets(fid);
    blockLength = textscan(headerLine, '%d');
    blockLength = blockLength{1};
    basepairs = zeros(2,blockLength);
    for i=1:blockLength
        lineData = textscan(fgets(fid), '%d %s %d %d %d %d %d %d');
        basepairs(1,i) = lineData{1};
        basepairs(2,i) = lineData{5};
    end
    ind = find(basepairs(2,:)~=0);
    tempMat = zeros(blockLength);
    tempMat(sub2ind(size(tempMat), basepairs(1,ind), basepairs(2,ind))) = 1;
    brackets{end+1} = rnaconvert(triu(tempMat));
end
fclose(fid);
%--------------------------------------------------------------------------
% Create output structure
%--------------------------------------------------------------------------
oligoData = sequences;
for i=1:length(energy)
    oligoData(i).dG = energy(i);
    oligoData(i).fold = brackets{i};
end

%--------------------------------------------------------------------------
% UNDER CONSTRUCTION: Cleanup Scratch Files
%--------------------------------------------------------------------------
dirData = dir(parameters.scratchPath);
filesToDelete = find(cellfun(@(x)~isempty(regexp(x, uid)), {dirData.name}));
if parameters.verbose
    display('---------------------------------------------------------');
    display('Deleting temporary files:');
end
for i=1:length(filesToDelete)
    delete([parameters.scratchPath dirData(filesToDelete(i)).name]);
    if parameters.verbose
        display(['   Deleted: ' parameters.scratchPath dirData(filesToDelete(i)).name]);
    end
end
