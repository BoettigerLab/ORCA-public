function [errStruct] = ParseCigarMDTag(MDs, Cigars, varargin)
%--------------------------------------------------------------------------
% [errStruct] = ParseCigarMDTag(MDs, Cigars, varargin)
% This function parses the MD tag from the sam format to determine the
% location and identity of deletions, insertions, and mutations
%--------------------------------------------------------------------------
% Necessary Inputs
% -MDTag/1 or 2D char array: A string or set of strings containing an MD tag. 
% -CigarTag/1 or 2D char array: A string or set of strings contaning a Cigar
%  string
%--------------------------------------------------------------------------
% Outputs
% -errStruct/struct: A structure containing fields that determine the
% location and identity of insertion, deletions, and mismatches
% 
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
% -verbose/boolean (true): Display function progress?
%
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% November 11, 2013
%
% Version 1.0
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA (2013)
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Hardcoded Variables
%--------------------------------------------------------------------------
quiet = 0;
flags = {'verbose'};

%--------------------------------------------------------------------------
% Global Variables
%--------------------------------------------------------------------------
global defaultDataPath;

%--------------------------------------------------------------------------
% Default Variables
%--------------------------------------------------------------------------
verbose = true;
recordPad = true;
numDisplay = 10000;
compact = false; % Determine if the output is an array of structures or a structure of arrays (compact)

%--------------------------------------------------------------------------
% Parse variable input
%--------------------------------------------------------------------------
if nargin < 2
    error('Please provide both an MD and Cigar string');
end

if nargin > 2
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;

    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            case 'verbose'
                verbose = CheckParameter(parameterValue,'boolean',parameterName);
            case 'recordPad'
                recordPad = CheckParameter(parameterValue, 'boolean', parameterName);
            case 'compact'
                compact = CheckParameter(parameterValue, 'boolean', parameterName);
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

%--------------------------------------------------------------------------
% Determine number of input strings
%--------------------------------------------------------------------------
if iscell(MDs)
    MDs = char(MDs);
end
if iscell(Cigars)
    Cigars = char(Cigars);
end
MDDim = size(MDs);
CigarDim = size(Cigars);
if MDDim(1) ~= CigarDim(1)
    error('matlabFunctions::incorrectInput', 'The number of Cigar and MD strings must be the same');
end
numEntries = MDDim(1);

%--------------------------------------------------------------------------
% Initialize Error Structure
%--------------------------------------------------------------------------
if ~compact
    errStruct.padInds = [0 0];
    errStruct.numDel = 0;
    errStruct.numIns = 0;
    errStruct.numMut = 0;
    errStruct.delInds = [];
    errStruct.insInds = [];
    errStruct.mutInds = [];

    errStruct = repmat(errStruct, [1 numEntries]);
else
    errStruct.padInds = zeros(1, 2*numEntries);
    errStruct.numDel = zeros(1, numEntries);
    errStruct.numIns = zeros(1, numEntries);
    errStruct.numMut = zeros(1, numEntries);
    errStruct.delInds = []; % Indices are stored differently here
    errStruct.insInds = [];
    errStruct.mutInds = [];
end

%--------------------------------------------------------------------------
% Parse Strings
%--------------------------------------------------------------------------
for i=1:length(MDs)
    % Index to remove space and empty characters (char(0))
    localMD = MDs(i,MDs(i,:) ~= ' ' & MDs(i,:) ~= char(0));
    localCigar = Cigars(i,Cigars(i,:) ~= ' ' & Cigars(i,:) ~= char(0));
    
    %----------------------------------------------------------------------
    % Parse Cigar to find indels
    %----------------------------------------------------------------------   
    [numStart, numEnd] = regexp(localCigar, '[0-9]+');
    localIndex = 0;
    for j=1:length(numStart)
        number = str2num(localCigar(numStart(j):numEnd(j)));
        
        if j==length(numStart)
            intermediateString = localCigar(numEnd(j)+1:end);
        else
            intermediateString = localCigar(numEnd(j)+1 : numStart(j+1)-1);
        end
        
        switch intermediateString
            case 'M'
                localIndex = localIndex + number;
            case 'D'
                if ~compact
                    errStruct(i).numDel = errStruct(i).numDel + number;
                    errStruct(i).delInds = [errStruct(i).delInds (localIndex+1):(localIndex+number)];
                else
                    errStruct.numDel(i) = errStruct.numDel(i) + number;
                    % Save deletion inds and the index of the sequence
                    localInds = (localIndex+1):(localIndex+number);
                    localIndsAndSequenceID = i*ones(1,2*length(localInds));
                    localIndsAndSequenceID(1:2:end) = localInds;
                    errStruct.delInds = [errStruct.delInds localIndsAndSequenceID];
                end
                localIndex = localIndex + number;
            case 'I'
                if ~compact
                    errStruct(i).numIns = errStruct(i).numIns + 1;
                    errStruct(i).insInds = [errStruct(i).insInds (localIndex+1):(localIndex+number)];
                else
                    errStruct.numIns(i) = errStruct.numIns(i) + 1;
                    % Save insertion inds and the index of the sequence
                    localInds = (localIndex+1):(localIndex+number);
                    localIndsAndSequenceID = i*ones(1,2*length(localInds));
                    localIndsAndSequenceID(1:2:end) = localInds;
                    errStruct.insInds = [errStruct.insInds localIndsAndSequenceID];
                end
                % No need to update local index
            case 'S' % Handle initial skip
                if j==1 && recordPad
                    if ~compact
                        errStruct(i).padInds(1) = number;
                    else
                        errStruct.padInds(2*i-1) = number;
                    end
                elseif j==length(numStart) && recordPad
                    if ~compact
                        errStruct(i).padInds(2) = number;
                    else
                        errStruct.padInds(2*i) = number;
                    end
                end
            otherwise 
        end
    end
    %----------------------------------------------------------------------
    % Parse MD to find mutations
    %----------------------------------------------------------------------  
    [numStart, numEnd] = regexp(localMD, '[0-9]+');
    localIndex = 0;
    for j=1:length(numStart)
        number = str2num(localMD(numStart(j):numEnd(j))); 
        if j==length(numStart)
            intermediateString = localMD(numEnd(j)+1 : end);
        else
            intermediateString = localMD(numEnd(j)+1 : numStart(j+1)-1);
        end
        
        if isempty(regexp(intermediateString, '\^', 'once')) % Skip deletions
            numMut = length(intermediateString);
            localIndex = localIndex + number;
            if ~compact
                errStruct(i).numMut = errStruct(i).numMut + numMut;
                errStruct(i).mutInds = [errStruct(i).mutInds (localIndex + 1) : (localIndex + numMut)];
            else
                errStruct.numMut(i) = errStruct.numMut(i) + numMut;
                % Save mutation inds and the index of the sequence
                localInds = (localIndex + 1) : (localIndex + numMut);
                localIndsAndSequenceID = i*ones(1,2*length(localInds));
                localIndsAndSequenceID(1:2:end) = localInds;
                errStruct.mutInds = [errStruct.mutInds localIndsAndSequenceID];
            end
            localIndex = localIndex + numMut;
        else
            localIndex = localIndex + number + length(intermediateString) - 1; % 
        end
    end
    
    if verbose
        if ~mod(i, numDisplay)
            display(['Completed ' num2str(i) ' sequences']);
        end
    end
    
end