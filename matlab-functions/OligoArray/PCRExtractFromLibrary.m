function PCRProducts = PCRExtractFromLibrary(library, primers, varargin)
%--------------------------------------------------------------------------
% PCRProducts = PCRExtractFromLibrary(library, primers, varargin)
% This function identifies priming regions within a library file and
% returns a structure array, foundMembers, that contains the name and
% sequence of all library members that would be amplified with the primers
% included in the list.  
%
%--------------------------------------------------------------------------
% Outputs:
%
% foundMembers/structure array: This array contains a structure for each
% PCR product identified from each primer combination amplified from the
% sequences provided in library. It contains the following fields
%   -- Header: header of the input library member from which the sequence
%   is amplified
%   -- Sequence: amplified sequence
%   -- ForwardPrimer: Index of the primer that amplified the sequence in
%   the forward direction.
%   -- ForwardPrimer: Index of the primer that amplified the sequence in
%   the reverse direction.
%
%--------------------------------------------------------------------------
% Inputs:
%
% library/structure: A structure array containing the following fields
%   -- Sequence: A nt sequence of the library member
%   -- Header: name of the primer
% primers/structure array: A structure array containing the following the
% same fields as above but in which the sequences describe primers. 
%--------------------------------------------------------------------------
% Variable Inputs:
%
% 'verbose'/boolean (true): Display or hide function progress
%  
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% December 16, 2013
%
% Version 1.0
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Hardcoded variables
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Global Variables
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Define default parameters
%--------------------------------------------------------------------------
verbose = false;
flags = {'verbose'};

%--------------------------------------------------------------------------
% Parse required inputs
%--------------------------------------------------------------------------
if ~strcmp(class(library), 'struct')
    error('Library must be a structure');
end
if ~strcmp(class(primers), 'struct')
    error('Primers must be a structure');
end

if ~all(ismember({'Header', 'Sequence'}, fieldnames(library)))
    error('Library is missing important fields');
end
if ~all(ismember({'Header', 'Sequence'}, fieldnames(primers)))
    error('Primers is missing important fields');
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
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

%--------------------------------------------------------------------------
% Identify sequences that contain primers
%--------------------------------------------------------------------------
containsForwardPrimer = false(length(primers), length(library));
containsReversePrimer = false(length(primers), length(library));

for i=1:length(primers)
    containsForwardPrimer(i,:) = cellfun(@(x) ~isempty(regexp(x, primers(i).Sequence)), {library.Sequence});
    containsReversePrimer(i,:) = cellfun(@(x) ~isempty(regexp(x, seqrcomplement(primers(i).Sequence))), {library.Sequence});
end

if length(primers) > 1 % Check for single primer input
    willBeAmplified = any(containsForwardPrimer) & any(containsReversePrimer);
else
    willBeAmplified = containsForwardPrimer & containsReversePrimer;
end

if verbose
    display('---------------------------------------------------------------------');
    display(['Found ' num2str(sum(willBeAmplified)) ' of ' num2str(length(library)) ' library members that will be amplified']);
end

%--------------------------------------------------------------------------
% Extract amplified sequences
%--------------------------------------------------------------------------
amplifiedInds = find(willBeAmplified);
PCRProducts = [];

count = 1; 
for i=1:length(amplifiedInds)
    localLibraryInd = amplifiedInds(i);
    forwardPrimerInds = find(containsForwardPrimer(:, localLibraryInd));
    reversePrimerInds = find(containsReversePrimer(:, localLibraryInd));

    % Loop overall possible configurations
    for j=1:length(forwardPrimerInds)
        for k=1:length(reversePrimerInds)
            forwardLocation = regexp(library(localLibraryInd).Sequence, primers(forwardPrimerInds(j)).Sequence);
            
            reverseLocation = regexp(library(localLibraryInd).Sequence, seqrcomplement(primers(reversePrimerInds(k)).Sequence));
            reverseLocation = reverseLocation + length(primers(reversePrimerInds(k)).Sequence) - 1; % Add offset to start location
            
            if reverseLocation > forwardLocation % Viable amplification
                PCRProducts(count).Header = library(localLibraryInd).Header;
                PCRProducts(count).Sequence = library(localLibraryInd).Sequence(forwardLocation:reverseLocation);
                PCRProducts(count).ForwardPrimer = forwardPrimerInds(j);
                PCRProducts(count).ReversePrimer = reversePrimerInds(k);
                count = count + 1;
            end
        end
    end
end


