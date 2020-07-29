function foundSites = IdentifyBindingSites(library, oligos, varargin)
%--------------------------------------------------------------------------
% foundSites = IdentifyBindingSites(library, oligos, varargin)
% This function identifies binding sites of oligos in members of the
% provided library. 
%
%--------------------------------------------------------------------------
% Outputs:
%
% foundSites/structure array: This array contains information on the
% location and orientation of each of oligos within the library. It
% includes the following fields.
%   --libraryHeader
%   --bindingSite
%   --orientation
%   --oligoHeader
%   --libraryID
%   --oligoID
%--------------------------------------------------------------------------
% Inputs:
%
% library/structure: A structure array containing the following fields
%   -- Sequence: A nt sequence of the library member
%   -- Header: name of the primer
% oligos/structure array: A structure array containing the
% same fields as above but in which the sequences describe oligo. 
%--------------------------------------------------------------------------
% Variable Inputs:
%
% 'verbose'/boolean (true): Display or hide function progress
%  
%--------------------------------------------------------------------------
% Jeffrey Moffitt
% jeffmoffitt@gmail.com
% January 14, 2014
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
if ~strcmp(class(oligos), 'struct')
    error('Oligos must be a structure');
end

if ~all(ismember({'Header', 'Sequence'}, fieldnames(library)))
    error('Library is missing important fields');
end
if ~all(ismember({'Header', 'Sequence'}, fieldnames(oligos)))
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
% Identify sequences that contain oligos
%--------------------------------------------------------------------------
foundSites = [];
count = 1;
for i=1:length(oligos)
    if verbose
        display(['Searching library for ' oligos(i).Header ': ' oligos(i).Sequence]);
    end
    for j=1:length(library)
        forLocations = regexp(library(j).Sequence, oligos(i).Sequence);
        revLocations = regexp(library(j).Sequence, seqrcomplement(oligos(i).Sequence));
        isForward = [true(1, length(forLocations)) false(1, length(revLocations))];
        locations = [forLocations revLocations];
        for k=1:length(isForward)
            foundSites(count).libraryHeader = library(j).Header;
            foundSites(count).oligoHeader = oligos(i).Header;
            foundSites(count).libraryID = j;
            foundSites(count).oligoID = i;
            foundSites(count).bindingSite = [locations(k) locations(k) + length(oligos(i).Sequence)];
            foundSites(count).isForward = isForward(k);
            count = count + 1;
        end
    end
end