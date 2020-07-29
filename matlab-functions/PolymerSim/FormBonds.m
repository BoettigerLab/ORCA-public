function Attached = FormBonds(B,Attached,kjoin,kbreak,varargin)
%--------------------------------------------------------------------------
% Attached = FormBonds(B,kjoin,kbreak);
% Attached = FormBonds(B,kjoin,kbreak,'BindPotential',BP,'DistMatrix',D);
%
%--------------------------------------------------------------------------
% Outputs
% Attached   - NxN matrix 
% 
%--------------------------------------------------------------------------
% Required Inputs
% B, Nx3 matrix,    - specifies the 3D coordinates of all the links in the
%                   polymer chain
% kjoin, scalar,    - specifies the probability [0-1] that two neighboring
%                   bonds will be joined 
% kbreak, scalar    - specifies the probaility [0-1] that two currently
%                   bound links in the chain will break
%
%--------------------------------------------------------------------------
% Optional Inputs
% ('Flag', data-type, default)
% 'BindPotential',Nx1 vector,[] - vector recording the binding potential
%                                   of each unit
% 'OpenSites', Nx1 vector, [],  - vector specifying which sites remain
%                                  open/unbound and will never form bonds
%                                  with other sites. 
% 'DistMatrix', NxN matrix, [] - optional initial Distance Matrix, this
%                                takes some time to generate, so if the
%                                function is to be called in a loop it is
%                                better to create it first and pass to the
%                                function, especially if N is large.  
%
%--------------------------------------------------------------------------
% Alistair Boettiger

%% Fixed Parameters
N = size(B,1); 

%% Default Parameters
DistMatrix = [];
BindPotential = [];
OpenSites = [];
BP = [];
%--------------------------------------------------------------------------
%% Parse variable input
%--------------------------------------------------------------------------
if nargin > 4
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;
    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            case 'DistMatrix' 
                DistMatrix = parameterValue;  % CheckParameter is slow
            case 'BindPotential'
                BindPotential = parameterValue;
            case 'OpenSites'
                OpenSites = parameterValue;
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

if isempty(DistMatrix)
    DistMatrix = zeros(N,N);
    DistMatrix(triu(true(N),1)) = 100; % 
end

if isempty(BindPotential) && ~isempty(OpenSites)    
    [h,w] = size(OpenSites);
    if w>h
        OpenSites = OpenSites';
    end
    BP = OpenSites*OpensSites';
end

%% Main Function
DistMatrix( tril(true(N),-1) ) = pdist(B); % pairwise distances      
DistMatrix( diag(true(N-1,1),-1) ) = 0; 

Neighbors = DistMatrix <= sqrt(12); 
Neighbors(diag(true(N,1))) = 0; 
Neighbors(diag(true(N-1,1),1)) = 0; 
Neighbors(diag(true(N-1,1),-1)) = 0; 
% figure(2); clf; imagesc( Neighbors);  

% Give new neighbors chance to bind and existing bonds a chance to break
NeighborIdx = find(Neighbors);
currLinkIdx = find(Attached); 
if ~isempty(BindPotential)
    [Bi,Bj] = meshgrid(BindPotential);
    BP = min(cat(3,Bi,Bj),[],3);  % A bit slow,
end
if ~isempty(BP)
    newLinks = rand(length(NeighborIdx),1) < kjoin*BP(NeighborIdx);
    brokeLinks = rand(length(currLinkIdx),1) < kbreak./BP(currLinkIdx);
else
    newLinks = rand(length(NeighborIdx),1) < kjoin;
    brokeLinks = rand(length(currLinkIdx),1) < kbreak;
end

% Update attachment matrix 
NewLinks = false(N,N);
BrokeLinks = false(N,N);
NewLinks(NeighborIdx(newLinks)) = true;     
BrokeLinks(currLinkIdx(brokeLinks)) = true;
Attached = (Attached | NewLinks) & ~BrokeLinks;
Attached = Attached | Attached'; % enforce symmetric  binding