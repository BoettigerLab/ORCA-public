function B2 = PivotMove(B,varargin)  
%--------------------------------------------------------------------------
% B2 = PivotMove(B)
% B2 = PivotMove(B,c,O)
%--------------------------------------------------------------------------
% Outputs 
% B2 - Nx3 Potentail 3D Polymer Matrix pivoted around c,
%           a randomly chosen link or user specified link. 
%      Note: B2 need not be a self-avoiding polymer.  To test that use
%      vertices = GetVertices and see if the vertices are unique. 
%--------------------------------------------------------------------------
% Inputs 
% B - N x 3 Polymer matrix (N links, 3 dimenions).
% c - link to pivot around
% O - octahedral matrix.  Will be computed de novo if not passed.
%     Precompute this using O=Octehedral group and pass to this function 
%     for speed. 
%--------------------------------------------------------------------------

N = size(B,1);
if nargin == 1
    c = randi([2,N]);
elseif nargin == 2
    c = varargin{1};
    O = OctehadralGroup;
elseif nargin == 3
    c = varargin{1};
    O = varargin{2};
end

%% Main Function
   
pivot=B(c,:);
idT=randi(48); % Choose a transformation from the octohedral group
Matrix=O(:,:,idT);

Btemp = (B(c+1:N,:) - repmat(pivot,N-c,1))*Matrix' + repmat(pivot,N-c,1);
B2 = [B(1:c,:); Btemp];    
    
    