function [B,Rg,Attached] = StickyPolymerSim(P,BindPars,varargin)
% % Basic simulation for sticky polymer
% % Binding Chemistry Parameters
% BindPars.TotMols - number of freely diffusing binding molecules
% BindPars.ka -   association binding rate for free molecules
% BindPars.kd -  disassociation rate from site
% BindPars.kbreak - chromatin cluster disassociation rate 
% BindPars.kjoin - chromatin cluster association rate 
% BindPars.EBond - Boltzmann factor associated with each bond
% BindPars.openChromatin;%  = openChromatin;


%% Default Optional Variables
N = size(P,1); 
T = N^2/4; % time for simulation to equilibrate
updateRate = 10; 
showmovie = true;
figHandle = []; 

% Binding Chemistry Parameters
% BindPars.TotMols =  3*N;%  600; % number of freely diffusing Ph molecules
% BindPars.ka  = .03;  % .03 .002 association binding rate for free Ph
% BindPars.kd = .02; % .02 % disassociation rate from site
% BindPars.kbreak  = 0.08; % .01; % chromatin cluster disassociation rate 
% BindPars.kjoin  = .5; % chromatin cluster association rate 
% BindPars.EBond =.8; % .1 % Boltzmann factor associated with each bond

Attached = false(N,N); 
BindPotential = ones(N,1); 

%--------------------------------------------------------------------------
%% Parse variable input
%--------------------------------------------------------------------------
if nargin > 3
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;
    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            case 'T'
                T = CheckParameter(parameterValue,'positive','T');
            case 'updateRate'
                updateRate = CheckParameter(parameterValue,'positive','updateRate');
            case 'showmovie'
                showmovie = CheckParameter(parameterValue,'boolean','showmovie');
            case 'figHandle'
                figHandle = CheckParameter(parameterValue,'handle','figHandle');
            case 'BindPotential'
                BindPotential= CheckParameter(parameterValue,'nonnegative','BindPotential');
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

%% Parse Input
%-------------------------------------------------

% Binding Chemistry Parameters
% TotMols = BindPars.TotMols;%  =  3*N;%  600; % number of freely diffusing Ph molecules
% ka = BindPars.ka;%  = .03;  % .03 .002 association binding rate for free Ph
% kd = BindPars.kd;%  = .02; % .02 % disassociation rate from site
kbreak = BindPars.kbreak;%  = 0.08; % .01; % chromatin cluster disassociation rate 
kjoin = BindPars.kjoin;%  = .5; % chromatin cluster association rate 
EBond = BindPars.EBond;% =.8; % .1 % Boltzmann factor associated with each bond



%% Initialize arrays and ICs
%-------------------------------------
M  = 2*N+2; % linear dimension of lattice (Must be EVEN)

% Initialize some useful arrays for computation
DistMatrix1 = zeros(N,N);
DistMatrix2 = zeros(N,N);
DistMatrix2(triu(true(N),1)) = 100; % 
O = OctehadralGroup; % Precomputed Octahedral Rotation groups for speed
nMoves = 0;  % counter for the number of moves
BP = BindPotential*BindPotential';

% Data Storage Structure
Rg = zeros(T,1); % vector to record the Radius of Gyration. 

% Initial Conditions for simulation
B = P*2 - .5; 
occupiedVertices = GetVertices(B,M); 

%% Begin Simulation
for t = 1:T      
     Rg(t) = RadiusOfGyration(B); % Compute Radius of Gyration 
    c = randi([2,N]); % Choose a Position at Random
    idModel=1; % Use Pivot Method for first attempt to compute a move 

   %% The Pivot Method:     
     if idModel==1   
        % Propose a Pivot move around link c. 
        B2 = PivotMove(B,c,O);

        % Determine if self-avoidance is satisfied 
        newVertices = GetVertices(B2((c+1):N,:),M);
        newVerts = [occupiedVertices(1:c,:); newVertices];
        collide = length(unique(newVerts(:)))==8*N;    

        % Determine number of bonds which break to accomodate move
        DistMatrix1( tril(true(N),-1) ) = pdist(B2);
        BrokenPairs = DistMatrix1 > sqrt(12); % & DistMatrix1~=0; 
        BrokenPairs = Attached & BrokenPairs;
        BrokenBonds = BP(BrokenPairs); % number of bonds broken in each pair
        totalBroken = sum(BrokenBonds(:));

        moveCond=collide && rand < EBond^totalBroken;

        % specific initial move
        if moveCond
            nMoves = nMoves + 1;  % Record that we met the criteria to make a move          
            B = B2;    % Move the particle in the position list
            Attached = Attached & ~BrokenPairs & ~BrokenPairs';
        end     
     end
     if ~moveCond
         idModel=2;
     end

     %% The Bond Fluctuation Method
     if idModel==2
        % Propose a move based on the Bond Flucutation Method
         [newlink,validMove] = BFMmove(B,c);

        % Determine if self avoidance is satisfied
        newVertex = GetVertices(newlink,M);
        collide = sum(ismember(newVertex(:),occupiedVertices(:)))>1;

        % Determine number of bonds broken to accomodate the move
        neigh=find(Attached(c,:));
        bondlength2=zeros(length(neigh),1);
        bondstrength=zeros(length(neigh),1);
        for i=1:length(neigh)
            bondlength2(i) =  ( newlink(1,1) - B(neigh(i),1) )^2 + ...
                              ( newlink(1,2) - B(neigh(i),2) )^2 + ...
                              ( newlink(1,3) - B(neigh(i),3) )^2;
            bondstrength(i) = min(BindPotential(c),BindPotential(neigh(i))); 
        end
        totalBroken=sum((bondlength2>12).*bondstrength);

        % Update Molecule Positions if move is legal; 
        moveCond = validMove && ~collide && rand<EBond^totalBroken;
        if moveCond     
            nMoves = nMoves + 1;  % Record that we met the criteria to make a move          
            B(c,:) = newlink;    % Move the particle in the position list
            for i=1:length(neigh)
                if bondlength2(i)>12
                    Attached(c,neigh(i))=false;
                    Attached(neigh(i),c)=false;
                end
            end
            occupiedVertices(c,:) = newVertex;   
        end
     end
 
    if moveCond           
         %% Chemical interaction between factors bound to polymer:  
  %       Bound = BindMols(ka,kd,Bound,TotMols);

         Attached = FormBonds(B,Attached,kjoin,kbreak,...
                            'DistMatrix',DistMatrix2,...
                            'BindPotential',BindPotential);% ,...
%                                       'OpenSites',BP,

         %% Dynamic Plot of chromatin and Ph clusters with bonds
         if mod(nMoves,updateRate) == 0 && showmovie 
             if ~isempty(figHandle)
                 figure(figHandle); clf;
             else
                 figHandle = figure; clf;
             end
             PlotPolymer(B,'Attached',Attached);
         end                        
                                
    end  % End update molecule positions
    %%
end
