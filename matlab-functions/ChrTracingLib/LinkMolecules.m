function [chainPts,chainQuality,chainIdx] = LinkMolecules(xy,f,varargin)
% chainPts = LinkMolecules(xy,f,varargin)
% link molecules 
% -------------------------------------------------------------------------
%%  Inputs
% -------------------------------------------------------------------------
%   xy is a Nx2 or Nx3 array of molecule positions 
%   f is an Nx1 array of the frames in which these positions were observed.
%   values range from 1 to F total frames (F may larger or smaller than N). 
% 
% -------------------------------------------------------------------------
%% Outputs
% -------------------------------------------------------------------------
%   chainPts
%   Function returns a cell array of chains (chainPts) in which each chain
%   contains only 1 point per frame.  Each cell is a Fx2 or Fx3 matrix of
%   the xy(z) coordinates of molecules in that chain.  
%   The number of chains and their starting positions are determined by
%   clustering the data based on a cutoff distance 'maxDistFromStart' and
%   keeping only clusters which contain data from at least 'minFracHybes'
%   fraction of the total number of hybes. 
%   This determines the number of total chains.  All molecules in subsequent
%   frames that are within a cut-off distance 'maxDistFromStart' from the
%   starting point are included in the chain.  If there is more than one
%   molecule in the same frame the one closest to the previous point in the
%   chain is chosen. This molecule must also be closer than 'maxJump' from
%   the last point.  
%
%   chainQuality
%   A vector of length F frames 
%   if a choice was made between K points, the chain quality for that frame
%   is K/2. 
% 
%   chainIdx
%   A cell array of vectors of length F, recording the indices in the
%   original N points which are associated with each chain.  
% 
%% Optional Parameters %
% -------------------------------------------------------------------------
% defaults(end+1,:) = {'showPlots', 'boolean', false}; % 
% defaults(end+1,:) = {'maxJump', 'positive', 2}; % in pixels
% defaults(end+1,:) = {'maxDistFromStart', 'positive', 5}; % in pixels
% defaults(end+1,:) = {'startPos', 'array', []}; % in pixels
% defaults(end+1,:) = {'brightness', 'array', []}; % in au, brightness of spot 
% defaults(end+1,:) = {'minFracHybes','fraction',.1};
% -------------------------------------------------------------------------
% Alistair Boettiger (boettiger@stanford.edu)
% February 10, 2017
% Revised algorithm, Feb 15 2017. 
%       Added maxDistFromStart and optional startPos focus.  Fixed 
%       potential bug for missing data in first frame. 
% Revised algorithm March 9 2017.
%       Added tracking of index values relative to original chain. These
%       can be used to grab associated data entries 
% 
% Copyright: CC BY NC
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'showPlots', 'boolean', false}; % 
defaults(end+1,:) = {'maxJump', 'positive', 2}; % in pixels
defaults(end+1,:) = {'maxDistFromStart', 'positive', 5}; % in pixels
defaults(end+1,:) = {'startPos', 'array', []}; % in pixels
defaults(end+1,:) = {'brightness', 'freeType', []}; % in au, brightness of spot 
defaults(end+1,:) = {'minFracHybes','fraction',.1};


% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);
% parameters = ParseVariableArguments([], defaults, mfilename);

% -------------------------------------------------------------------------
% Begin main function
% -------------------------------------------------------------------------
numHybes = max(f);
dim = size(xy,2);
if dim == 2 
    is2D = true;
elseif dim == 3
    is2D = false;
else
    error('xy must be a Nx2 or Nx3 position matrix');
end
    
frLabel = cellstr(num2str(f)) ;


if parameters.showPlots
    if is2D
        plot(xy(:,1),xy(:,2),'k.')
        text(xy(:,1),xy(:,2),frLabel);
    else
        plot3(xy(:,1),xy(:,2),xy(:,3),'k.')
        text(xy(:,1),xy(:,2),xy(:,3),frLabel);
    end
end

%%
% xy = xyz_cy3; f = hybe_cy3;



% Find start position
searchFromStart = false;
if isempty(parameters.startPos)
    % cluster based on the number of hybes represented within cutoff distance
    m = squareform( pdist(xy));
    numNeighbors = sum(m<parameters.maxDistFromStart);
    numInCluster = occurrences(numNeighbors);
    xy_1 = [];
    for i=1:length(numInCluster)
        inCluster = numNeighbors == numInCluster(i); 
        hybeInCluster = f(inCluster);
        numHybesInCluster = length(unique(hybeInCluster));
        if numHybesInCluster > parameters.minFracHybes*numHybes
            xy_1 =[xy_1; mean(xy(inCluster,:),1)];
        end
    end  
    startFrame = 1;
    if isempty(xy_1)
        searchFromStart = true;
    end
    if searchFromStart
        % chose a point to trace frames from
        startFrame = 1; % let's start at the first frame. 
        xy_1 = xy(f==startFrame,:); % list of starting points for all frames
        % if the first frame is empty, don't give up, try the the next.
        while isempty(xy_1) && startFrame<numHybes
           startFrame = startFrame+1;
           xy_1 = xy(f==f(startFrame),:); % advance a frame
        end
    end
else
    startFrame = 1;
    xy_1 = parameters.startPos;
   
end

numChains = size(xy_1,1);
try
    distFromStart = pdist2(xy_1,xy);
catch
    error('error running pdist2 in LinkMolecules');
end
s = (1:size(xy,1))'; % an unique index value for each point in the orig list
xyfAll = [xy,f,s];

if ~isempty(parameters.brightness)
    xyfAll = [xy,f,s,parameters.brightness];
end


chainPts = cell(numChains,1);
chainQuality = cell(numChains,1);
chainIdx = cell(numChains,1); 
for n=1:numChains
    linkPts = NaN(numHybes,dim); % store the xy(z) chain 
    linkIdx = zeros(numHybes,1); % store the indices for the chain
    quality = ones(numHybes,1);
    xyf = xyfAll(distFromStart(n,:) < parameters.maxDistFromStart,:); % all the points near my start
    % now the trick is just to select the best choice between duplicates
    % and add NaNs for hybes/frames/links missing data
    
    xy_h = xyf(xyf(:,dim+1)==startFrame,1:dim); % xy(z) of all spots in hybeH
    xy_h_i = xyf(xyf(:,dim+1)==startFrame,dim+2); % unique index values of all spots in hybeH
    
    while isempty(xy_h) && startFrame<numHybes
        startFrame = startFrame +1;
        xy_h = xyf(xyf(:,dim+1)==startFrame,1:dim); 
        xy_h_i = xyf(xyf(:,dim+1)==startFrame,dim+2); % unique index values of all spots in hybeH
    end
    [idx,~]= knnsearch(xy_h, xy_1(n,:));
    linkPts(startFrame,:) = xy_h(idx,:);
    linkIdx(startFrame,:) = xy_h_i(idx); 
    lastLink = linkPts(startFrame,:); % initialize coordiantes of last link
    for h=startFrame+1:numHybes
        xy_h = xyf(xyf(:,dim+1)==h,1:dim); % xy(z) all spots in hybe h
        xy_h_i =xyf(xyf(:,dim+1)==h,dim+2);% the unique index values associated with coordinates above
        xy_h_b = xyf(xyf(:,dim+1)==h,dim+3);% the brightness values associated with coordinates above
        if size(xy_h,1) == 1
            dist = pdist([xy_h; lastLink]);
            if dist<parameters.maxJump
                linkPts(h,:) = xy_h; % if it's unique we keep it
                linkIdx(h) = xy_h_i;
                lastLink = xy_h;
            else
                quality(h) = 0;
                continue % if we have none to chose from skip to next hybe    
            end
        elseif size(xy_h,1) == 0
            quality(h) = 0;
            continue % if we have none to chose from skip to next hybe
        else
            % if we have multiple to chose from, we'll keep the one closest
            % to the previous hybe. 
            if isempty(parameters.brightness)
                [idx,dist] = knnsearch(xy_h,lastLink);
            else
                [~,dist] = knnsearch(xy_h,lastLink);
                [~,idx] = max(xy_h_b);
            end                
             quality(h) = 1/size(xy_h,1); % record the presence of an ambiguous choice
             if dist<parameters.maxJump 
                linkPts(h,:) = xy_h(idx,:); % keep the closest point
                linkIdx(h) = xy_h_i(idx); 
                lastLink = xy_h(idx,:);
             else % except if it's too far away, we toss it. 
                 continue
             end
        end
    end
    if parameters.showPlots
        if is2D
            hold on; plot(linkPts(:,1),linkPts(:,2),'r');
        else
            hold on; plot3(linkPts(:,1),linkPts(:,2),linkPts(:,3),'r');
        end
    end
    chainPts{n} = linkPts;
    chainQuality{n} = quality;
    chainIdx{n} = linkIdx;
end

