function [badFits,strays,strayID] = RemoveBadPts(allMaps,varargin)
% badFits = RemoveBadPts(allMaps,varargin)
% badFits is a logical array, N-monomers x N-monomers x N-observations
% strays is a logical array, N-monomers x Nobservations 
% 
% Input can be either a single map, and stack of maps, a single polymer or
% a stack of polymers. 
%
% Examples:
% [badFits,strays] = RemoveBadPts(allMaps)
% Remove bad points from the map: allMaps(badFits) = NaN
% Reomve bad points from the polymer list: polymerList() = NaN;

defaults = cell(0,3);
defaults(end+1,:) = {'maxJump','positive',400};
defaults(end+1,:) = {'maxDist','positive',1800};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'useScalingMatrix','boolean',false};
defaults(end+1,:) = {'scalingExponent','positive',.6};
defaults(end+1,:) = {'scalingMagnitude','positive',300};
pars = ParseVariableArguments(varargin,defaults,mfilename);

%% remove badPts;
[numReads,dim2,numSpots] = size(allMaps);

% if polymers were passed instead of maps, build maps
if dim2 < numReads && dim2 < 5 
    numSpots = size(allMaps,3);
    polydata = allMaps;
    allMaps = nan(numReads,numReads,numSpots);
    for n=1:numSpots
        allMaps(:,:,n) = squareform(pdist(polydata(:,1:3,n)));
    end
end

% otherwise we just continue with the maps
if pars.verbose
    disp(['Showing data from ',num2str(numSpots) ' cells']);
end

strays = false(numReads,numSpots);
if ~pars.useScalingMatrix
    badFits = false(numReads,numReads,numSpots);
    for s=1:numSpots
        badPts = false(numReads);
        distMap_s = allMaps(:,:,s);

        % remove points the jump out AND jump back more than "maxJump" in 1 step
        jumps = diag(distMap_s,1) > pars.maxJump;
        % stray = find([0; diff(jumps)==0] & jumps==1);
        stray = ([0; diff(jumps)==0] & jumps==1);
        badPts(stray,:) = true;
        badPts(:,stray) = true;
        strays(stray,s) = true;

        % remove points whose median distance from all others exceeds "maxDist" 
        % stray = find(nanmedian(distMap_s) > pars.maxDist);
        stray = (nanmedian(distMap_s) > pars.maxDist);  % boolean indexing is faster
        badPts(stray,:) = true;
        badPts(:,stray) = true;
        strays(stray,s) = true;

        % record 
        badFits(:,:,s) = badPts;
    end
else

    maxJumpMatrix = zeros(numReads,numReads); 
    for j=1:numReads
        d = pars.scalingMagnitude*j^pars.scalingExponent;
        maxJumpMatrix = maxJumpMatrix + d*diag(ones(numReads-j,1),j)+ d*diag(ones(numReads-j,1),-j);
    end
    % figure(1); clf; imagesc(maxJumpMatrix); colorbar; colormap(jet);

    badFits = false(numReads,numReads,numSpots);
    for s=1:numSpots
        distMap_s = allMaps(:,:,s);   
        badPts = distMap_s > maxJumpMatrix;
        badFits(:,:,s) = badPts;
    end

end

if nargout>2
     strayID = permute(repmat(strays,1,1,4),[1,3,2]);
end
