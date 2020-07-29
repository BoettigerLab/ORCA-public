function [badFits,strays] = RemoveBadPts(allMaps,varargin)
% badFits = RemoveBadPts(allMaps,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'maxJump','positive',400};
defaults(end+1,:) = {'maxDist','positive',1800};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'useScalingMatrix','boolean',false};
defaults(end+1,:) = {'scalingExponent','positive',.6};
defaults(end+1,:) = {'scalingMagnitude','positive',300};
pars = ParseVariableArguments(varargin,defaults,mfilename);

%% remove badPts;
[~,numReads,numSpots] = size(allMaps);
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
