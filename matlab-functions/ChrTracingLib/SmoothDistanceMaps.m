function [newMap,skip] = SmoothDistanceMaps(distMap,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'movingAveBoxWidth','integer',5};
defaults(end+1,:) = {'maxBlankFrac','fraction',.75};
defaults(end+1,:) = {'verbose','boolean',true};
pars = ParseVariableArguments(varargin,defaults,mfilename);


%% moving box instead of moving 
% this should be its own function
% and it should be rewitten for speed with less looping 
%   (or maybe par for at least for large distance maps?) 

if pars.verbose
    disp('filtering distance maps...');
end

[nReads,~,nSpots] = size(distMap);
newMap = nan(nReads,nReads,nSpots);
stp = pars.movingAveBoxWidth;
tic
if stp > 0 
for g=1:nSpots
    skip = false;
    currMap = distMap(:,:,g);
    for a=1:nReads
        for b=1:nReads
            if ~skip
                c1 = max(1,a-stp);
                c2 = min(nReads,a+stp);
                r1 = max(1,b-stp);
                r2 = min(nReads,b+stp);
                cutMap = currMap(c1:c2,r1:r2);
                fracBlank = sum(isnan(cutMap(:)))/(size(cutMap,1)*size(cutMap,2));
                if fracBlank < pars.maxBlankFrac
                    score = nanmean(cutMap(:));
                else
                    score = nan;
                    newMap(:,:,g) = nan(nReads,nReads);
                    skip = true;
                end
                newMap(a,b,g) = score;
            end
        end
    end
end
t = toc;
if pars.verbose
    disp(['filtering complete in ',num2str(t/60),' min']);
end
else
   newMap = distMap;  
end

% remove entries that are purely nan
skip = isnan(sum(reshape(newMap,nReads^2,nSpots)));
