function [boundaryMap,jumpPoints] = InsulationWindow(distMapFilt,varargin)


defaults = cell(0,3);
defaults(end+1,:) = {'w','positive',10};
defaults(end+1,:) = {'showPlot','boolean',false};
defaults(end+1,:) = {'threshold','float',100};

pars = ParseVariableArguments(varargin,defaults,mfilename);

[nReads,~,nSpots] = size(distMapFilt);
%% find TAD boundaries on individual maps
% this gives better results for population substructure
ds = zeros(nReads,nSpots);
us = zeros(nReads,nSpots);
w = pars.w; % window size in bins
for s = 1:nSpots
    map = distMapFilt(:,:,s);
    for r=1:nReads % r = 20
        st = max(1,r-w);
        en = min(nReads,r+w);
        downstream = map(r,st:r-1);
        upstream = map(r,r+1:en);
        ds(r,s) = nanmedian(downstream);
        us(r,s) = nanmedian(upstream);
    end
end
boundaryMap = (ds(:,:)-us(:,:))';
boundaryMap(isnan(boundaryMap)) = 0;
if pars.showPlot
 figure(4); clf; imagesc(boundaryMap);
 GetColorMap('RedWhiteBlue'); colorbar; caxis([-50,50]);
end

%%
map2s = boundaryMap; 
t = pars.threshold;
map2s(map2s<-1*t) = -1*t;
map2s(map2s> t) =  t;
map2s(map2s>-1*t & map2s<t) = 0;

jumpPoints = zeros(nSpots,1);
for s= 1:nSpots
    score = nan(nReads,1);
    for r=1:nReads
        score(r) = sum(map2s(s,1:r)) - sum(map2s(s,r+1:end));        
    end
    [~,jumpPoints(s)] = min(score);
end

% [mids,sortIdx] = sort(jumpPoints);
% mapOut = map2s(sortIdx,:);
