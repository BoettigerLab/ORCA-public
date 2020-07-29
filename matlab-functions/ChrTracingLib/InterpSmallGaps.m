function [smPoly,smMap,exceedMaxBlanks] = InterpSmallGaps(smPoly,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'maxGap','nonnegative',3};
defaults(end+1,:) = {'maxBlanks','nonnegative',.3};
defaults(end+1,:) = {'badHybes','freeType',[]};
defaults(end+1,:) = {'w','fraction',[]};
defaults(end+1,:) = {'method',{'lowess','rlowess','loess','rloess'},'rlowess'};
defaults(end+1,:) = {'showPlots','boolean',false};
defaults(end+1,:) = {'showEvery','integer',1};
pars = ParseVariableArguments(varargin,defaults,mfilename);

smPoly = smPoly(:,1:3,:);
nBins = size(smPoly,1);
nCells = size(smPoly,3);

% maxGap= 3;
if isempty(pars.w)
    pars.w =  3/nBins;% .1; % .08  .05
end

if pars.showPlots
    f1 = figure();
    f2 = figure();
end

smMap = nan(nBins,nBins,nCells);
exceedMaxBlanks = false(nCells,1);
for n=1:nCells
    smPoly(pars.badHybes,:,n) = nan;
    skip = isnan(smPoly(:,1,n));
    skipProp = regionprops(skip,'Area','PixelIdxList');
    bigGaps = [skipProp.Area] > pars.maxGap;
    gap = cat(1,skipProp(bigGaps).PixelIdxList);
    notGap = 1:nBins;
    notGap(gap) = [];
    m1 = nan(nBins,nBins);
    if sum(~skip) > nBins*pars.maxBlanks
        p1 =  [smooth(smPoly(:,1,n),pars.w,pars.method),smooth(smPoly(:,2,n),pars.w,pars.method),smooth(smPoly(:,3,n),pars.w,pars.method)];
        p1(gap,:) = nan;
        if nargout > 1
            m1(notGap,notGap) = squareform(pdist(p1(notGap,:))); 
        end
        % m1 = squareform(pdist(p1));
        if pars.showPlots && rem(n,pars.showEvery)==0
            figure(f1); clf;
            PlotPolymerTube(p1,'tubeRadius',2,'shortestPath',true,'showTube',true); 

            figure(f2); clf;
            m0 = nan(nBins,nBins);
            m0(~skip,~skip) = squareform(pdist( smPoly(~skip,:,n) ));
            subplot(1,2,1); imagesc(m0); %  caxis([50,450]);
            subplot(1,2,2); imagesc(m1); %  caxis([50,450]);
        end
        smPoly(:,:,n) = p1;    
    else
        exceedMaxBlanks(n) = true;
    end
    if nargout > 1
        smMap(:,:,n) = m1;
    end
end


% 
% figure(1); clf;
% subplot(1,2,1); imagesc(nanmedian(segDist,3)); colorbar;
% map2 = nanmedian(smMap,3);
% subplot(1,2,2); imagesc(map2); colorbar;
% 
% figure(2); clf;
% subplot(1,2,1); imagesc(log2(map1)); colorbar;
% map2 = ContactFrac(smMap);
% subplot(1,2,2); imagesc(log2(map2)); colorbar;
