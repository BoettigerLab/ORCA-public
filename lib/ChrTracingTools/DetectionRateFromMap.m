function [m,s] = DetectionRateFromMap(map,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'showPlot','boolean',true};
defaults(end+1,:) = {'minFracValidData','fraction',1/5};
pars = ParseVariableArguments(varargin,defaults,mfilename);

[nReads,~,nCells] = size(map);
detectCnt = squeeze(nansum(map>-1));
obsPerCell = sum(detectCnt,1);
detectCnt(:,obsPerCell<nReads*pars.minFracValidData*nReads) = []; % remove secondary points from save FOV and unfiltered background
detectRate =  detectCnt>1; % 
% figure(3); clf; violin(detectRate','method','pchip');
m = nanmean(detectRate,2);
s = nanstd(detectRate,[],2)./sqrt(sum(detectRate,2));
if pars.showPlot
    bar(m); 
    hold on; ploterr(1:length(m),m,[],s,'.','color','k');
    ylim([0,1]);
    set(gcf,'color','w');
    title(['nCells = ',num2str(nCells)]);
end