function [normData,corMap] = MapCompartment(medMap,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'showPlots','boolean',true};
defaults(end+1,:) = {'start','integer',4};
defaults(end+1,:) = {'stop','integer',8};
defaults(end+1,:) = {'method',{'powerlaw','median'},'powerlaw'};
pars = ParseVariableArguments(varargin,defaults,mfilename);


normMap = ComputeNormMap(medMap,'parameters',pars); % power law fit normalization
normData = medMap./normMap;
normData(isnan(normData)) = 0;
normData(isinf(normData)) = 0;
nReads = size(normData,1);
corMap = zeros(nReads,nReads);
for i=1:nReads
    for j=1:nReads
        corMap(i,j) = corr(normData(i,:)',normData(:,j));
    end
end

if pars.showPlots
    figure(1); clf; imagesc(normData); caxis([.7 1.3]);
    colormap(flipud(GetColorMap('RedWhiteBlueK'))); colorbar;
    set(gcf,'color','w');


    fOut = figure(2); clf; imagesc(corMap); colorbar;  caxis([-1 1]);
    colormap((GetColorMap('RedWhiteBlueK'))); colorbar;
   %  set(gca,'XTick',[1,nBins],'XTickLabel',coords,'YTick',[1,nBins],'YTickLabel',coords,'FontSize',15);
    set(gcf,'color','w');
    axis image; colorbar;
end