% Plots a distance map with spots selected based on RNA expression.
% Uses spots with RNA expression between quantileLow and quantileHigh.
% Throws out spots with RNA expression outside of minRNA and maxRNA.
% Returns number of distance maps compiled for final image.

function n = PlotRNAFilteredDistanceMap(rnaTables, distMapFilts, varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'FOVs','array',1};  % FOVs to analyze
defaults(end+1,:) = {'RNAIndex','positive',1};  % Column number of RNA column (usually AbdA=1, AbdB = 2)
defaults(end+1,:) = {'minRNA','nonnegative',0};  % Throws out RNA expression below minRNA and maxRNA (before finding quantile vals)
defaults(end+1,:) = {'maxRNA','nonnegative',realmax};
defaults(end+1,:) = {'quantileLow','nonnegative',0};  % Will filter out RNA outside of quantileLow and quantileHigh
defaults(end+1,:) = {'quantileHigh','nonnegative',1};
defaults(end+1,:) = {'badHybes','positive',[]};
defaults(end+1,:) = {'caxis','array',[50,450]};
pars = ParseVariableArguments(varargin,defaults,mfilename);


% Calculate threshold values
allRNA = [];
for f = pars.FOVs
    allRNA = [allRNA; rnaTables{f}.(pars.RNAIndex)];
end
allRNA(allRNA < pars.minRNA) = [];
allRNA(allRNA > pars.maxRNA) = [];
highRNAThresh = quantile(allRNA, pars.quantileHigh);
lowRNAThresh  = quantile(allRNA, pars.quantileLow);

% Gemerate maps for selected spots
mapRNA_allCells = [];
for f = pars.FOVs
    rna_f = rnaTables{f}.(pars.RNAIndex);
   
    rnaKeep = lowRNAThresh <= rna_f & rna_f <= highRNAThresh;
    mapRNAKeep = distMapFilts{f}(:,:,rnaKeep);
    mapRNA_allCells = cat(3, mapRNA_allCells, mapRNAKeep);
end

% Plot
medianMapRNA = nanmedian(mapRNA_allCells,3);
sx = size(medianMapRNA, 1);
keepHybes = 1:sx;
keepHybes(pars.badHybes) = [];
imagesc(medianMapRNA(keepHybes, keepHybes));
GetColorMap('RedWhiteBlueSat'); 
caxis(pars.caxis);
colorbar();

n = sum(lowRNAThresh <= allRNA & allRNA <= highRNAThresh);