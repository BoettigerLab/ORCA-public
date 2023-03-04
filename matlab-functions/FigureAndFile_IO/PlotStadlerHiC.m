function map = PlotStadlerHiC(locustxt,varargin)
% Plot the Stadler HiC data  eLife 2017
%  Note, data is aligned and mapped to dm3

global stadlerHiC


defaults = cell(0,3);
defaults(end+1,:) = {'dataPath', 'string', 'U:\GenomeData\ByPublication\Stadler2017\stage5_whole_rep1_5kb.txt'}; % stage5Stadler5kb
defaults(end+1,:) = {'fontSize', 'positive', 14};
defaults(end+1,:) = {'showplot','boolean',true};
defaults(end+1,:) = {'colormap','colormap',GetColorMap('whiteToRed')};
defaults(end+1,:) = {'units',{'genome','bins'},'genome'};
defaults(end+1,:) = {'binSize','positive',5000}; % rebin data
defaults(end+1,:) = {'resizeMethod',{'bilinear','nearest'},'bilinear'}; % rebin data
pars = ParseVariableArguments(varargin,defaults,mfilename);


if ~istable(stadlerHiC)
    stadlerHiC  = readtable(pars.dataPath);
end
%% Hi-C comparison
% locustxt = 'chr3R:2,447,305-2,920,782'; % ANTC  reg 1 - 94
hicTable = stadlerHiC;
stp = 5E3; % hi-C stepSize
[chr,st,en] = ParseLocusName(locustxt,'removeChr',false);

inChr = strcmp(hicTable.chrom1,chr) & strcmp(hicTable.chrom2,chr);
inBox = hicTable.start1 > st & hicTable.start1 < en & hicTable.start2 > st & hicTable.start2 < en;
localTable = hicTable(inChr & inBox,:);

nbins = ceil( (en-st)/stp);
map = zeros(nbins,nbins);
rows = ceil((localTable.start1 -st)/stp);
cols = ceil((localTable.start2 -st)/stp);
idx = sub2ind([nbins,nbins],rows,cols);
map(idx) = localTable.balanced;
map = map + map';


if nargout == 0
    imagesc(map); colorbar;
end