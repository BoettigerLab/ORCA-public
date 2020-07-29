function [normMap,rawMap] = PlotEmbryoHiC(locustxt,varargin)
% Plot the Schwartz and Cavalli 2014 5kb res Hi-C data

% data from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61471

defaults = cell(0,3);
defaults(end+1,:) = {'flyHiCbinCoordsPath', 'string', 'C:\Data\Fly\Embryo\Cavalli2014\GSE61471_nm_none_5000.bins.txt'};
defaults(end+1,:) = {'flyHiCbinContactPath', 'string', 'C:\Data\Fly\Embryo\Cavalli2014\GSE61471_nm_none_5000.n_contact.txt'};
defaults(end+1,:) = {'fontSize', 'positive', 14};
defaults(end+1,:) = {'showplot','boolean',true};
defaults(end+1,:) = {'colormap','colormap',GetColorMap('whiteToRed')};
defaults(end+1,:) = {'units',{'genome','bins'},'genome'};
defaults(end+1,:) = {'binSize','positive',5000}; % rebin data
defaults(end+1,:) = {'resizeMethod',{'bilinear','nearest'},'bilinear'}; % rebin data


pars = ParseVariableArguments(varargin,defaults,mfilename);

global flyHiCbinCoords flyHiCbinContact ;

if isempty(flyHiCbinCoords)
flyHiCbinCoords = readtable(pars.flyHiCbinCoordsPath);
flyHiCbinContact = readtable(pars.flyHiCbinContactPath);
end


% locustxt = 'chr3R:12635000-12765000'; % abd
% locustxt = 'chr2R:7330000-7490000'; % en
% locustxt = 'chr2R:7230000-7590000'; % en +/- 100 kb

% Parse Locus Name
[chr,locusStart,locusEnd] = ParseLocusName(locustxt);


onChr = StringFind(flyHiCbinCoords.chr,regexprep(chr,'chr',''),'boolean',true);
inRange = flyHiCbinCoords.from_coord >= locusStart & flyHiCbinCoords.to_coord <= locusEnd;
cbins = flyHiCbinCoords.cbin(onChr & inRange);
% clc;
% steps = flyHiCbinCoords.to_coord(onChr & inRange) - flyHiCbinCoords.from_coord(onChr & inRange);
% figure(7); clf; plot(steps,'.');
% figure(7); clf; plot(diff(cbins));
selectBins = ismember(flyHiCbinContact.cbin1, cbins) & ismember(flyHiCbinContact.cbin2, cbins);
locusSizeInBins = round((locusEnd - locusStart + 1)/5000);
blankMap = nan(locusSizeInBins);
idx = sub2ind(size(blankMap), flyHiCbinContact.cbin1(selectBins)-cbins(1)+1,flyHiCbinContact.cbin2(selectBins)-cbins(1)+1);
rawMap = blankMap;
expectMap = blankMap;
rawMap(idx) = flyHiCbinContact.observed_count(selectBins);
expectMap(idx) = flyHiCbinContact.expected_count(selectBins);

normMap = rawMap./expectMap; 

if pars.binSize ~= 1
    displayBins = length(cbins)*5000/pars.binSize;
    rawMap =  imresize(rawMap,[displayBins,displayBins],'bilinear'); 
    normMap = imresize(normMap,[displayBins,displayBins],'bilinear'); 
else
    displayBins = length(cbins);
end

if nargout == 0 || pars.showplot % figure(5);  clf;
    imagesc(normMap); colormap(pars.colormap); colorbar;
    if strcmp(pars.units,'genome')
        ticks = cellfun(@(x) num2str(x),{round(100*locusStart/1e6)/100,round(100*locusEnd/1e6)/100},'UniformOutput',false);
        set(gca,'Xtick',[1,displayBins],'XTickLabel',ticks,'FontSize',pars.fontSize);
        set(gca,'Ytick',[1,displayBins],'YTickLabel',ticks,'FontSize',pars.fontSize);
    end
end