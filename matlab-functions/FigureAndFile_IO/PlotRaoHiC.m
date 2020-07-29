function [hicLocal,ticks] = PlotRaoHiC(locusName,varargin)
% 
%%  Examples
%   Display map from this chr21 region in IMR90 at default 30 kb res
% figure(1); clf; 
% PlotRaoHiC('chr21:29372319-31372318');
%   Display the K562 map at 5 kb res
% figure(1); clf;
% PlotRaoHiC('chr21:29372319-31372318','displayRes',5E3,'dataset','K562');
% 
% 


% globals
global hicIMR90  hicK562 currentChr

% defaults
defaults = cell(0,3);
defaults(end+1,:) = {'hicDataFolder','string','U:\GenomeData\Human\Hi-C\'};
defaults(end+1,:) = {'dataset',{'IMR90','K562'},'IMR90'};
defaults(end+1,:) = {'mapRes','positive',5000}; % in kb.  must match map
defaults(end+1,:) = {'displayRes','positive',30E3};
defaults(end+1,:) = {'interpMethod',{'nearest','bilinear'},'nearest'};
defaults(end+1,:) = {'colormap','colormap',GetColorMap('default')};
defaults(end+1,:) = {'fontSize', 'positive', 14};
defaults(end+1,:) = {'units',{'genome','bins'},'genome'};
defaults(end+1,:) = {'showPlots','boolean',true};

pars = ParseVariableArguments(varargin,defaults,mfilename);

%% Comparison to Hi-C data
mapRes = 5000;
[chr,locusStart,locusEnd] = ParseLocusName(locusName);
if contains(pars.dataset,'IMR90')
    HiC_file =[pars.hicDataFolder, 'IMR90\5kb_resolution_intrachromosomal\',chr,'\MAPQG0\',chr,'_5kb.RAWobserved'];
    if isempty(hicIMR90) || ~strcmp(currentChr,chr)
        hicIMR90 = ReadRaoData(HiC_file,mapRes);
    end
    hicMap = hicIMR90;
elseif contains(pars.dataset,'K562')
    HiC_file =[pars.hicDataFolder, 'K562\5kb_resolution_intrachromosomal\',chr,'\MAPQG0\',chr,'_5kb.RAWobserved'];
    if isempty(hicK562) || ~strcmp(currentChr,chr)
        hicK562 = ReadRaoData(HiC_file,mapRes);
    end
    hicMap = hicK562;
else
    error([pars.dataset,' not currently a valid HiC-file']);
end



c0 = round(locusStart/mapRes);
c1 = round(locusEnd/mapRes);
displayBins = round((locusEnd-locusStart)/pars.displayRes);

% rebin at 'displayRes' resolution
if pars.displayRes ~= 5000
    hicLocal = imresize(double(hicMap(c0:c1,c0:c1)),[displayBins,displayBins],pars.interpMethod);%  ,'bilinear'
else
    hicLocal = double(hicMap(c0:c1,c0:c1)); 
end

if pars.showPlots
    hicFig = gcf; clf;
    imagesc(log10(hicLocal));
    cbar = colorbar;
    title('Hi-C raw reads');
    ylabel(cbar,'Read Counts');
    colormap(pars.colormap); 
    hicFig.Color = 'w';
    hicFig.Name = 'HiC_Fig';

    if strcmp(pars.units,'genome')
        ticks = cellfun(@(x) num2str(x),{round(100*locusStart/1e6)/100,round(100*locusEnd/1e6)/100},'UniformOutput',false);
        set(gca,'Xtick',[1,displayBins],'XTickLabel',ticks,'FontSize',pars.fontSize);
        set(gca,'Ytick',[1,displayBins],'YTickLabel',ticks,'FontSize',pars.fontSize);
        xlabel(chr); ylabel(chr);
    end
end


function hicMap = ReadRaoData(HiC_file,mapRes) 
    hicData = dlmread(HiC_file,'\t');
    lenChr = max(hicData(:,2))/mapRes; % in mapRes (5kb) intervals
    hicMap = zeros(lenChr);
    for n=1:length(hicData(:,1))
        i = round(hicData(n,1)/mapRes);
        j = round(hicData(n,2)/mapRes);
        hicMap(i,j) = hicData(n,3);
        hicMap(j,i) = hicData(n,3);
    end


