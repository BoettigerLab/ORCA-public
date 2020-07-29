function PlotConvGeneLocations(alignedIm,mRNAcents,libGenes,libCodes,spotIDs,varargin)


% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'mlists', 'cell', []};
defaults(end+1,:) = {'MarkerSize', 'positive', 10};
defaults(end+1,:) = {'FontSize', 'positive', 6};
defaults(end+1,:) = {'nativeColor', 'boolean', false};
defaults(end+1,:) = {'showNames', 'boolean', false};
defaults(end+1,:) = {'showLegend', 'boolean', true};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 4
    error('matlabSTORM:invalidArguments', 'requires alignedIm,mRNAcents,libGenes,libCodes');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

%%

% figure(7); clf; 
[numGenes,numHybes] = size(libCodes);
clrmap = jet(numHybes);

Ncolor(uint16(alignedIm),clrmap);  hold on;

if parameters.nativeColor
    clrmapGenes = zeros(numGenes,3);
else
    clrmapGenes = jet(numGenes);
    clrmapGenes = clrmapGenes(randperm(numGenes),:);
end

for n=1:numGenes
    if parameters.nativeColor
        clrmapGenes(n,:) = sum(clrmap(libCodes(n,:),:))/4;
    end
    plot(0,0,'o','color',clrmapGenes(n,:));
end

for n=1:numGenes
    geneSpots = spotIDs == n;
    plot(mRNAcents(geneSpots,1),mRNAcents(geneSpots,2),'o','color',clrmapGenes(n,:));
    if parameters.showNames
      text(mRNAcents(geneSpots,1)+1,mRNAcents(geneSpots,2),libGenes{n},'color',clrmapGenes(n,:),'FontSize',parameters.FontSize);
    end
end
if parameters.showLegend
legend(libGenes,'Location','EastOutside');
end

if ~isempty(parameters.mlists)
    for h=1:numHybes
        plot(parameters.mlists{h}.xc,parameters.mlists{h}.yc,'.',...
            'color',clrmap(h,:),'MarkerSize',parameters.MarkerSize);
    end
end

axis image;

