function [correctClusters,spotMatrix] = PlotGeneLocationsConv(mRNAcents,Mcents,uniqueMsg,hammingDis,flists,libGenes,maxDtoCentroid,varargin)


% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'showbkd', 'boolean', true};
defaults(end+1,:) = {'showplots', 'boolean',true};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 7
    error('matlabSTORM:invalidArguments', 'requires mRNAcents,Mcents,uniqueMsg,hammingDis,flists,libGenes');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

showbkd = parameters.showbkd;
showplots = parameters.showplots; 
numGenes = length(libGenes);
numHybes = size(Mcents,2);
clrmap = hsv(numHybes); 


clear spotMatrix;
cntsPerSpot = cell(numGenes,1);
spotMatrix = cell(numGenes,1); % figure(40); clf;
for h = 1:numHybes
    plot(0,0,'.','color',clrmap(h,:),'MarkerSize',10); hold on;
end
colormap(clrmap);
legend(num2str( (1:numHybes)' )); hold on; 
 

correctClusters = cell(numGenes,1);
for i=1:numGenes
    correctClusters{i} = ismember(Mcents,...
            uniqueMsg((hammingDis(:,i)==1|hammingDis(:,i)==0),:),'rows');
end


for h = 1:numHybes
    if showbkd && showplots
        plot(flists{h}.xc,flists{h}.yc,'.','color',clrmap(h,:),'MarkerSize',1); hold on;
    end
    for n=1:numGenes
        [idx,dist] = knnsearch( mRNAcents(correctClusters{n},:) , [flists{h}.xc,flists{h}.yc] );
        inRange = dist<maxDtoCentroid;  
        cntsPerSpot{n}{h} = hist(idx(inRange),1:size(mRNAcents(correctClusters{n},:),1) );
        spotMatrix{n} = cat(1,spotMatrix{n},cntsPerSpot{n}{h});
        
        if showplots
            plot(flists{h}.xc(inRange),flists{h}.yc(inRange),'.','color',clrmap(h,:)); hold on;
            text(mRNAcents(correctClusters{n},1),mRNAcents(correctClusters{n},2),libGenes{n});
        end
    end
end