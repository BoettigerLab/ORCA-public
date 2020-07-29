function [mRNAcents,stainCents,stainCounts,clusterStats] = FindMRNA(xf,yf,varargin)
% Clusters mRNA and returns centroids
% cents = FindMRNA(xf,yf,'binSize',binSize,'minDotPerBin',minDotPerBin,...
%         'minLocsPerDot',minLocsPerDot,'minArea',minArea,'maxArea',maxArea);
% 
%--------------------------------------------------------------------------
% Required parameters
% xf and yf -- cell arrays of the filtered x and y STORM localizations.
% 
%--------------------------------------------------------------------------
% Optional Parameters
% 
% binSize = 16;  % size of bins in nm
% minDotPerBin = 2;  % min number of localizations to call a bin occupied
% minLocsPerDot = 0; % min number of localization in all bins assigned to a cluster to be called an mRNA
% minArea = 0;   % min area in bins to be called a cluster of localization
% maxArea = 300; % max area in bins to be called a cluster of localization
% 
%--------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'binSize', 'positive', 16};
defaults(end+1,:) = {'minDotPerBin', 'nonnegative',2 };
defaults(end+1,:) = {'minLocsPerDot', 'positive',0 };
defaults(end+1,:) = {'minArea', 'nonnegative',0 };
defaults(end+1,:) = {'maxArea', 'nonnegative',300 };
defaults(end+1,:) = {'showPlots', 'boolean', true};
defaults(end+1,:) = {'verbose', 'boolean', true};
defaults(end+1,:) = {'clusterFig', 'handle',[] };
defaults(end+1,:) = {'histFig', 'handle',[] };
defaults(end+1,:) = {'imageSize', 'array',[256,256] };
defaults(end+1,:) = {'npp', 'positive',160 };

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2
    error('matlabSTORM:invalidArguments', 'two nx2 vectors of points are required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);
% parameters = ParseVariableArguments('', defaults, mfilename);

H = parameters.imageSize(1);
W = parameters.imageSize(2);

%%
D = length(xf); 

% Identify "true" mRNA (based on minLocsPerDot)
% Clusters localization by combining all hybes
xall = cat(1,xf{:});
yall = cat(1,yf{:});
cluster_scale = parameters.npp/parameters.binSize; 
M = hist3([yall,xall],{1/cluster_scale:1/cluster_scale:H,1/cluster_scale:1/cluster_scale:W});
P = regionprops(M>=parameters.minDotPerBin,M,'PixelValues','PixelIdxList','Centroid','Area');
mRNAcents = cat(1,P.Centroid)/cluster_scale;
clusterLocs = cellfun(@sum,{P.PixelValues});
clusterarea = [P.Area];
goodClusters = clusterLocs>parameters.minLocsPerDot & clusterarea > parameters.minArea & clusterarea < parameters.maxArea;
mRNAcents = mRNAcents(goodClusters,:);

if parameters.verbose
    tooFewLocs = sum(clusterLocs < parameters.minLocsPerDot)/length(clusterLocs)*100;
    tooSmallSpot = sum(clusterLocs < parameters.minArea)/length(clusterLocs)*100;
    tooBigSpot = sum(clusterLocs > parameters.maxArea)/length(clusterLocs)*100; 
    disp(['Rejected ',num2str(tooFewLocs),'% of clusters for too few localizations']);
    disp(['Rejected ',num2str(tooSmallSpot),'% of clusters for too small a cluster']);
    disp(['Rejected ',num2str(tooBigSpot),'% of clusters for too large a cluster']);
end

if parameters.showPlots
    if ~isempty(parameters.histFig)
    %     figure(1); clf; imagesc(M); caxis([0,30]); colormap(jet); colorbar; 
    trueMRNA = bwlabel(M>=parameters.minDotPerBin); 
    noise = ~(goodClusters);
    allPix = {P.PixelIdxList};
    allPix = cell2mat(cat(1,allPix(noise)'));
    trueMRNA(allPix) = 0;
    figure(parameters.histFig); clf; imagesc(trueMRNA); pause(.1); 
    end
    
   if ~isempty(parameters.clusterFig)       
       figure(parameters.clusterFig); 
       plot(mRNAcents(:,1),mRNAcents(:,2),'ko','MarkerSize',14); hold on;
   end
end

%----------------- Cluster localizations in each individual hybe
if nargout > 1 

    Mp = cell(D,1);
    stainCents = cell(D,1);
    stainCounts = cell(D,1);
    for d=1:D
        Mp{d} = hist3([yf{d},xf{d}],...
           {1/cluster_scale:1/cluster_scale:H,1/cluster_scale:1/cluster_scale:W});
        P = regionprops(Mp{d}>0,Mp{d},'PixelValues','Centroid');
        stainCents{d} = cat(1,P.Centroid)/cluster_scale;
        stainCounts{d} = cellfun(@sum,{P.PixelValues});

     % Plots for troubleshooting   
    %      positionsInStain = stainCents{d}(stainCounts{d} > minLocsPerStain,:);
    %     figure(5); clf; plot(xf{d},yf{d},'r.');
    %     hold on; plot(positionsInStain(:,1),positionsInStain(:,2),'k.');
    %     
    end

    % multicolor image on binSize scale of different clusters    
    if ~isempty(parameters.histFig);
        [h,w] = size(M);
        I = zeros(h,w,D);
        for d=1:D
            I(:,:,d) = Mp{d}.*trueMRNA;
        end
        figure(parameters.histFig); clf; Ncolor(uint8(100*I),hsv(D));
    end
end
% 
% legend(num2str([0:D]'))

if nargout == 4
    clusterStats.localizations = clusterLocs;
    clusterStats.area = clusterarea;
    clusterStats.goodClusters = goodClusters;
end







