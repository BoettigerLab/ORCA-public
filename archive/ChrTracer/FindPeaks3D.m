function dTable = FindPeaks3D(Im1,varargin)

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
% key parameters
defaults(end+1,:) = {'minPeakHeight', 'positive', 1000};
defaults(end+1,:) = {'cameraBackground', 'nonnegative', 0};
defaults(end+1,:) = {'peakBlur', 'nonnegative', .5}; % gaussian smoothing before initial find max 
defaults(end+1,:) = {'maxFitWidth', 'positive', 8};
defaults(end+1,:) = {'keepBrightest','integer',inf};
defaults(end+1,:) = {'relativeHeight','fraction',0};
% less common parameters
defaults(end+1,:) = {'minSep','positive',3}; % minimum separation in pixels between initial maxima
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'veryVerbose','boolean',false};
defaults(end+1,:) = {'troubleshoot','boolean',true};

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

warning('off','MATLAB:singularMatrix');

% -------------------------------------------------------------------------
% Main Function
% -------------------------------------------------------------------------
w = floor(parameters.maxFitWidth/2); % short hand for crop;
if isinf(parameters.minPeakHeight)
    parameters.minPeakHeight = quantile(Im1(:),.99);
end
% find local maxima to pixel precision
Im3 = Im1(w+1:end-w+1,w+1:end-w+1,w+1:end-w+1); % exclude points on border
dtype = class(Im3);
Im3 = Im3 - cast(parameters.cameraBackground,dtype);
Im3 = imgaussfilt3(Im3,parameters.peakBlur);
bw = imregionalmax(Im3);
if ~isempty(parameters.relativeHeight) && parameters.relativeHeight ~= 0
    maxHeight = max(Im3(:));
    bw(Im3< parameters.relativeHeight*maxHeight) = 0;
end
bw(Im3<parameters.minPeakHeight) = 0; % only keep peaks above min Threshold
connMap = bwconncomp(bw);


nMax = connMap.NumObjects;
xyzh = zeros(nMax,4);
for n=1:nMax
    [y,x,z] = ind2sub(size(bw),connMap.PixelIdxList{n});
    xyzh(n,:) = mean([x,y,z,Im3(connMap.PixelIdxList{n})],1);
end

nMax = min(nMax,parameters.keepBrightest);
if nMax > 0
    [~,i] = sort(xyzh(:,4),'descend');
    keepSpots = i(1:nMax);
    % keepSpots = xyzh(:,4) >= max(xyzh(:,4))*parameters.keepBrightest;
    xyzh = xyzh(keepSpots,:); 
end
nMax = size(xyzh,1);

if parameters.minSep > 0
    % deal with split maxima by combining spots within a pixel distance of 3
    % in a perfect function this would be iterated, since after averaging the
    % points may again be within the min cut-off distance from one another.
    % in practice I'm mostly dealing with 2 points to start with. 
    % xyzh = rand(9,4)*6; % some test points
    toDelete = [];
    ds = squareform(pdist(xyzh(:,1:3)));
    for n=1:nMax-1
       tooClose = n+find(ds(n,n+1:end) < parameters.minSep); % indices of points within minSep from point n.  
       combine = [n,tooClose];
       xyzh(n,:) = mean(xyzh(combine,:),1); % average over rows
       toDelete = [toDelete,tooClose]; %#ok<AGROW>
       if parameters.veryVerbose
          disp('combining spots ');
          disp(combine); 
       end
    end
    xyzh(toDelete,:) = [];
end

x = round(xyzh(:,1)+w); % we did some averaging, but we need this to be to pixel precision only.  
y = round(xyzh(:,2)+w);
z = round(xyzh(:,3)+w); 
h = xyzh(:,4); % so variable names will be auto imported by table; 
dTable = table(x,y,z,h);

if parameters.troubleshoot
    % plot local maxima to pixel precision
    subplot(1,2,1); imagesc(max(Im1,[],3)); hold on; plot(dTable.x,dTable.y,'o','color',[1 .5 .5]);  % plots x vs y
    % squashed the columns (y-dim); put the stack (z-dim) on the y-axis (new columns);  plots x vs z;
    subplot(1,2,2); imagesc(max(permute(Im1,[3,2,1]),[],3)); hold on; plot(dTable.x,dTable.z,'o','color',[1 .5 .5]); 
end

