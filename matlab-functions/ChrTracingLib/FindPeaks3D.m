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
defaults(end+1,:) = {'maxFitZdepth', 'positive', 14};
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
pars = ParseVariableArguments(varargin, defaults, mfilename);

warning('off','MATLAB:singularMatrix');

% -------------------------------------------------------------------------
% Main Function
% -------------------------------------------------------------------------
[rows,cols,stcks] = size(Im1);
%  THIS doesn't seem right
wb = min([floor(pars.maxFitWidth/2),floor(rows/2),floor(cols/2)]); % short hand for crop;
wz = min([floor(pars.maxFitZdepth/2),floor(stcks/2)])-1; % short hand for crop;

% defaults to maxFitWidth 
% wb = floor(rows/2)-floor(pars.maxFitWidth/2);
wb = min([wb,rows,cols]); % threshold at max rows or max cols if wb>rows
wb = max([0,wb]); % threshold at min of 0
% wz = floor(rows/2)-floor(pars.maxFitZdepth/2);
wz = min([wz,rows,cols]); % threshold at max rows or max cols if wb>rows
wz = max([0,wz]); % threshold at min of 0

if isinf(pars.minPeakHeight)
    pars.minPeakHeight = quantile(Im1(:),.99);
end
% find local maxima to pixel precision
% Im3 = Im1(wb+1:end-wb+1,wb+1:end-wb+1,wz+1:end-wz+1); % exclude points on border OLD version  
Im3 = Im1(wb+1:end-wb,wb+1:end-wb,wz+1:end-wz); % exclude points on border 
dtype = class(Im3);
Im3 = Im3 - cast(pars.cameraBackground,dtype);
if pars.peakBlur ~= 0
    Im3 = imgaussfilt3(Im3,pars.peakBlur);
end
bw = imregionalmax(Im3);
if ~isempty(pars.relativeHeight) && pars.relativeHeight ~= 0
    maxHeight = max(Im3(:));
    bw(Im3< pars.relativeHeight*maxHeight) = 0;
end
% figure(100); clf; imagesc(max(Im1,[],3)); colorbar;
% figure(100); clf; imagesc(max(Im3,[],3)); colorbar;
% figure(101); clf; imagesc(max(bw,[],3));
bw(Im3<pars.minPeakHeight) = 0; % only keep peaks above min Threshold
connMap = bwconncomp(bw);


nMax = connMap.NumObjects;
xyzh = zeros(nMax,4);
for n=1:nMax
    [y,x,z] = ind2sub(size(bw),connMap.PixelIdxList{n});
    xyzh(n,:) = mean([x,y,z,Im3(connMap.PixelIdxList{n})],1);
end

nMax = min(nMax,pars.keepBrightest);
if nMax > 0
    [~,i] = sort(xyzh(:,4),'descend');
    keepSpots = i(1:nMax);
    % keepSpots = xyzh(:,4) >= max(xyzh(:,4))*parameters.keepBrightest;
    xyzh = xyzh(keepSpots,:); 
end
nMax = size(xyzh,1);

if pars.minSep > 0
    % deal with split maxima by combining spots within a pixel distance of 3
    % in a perfect function this would be iterated, since after averaging the
    % points may again be within the min cut-off distance from one another.
    % in practice I'm mostly dealing with 2 points to start with. 
    % xyzh = rand(9,4)*6; % some test points
    toDelete = [];
    ds = squareform(pdist(xyzh(:,1:3)));
    for n=1:nMax-1
       tooClose = n+find(ds(n,n+1:end) < pars.minSep); % indices of points within minSep from point n.  
       combine = [n,tooClose];
       xyzh(n,:) = mean(xyzh(combine,:),1); % average over rows
       toDelete = [toDelete,tooClose]; %#ok<AGROW>
       if pars.veryVerbose
          disp('combining spots ');
          disp(combine); 
       end
    end
    xyzh(toDelete,:) = [];
end

x = round(xyzh(:,1)+wb); % we did some averaging, but we need this to be to pixel precision only.  
y = round(xyzh(:,2)+wb);
z = round(xyzh(:,3)+wz); 
h = xyzh(:,4); % so variable names will be auto imported by table; 
dTable = table(x,y,z,h);

if pars.troubleshoot
    % plot local maxima to pixel precision
    subplot(2,2,1); imagesc(max(Im1,[],3)); hold on; plot(dTable.x,dTable.y,'o','color',[1 .5 .5]);  % plots x vs y
    % squashed the columns (y-dim); put the stack (z-dim) on the y-axis (new columns);  plots x vs z;
    subplot(2,2,2); imagesc(max(permute(Im1,[3,2,1]),[],3)); hold on; plot(dTable.x,dTable.z,'o','color',[1 .5 .5]); 
end

