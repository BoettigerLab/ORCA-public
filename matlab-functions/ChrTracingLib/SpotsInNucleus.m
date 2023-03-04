function [fovTable,croppedSpots] = SpotsInNucleus(im3D,spotXY,cellBorders,varargin)
%
%% inputs
% im3D a 3D image  (uncropped)
% spotXY - Nx2xC array of all spots found in the image
% cellBorders - cell array of XY coordinateds of cell borders
%% Outputs
% fovTable - table of 3DPSF fitted spots

defaults = cell(0,3);
defaults(end+1,:) = {'maxSpotsPerCell','integer',2};
defaults(end+1,:) = {'xyCropWindow','integer',6}; % expand by this sq-radius the area cropped
defaults(end+1,:) = {'brightWindow','integer',1}; % expand by this sq-radius the area used for brightness peak
defaults(end+1,:) = {'zCropWindow','integer',10};
defaults(end+1,:) = {'troubleshoot','boolean',false};
defaults(end+1,:) = {'refineFit','boolean',true};
% fit pars
defaults(end+1,:) = {'keepBrightest', 'integer', 1};
defaults(end+1,:) = {'maxFitWidth', 'positive', 8};
defaults(end+1,:) = {'maxFitZdepth', 'positive', 14};
defaults(end+1,:) = {'initSigmaXY','positive',1.25};
defaults(end+1,:) = {'initSigmaZ','positive',2.5};
defaults(end+1,:) = {'minHBratio','nonnegative',1.2}; % peak value over background value
defaults(end+1,:) = {'minAHratio','nonnegative',.25}; % fitted height over background vs peak value
defaults(end+1,:) = {'maxUncert','nonnegative',2}; % pixels
pars = ParseVariableArguments(varargin,defaults,mfilename);

% assign to nuclei
nCells = length(cellBorders); 
fovTable = cell(nCells,1);
croppedSpots = cell(nCells,pars.maxSpotsPerCell);
[ys,xs,zs] = size(im3D);
im_max = max(im3D,[],3);
for c=1:nCells % c=11
    in = inpolygon(spotXY(:,1),spotXY(:,2),...
                   cellBorders{c}(:,1),cellBorders{c}(:,2));
    currSpots = spotXY(in,:);
    if pars.troubleshoot
        figure(4); clf; plot(spotXY(:,1),spotXY(:,2),'r.'); hold on;
        plot(cellBorders{c}(:,1),cellBorders{c}(:,2),'k.-');
        pause(.01);
    end
    % dt = delaunay(cellBorders{c}(:,1),cellBorders{c}(:,2),cellBorders{c}(:,3));
    % inConvexHull = tsearchn(cellBorders{c},dt,spotsXY);
    
    % scan current spots, keep only brightest from 2D projection
    numSpots = size(currSpots,1);
    brightness = zeros(numSpots,1);
    for s=1:numSpots
        % take a small window
        y1 = min(ys,max(1,-pars.brightWindow + round(currSpots(s,2))));
        y2 = min(ys,max(1, pars.brightWindow + round(currSpots(s,2))));
        x1 = min(xs,max(1,-pars.brightWindow + round(currSpots(s,1))));
        x2 = min(xs,max(1,-pars.brightWindow + round(currSpots(s,1))));
        brightBox = im_max(y1:y2,x1:x2);
        brightness(s) = max(brightBox(:));
    end
    [~,idx] = sort(brightness,'descend');
    spotsPerCell = min([numSpots,pars.maxSpotsPerCell]);
    selSpots = idx(1:spotsPerCell);
    bright = brightness(selSpots);
    x = currSpots(selSpots,1);
    y = currSpots(selSpots,2);
    % crop the 3D image to record the z-position and do fine fit
    z = zeros(spotsPerCell,1);
    fineFit = false(spotsPerCell,1); 
    for s=1:spotsPerCell
        x1 = max(1 ,floor(x(s)-pars.xyCropWindow+1));
        x2 = min(xs,floor(x(s)+pars.xyCropWindow));
        y1 = max(1 ,floor(y(s)-pars.xyCropWindow+1));
        y2 = min(ys,floor(y(s)+pars.xyCropWindow+1));
        crop = im3D(y1:y2,x1:x2,:);
         
        % just for troubleshooting
        if pars.troubleshoot
           figure(3); clf; 
           ProjectIm3D(crop,'showPlots',true);
           colormap(gray);
           pause;
        end
        % get peak in z
        yZstack = min(ys,max(1,round(y(s))));
        xZstack = min(xs,max(1,round(x(s))));
        zPixStack = squeeze(im3D(yZstack,xZstack,:));
        % figure(10); clf; plot(zPixStack);
        [~,z(s)] = max(zPixStack);
        z1 = max(1,floor(z(s)-pars.zCropWindow+1));
        z2 = min(zs,floor(z(s)+pars.zCropWindow));
        croppedSpots{c,s} = crop(:,:,z1:z2);
        
        % fit PSF-3D
        try
        if pars.refineFit         
%              
%             figure(1); 
%             subplot(3,2,1); cla; imagesc(xyC); colorbar;
%             subplot(3,2,2); cla; imagesc(xzC); colorbar;
%             [xyC,xzC] = ProjectIm3D(croppedSpots{c,s});
%             subplot(3,2,3); cla; imagesc(xyC); colorbar;
%             subplot(3,2,4); cla; imagesc(xzC); colorbar;
            dTable = FitPsf3D(croppedSpots{c,s},'parameters',pars,'troubleshoot',false); 
            % adjust to FOV coordinates
             % keep units in pixels, and lets record the pix to nm in the table.
             if ~isempty(dTable)
                x(s) = dTable.x + x1-1; 
                y(s) = dTable.y + y1-1; 
                z(s) = dTable.z + z1-1; 
                bright(s) = dTable.h;
                fineFit(s) = true;
%                 [xyC,xzC] = ProjectIm3D(im3D);
%                 figure(11); clf; 
%                 subplot(1,2,1); imagesc(xyC);
%                 hold on; plot(x(s),y(s),'r.');
%                 subplot(1,2,2); imagesc(xzC);
%                 hold on; plot(x(s),z(s),'r.');
             end
        end
        catch er
            disp(er.getReport);
            disp('debug here');
        end
    end
%     figure(4); clf; imagesc(IncreaseContrast(im_max,'high',.999)); colormap(gray); hold on; colorbar;
%     plot(x,y,'r*');
    
    cellNum = ones(spotsPerCell,1)*c;
    cellTable = table(x,y,z,bright,cellNum,fineFit);
%  additional data best recorded in the parent function: fov, hyb
%     fovNum = ones(spotsPerCell,1)*fov;
%     hybNum = ones(spotsPerCell,1)*hyb;
%     cellTable = table(x,y,z,bright,cellNum,fovNum,hybNum);
    fovTable{c} = cellTable;
end
fovTable = cat(1,fovTable{:});

