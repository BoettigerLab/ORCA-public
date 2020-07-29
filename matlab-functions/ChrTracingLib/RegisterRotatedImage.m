function [rotImageOut,tform] = RegisterRotatedImage(refImage,rotImage,varargin)
% [rotImageOut,tform] = RegisterRotatedImage(refImage,rotImage)
% 

defaults = cell(0,3);
defaults(end+1,:) = {'maxD','positive',15};
defaults(end+1,:) = {'useCorrAlign','boolean',false};
defaults(end+1,:) = {'showPlots','boolean',false};
defaults(end+1,:) = {'warp',{'NonreflectiveSimilarity','lwm','polynomial','pwl'},'NonreflectiveSimilarity'};

pars = ParseVariableArguments(varargin,defaults,mfilename);

% Find spots in both images
spots1 = FindSpots(refImage); 
spots2 = FindSpots(rotImage); 

if pars.showPlots
    subplot(1,3,1); 
    Ncolor(cat(3,refImage,rotImage));
    hold on; 
    plot(spots1(:,1),spots1(:,2),'yo');
    plot(spots2(:,1),spots2(:,2),'r>');
end

% Match spots and compute warp. 
[matched1,matched2] = MatchFiducials(spots1,spots2,'useCorrAlign',pars.useCorrAlign,'maxD',pars.maxD,'showPlots',pars.showPlots,'verbose',false);

if strcmp(pars.warp,'lwm')
    tformF = fitgeotrans(spots2(matched2,:),spots1(matched1,:),'lwm',6); % compute warp. uses 12 points per area
elseif strcmp(pars.warp,'polynomial')
    tformF = fitgeotrans(spots2(matched2,:),spots1(matched1,:),pars.warp,3); % compute warp. uses 12 points per area
else
    tformF = fitgeotrans(spots2(matched2,:),spots1(matched1,:),pars.warp); % compute warp
end
% pre-2013b version
% tform2Dorig = cp2tform( spots1(matched1,:),spots2(matched2,:),'polynomial',2); % compute warp
% warped = tforminv(tform2D, spots2(matched2,:)); % apply warp

% Apply the rotation
rotImageOut = imwarp(rotImage,tformF,'OutputView',imref2d(size(refImage))); 

% Compute rotation angle and rotation scale. 
[x,y] = transformPointsInverse(tformF,[0 1],[0 0]);
dx = x(2) - x(1); 
dy = y(2) - y(1); 
tform.rotationAngle =  (180/pi) * atan2(dy, dx);
tform.scaleImage = 1 / sqrt(dx^2 + dy^2);
tform.tobj = tformF;

% for troubleshooting create plot
if pars.showPlots
    subplot(1,3,2); % plot reference points 
    
    if strcmp(pars.warp,'lwm')
        tformI = fitgeotrans(spots1(matched1,:),spots2(matched2,:),'lwm',6); % compute warp. uses 12 points per area
    elseif strcmp(pars.warp,'polynomial')
        tformI = fitgeotrans(spots1(matched1,:),spots2(matched2,:),pars.warp,3); % compute warp. uses 12 points per area
    else
        tformI = fitgeotrans(spots1(matched1,:),spots2(matched2,:),pars.warp); % compute warp
    end
    
    
    warped = transformPointsInverse(tformI,spots2(matched2,:));
    plot(spots1(matched1,1),spots1(matched1,2),'k.'); hold on;
    plot(spots2(matched2,1),spots2(matched2,2),'ro'); 
    plot(warped(:,1),warped(:,2),'bo'); 

    subplot(1,3,3); 
	Ncolor(3*cat(3,refImage,rotImageOut));
end