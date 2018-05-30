function [spots,im,pars] = ChrTracer_ManualSelectLocus(fiducialFrames,varargin)
% fiducialFrames is a cell of H hybes, each containing a 3-D image.
% b

% defaults
defaults = cell(0,3);
defaults(end+1,:) = {'displayContrastLow','fraction',0};
defaults(end+1,:) = {'displayContrastHigh','fraction',.999};
defaults(end+1,:) = {'gain','nonnegative',1};
defaults(end+1,:) = {'goodHybes','array',[]};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'fov','integer',1};
defaults(end+1,:) = {'lociXY','freeType',[]};
defaults(end+1,:) = {'im','freeType',[]};
defaults(end+1,:) = {'badSpots','freeType',[]};
pars = ParseVariableArguments(varargin,defaults,mfilename);


% auto complete a few parameters if missing
numHybes = length(fiducialFrames);
if isempty(pars.goodHybes)
    pars.goodHybes = true(1,numHybes);
else
    pars.goodHybes = logical(pars.goodHybes);
end


if isempty(pars.im)
    if iscell(fiducialFrames)
        im = cat(4,fiducialFrames{pars.goodHybes});
        im = pars.gain*squeeze(max(im,[],3));
    else
       im = fiducialFrames; 
    end
else
    im = pars.im;
end


if pars.showPlots
    nChns = size(im,3);
    imRGB = Ncolor(2*pars.gain/nChns*IncreaseContrast(im,'low',pars.displayContrastLow,'high',pars.displayContrastHigh));
    overlayFig = figure(1); clf; 
    imagesc(imRGB); hold on;
end

if ~isempty(pars.lociXY)
    plot(pars.lociXY(:,1),pars.lociXY(:,2),'yo');
    spotNums = cellstr(num2str((1:size(pars.lociXY,1))'));
    text(pars.lociXY(:,1)+1,pars.lociXY(:,2),spotNums,'color','w');
end
if ~isempty(pars.badSpots)
    plot(pars.badSpots(:,1),pars.badSpots(:,2),'ro');
    spotNums = cellstr(num2str((1:size(pars.badSpots,1))'));
    text(pars.badSpots(:,1)+1,pars.badSpots(:,2),spotNums,'color','w');
end

% prompt user to select locus of interest
dcmObj = datacursormode(overlayFig);
set(dcmObj,'DisplayStyle','datatip','SnapToDataVertex','off','Enable','on')
disp('Select a spot with the Data Cursor tool, then press Return.')
pause  % Wait while the user does this.hand
c_info = getCursorInfo(dcmObj);
if pars.verbose
    disp('Position captured'); 
end
spots = c_info.Position;