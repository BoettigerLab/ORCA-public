function [spots,im,pars] = ChrTracer_AutoSelectLoci(fiducialFrames,varargin)
% [spots,im,pars] = ChrTracer_AutoSelectLoci(fiducialFrames)
% fiducialFrames is a cell of H hybes, each containing a 3-D image.

% defaults
defaults = cell(0,3);
defaults(end+1,:) = {'mergeFrames','string','median'};
defaults(end+1,:) = {'goodHybes','array',[]};
defaults(end+1,:) = {'showPlots','boolean',true};
defaults(end+1,:) = {'saveData','boolean',true};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'fov','integer',1};
defaults(end+1,:) = {'gain','positive',1};
defaults(end+1,:) = {'im','nonnegative',[]};
defaults(end+1,:) = {'displayContrastLow','fraction',0};
defaults(end+1,:) = {'displayContrastHigh','fraction',.999};
% parameters for AutoSelectSpots
defaults(end+1,:) = {'autoSelectThreshold','fraction',.992};
defaults(end+1,:) = {'autoSelectDownsample','positive',3};
defaults(end+1,:) = {'border','nonnegative',10};
defaults(end+1,:) = {'coarseBackgroundScale','nonnegative',0}; % 50
defaults(end+1,:) = {'fineBackgroundScale','nonnegative',0};   % 5
defaults(end+1,:) = {'showPlots','boolean',true};
defaults(end+1,:) = {'numberSpots','boolean',true};
defaults(end+1,:) = {'saveData','boolean',false};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'gain','positive',1};
defaults(end+1,:) = {'fov','nonnegative',1};
defaults(end+1,:) = {'laplaceFilter','boolean',true};
defaults(end+1,:) = {'edgeThresh','nonnegative',0}; %.0001
defaults(end+1,:) = {'edgeSigma','nonnegative',1.2};
defaults(end+1,:) = {'overlayIm','array',[]};
defaults(end+1,:) = {'removeEdgeStack','nonnegative',20};

% parse variable arguments
pars = ParseVariableArguments(varargin, defaults, mfilename);
pars.showPlots = true; % hack

% auto complete a few parameters if missing
numHybes = length(fiducialFrames);
if isempty(pars.goodHybes)
    pars.goodHybes = true(1,numHybes);
else
    pars.goodHybes = logical(pars.goodHybes);
end

% just keeping backward compatibility here
%    in future release we will clean this up and break back-compatability
if isempty(pars.im)
    if iscell(fiducialFrames)
        try
            im = cat(4,fiducialFrames{pars.goodHybes});
            im = squeeze(max(im,[],3));
        catch
            im = cat(3,fiducialFrames{:});
        end
    else
        im = fiducialFrames;
    end
else
    im = pars.im;
end


if pars.showPlots   
    nChns = size(im,3);
    if nChns > 1
    imRGB = Ncolor(2*pars.gain/nChns*IncreaseContrast(im,'low',pars.displayContrastLow,'high',pars.displayContrastHigh));
    else
        imRGB = im; colormap(gray);
    end
    overlayFig = figure(1); clf; 
    imagesc(imRGB); hold on;
end

% Auto Select ROI  
switch pars.mergeFrames
    case 'median'
        imAve1 = nanmedian(im,3); 
    case 'average'
        imAve1 = nanmean(im,3); 
    case 'max'
        imAve1 = max(im,[],3);
end
% figure(10); imagesc(imAve1); 

[spots,pars] = AutoSelectSpots(imAve1,'parameters',pars,'overlayIm',imRGB);

