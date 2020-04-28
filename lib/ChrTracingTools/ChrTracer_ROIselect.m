function [spots,im,pars] = ChrTracer_ROIselect(fiducialFrames,varargin)
% fiducialFrames is a cell of H hybes, each containing a 3-D image.

% defaults
defaults = cell(0,3);
% fov pars
defaults(end+1,:) = {'zstart','positive',6};
defaults(end+1,:) = {'zstop','positive',55};
defaults(end+1,:) = {'goodHybes','array',[]};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'fov','integer',1};
% ROI select pars
defaults(end+1,:) = {'ROIselectMode',{'manual','auto'},'manual'};
defaults(end+1,:) = {'autoSelectThreshold','fraction',.992};
defaults(end+1,:) = {'autoSelectDownsample','positive',3};
defaults(end+1,:) = {'mergeFrames','string','median'};
defaults(end+1,:) = {'border','nonnegative',10};
defaults(end+1,:) = {'coarseBackgroundScale','nonnegative',0}; % 50
defaults(end+1,:) = {'fineBackgroundScale','nonnegative',0};   % 5
defaults(end+1,:) = {'gain','postive',1};
defaults(end+1,:) = {'displayContrastLow','fraction',0};
defaults(end+1,:) = {'displayContrastHigh','fraction',.9999};
% 
defaults(end+1,:) = {'im','freeType',[]};
defaults(end+1,:) = {'badSpots','freeType',[]};
pars = ParseVariableArguments(varargin,defaults,mfilename);

if strcmp(pars.ROIselectMode,'manual')
    [spots,im] = ChrTracer_ManualSelectLocus(fiducialFrames,'parameters',pars);
else
    [spots,im] = ChrTracer_AutoSelectLoci(fiducialFrames,'parameters',pars);
end