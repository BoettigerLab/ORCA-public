function [overlays,success] = SaveFidAlignFrames(fiducialAlignFrames,varargin)
 
global figureSavePath;

defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true}; 
defaults(end+1,:) = {'fov', 'integer', 1}; 
defaults(end+1,:) = {'saveFigFile', 'boolean', true}; 
defaults(end+1,:) = {'saveFolder', 'string', figureSavePath}; 
defaults(end+1,:) = {'saveFormats', 'cell',{'fig','png'}}; %  {'fig','png'}, {}
defaults(end+1,:) = {'overwrite', 'boolean', false}; 
defaults(end+1,:) = {'contrastLow','fraction',.5};
defaults(end+1,:) = {'contrastHigh','fraction',.9995};
pars = ParseVariableArguments(varargin,defaults,mfilename);

f = pars.fov;
overlays = [];
figName = ['fov',num2str(f,'%03d'),'_fid_overlayFig'];
% if file isn't already written and there's no data to write it, SKIP
if ~exist([pars.saveFolder,figName,'.fig'],'file') && isempty(fiducialAlignFrames)
    success = false;
    if pars.verbose
        disp(['no aligned data preview found for FOV ',num2str(f),' skipping this FOV']);
    end
% if there is data to write, and the file doesn't exist or overwrite is on:
elseif ~exist([pars.saveFolder,figName,'.fig'],'file') || pars.overwrite
    im = cat(3,fiducialAlignFrames{:}); 
    nChns = size(im,3); 
    im = 2/nChns*IncreaseContrast(im,'low',pars.contrastLow,'high',pars.contrastHigh);
    imO = Ncolor(im); 
    overlays = imO;
    overlayFig = figure(1); clf; 
    imagesc(overlays);
    overlayFig.Name = ['fov',num2str(f,'%03d'),'_fid_overlayFig'];
    SaveFigure(overlayFig,'formats',pars.saveFormats,'overwrite',true,'saveData',pars.saveFigFile,'verbose',pars.verbose);
    success = true;
elseif exist([pars.saveFolder,figName,'.fig'],'file')
    success = true;
    if pars.verbose
       disp(['found existing ', figName, '.']);
    end
end
%----