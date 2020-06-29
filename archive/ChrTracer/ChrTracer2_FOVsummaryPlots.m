function im = ChrTracer2_FOVsummaryPlots(fiducialAlignFrames,varargin)

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
% fov parameters
defaults(end+1,:) = {'verbose', 'boolean', true}; 
defaults(end+1,:) = {'showPlots', 'boolean', true}; 
defaults(end+1,:) = {'saveData', 'boolean', false}; 
defaults(end+1,:) = {'fov', 'integer', 1}; 
defaults(end+1,:) = {'stopOnError','boolean',false};
defaults(end+1,:) = {'eTable','freeType',[]};
% summary plot specific parameters
defaults(end+1,:) = {'contrastLow','fraction',0};
defaults(end+1,:) = {'contrastHigh','fraction',.999};
defaults(end+1,:) = {'gain','nonnegative',1};
defaults(end+1,:) = {'showFidXY','boolean',true};
defaults(end+1,:) = {'rescaleTile','positive',.1};
defaults(end+1,:) = {'colormap','colormap','jet'};
defaults(end+1,:) = {'saveFigFile', 'boolean', true}; 
pars = ParseVariableArguments(varargin,defaults,mfilename);

%% Main function

[tileLabels_fid,tileLabels_dat] = TileLabelsFromEtable(pars.eTable);

if pars.saveFigFile
    saveFormats = {'png','fig'};
else
    saveFormats = {'png'};
end

% Plot and save data
if pars.showFidXY
    fids = cat(3,fiducialAlignFrames{:}); 
    im =  pars.gain*fids; % figure(100); clf; imagesc(im(:,:,1));
    nChns = size(im,3); 
    imO = Ncolor(2/nChns*IncreaseContrast(im,'low',pars.contrastLow,'high',pars.contrastHigh)); 
    overlayFig = figure(1); clf; imagesc(imO);
    overlayFig.Name = ['fov',num2str(pars.fov,'%03d'),'_fid_overlayFig'];
    SaveFigure(overlayFig,'formats',saveFormats,'overwrite',true,'saveData',pars.saveData);

    tileFig = figure(2); clf; 
    TileImageStack(imresize(im,pars.rescaleTile),'tileLabels',tileLabels_fid,'colormap',pars.colormap);
    tileFig.Name = ['fov',num2str(pars.fov,'%03d'),'_fid_tileFig'];
    SaveFigure(tileFig,'formats',saveFormats,'overwrite',true,'saveData',pars.saveData);
end



