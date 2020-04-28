function im = ChrTracer_FOVsummaryPlots(fiducialAlignFrames,dataAlignFrames,varargin)

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
% key parameters
defaults(end+1,:) = {'verbose', 'boolean', true}; 
defaults(end+1,:) = {'showPlots', 'boolean', true}; 
defaults(end+1,:) = {'saveData', 'boolean', false}; 
defaults(end+1,:) = {'fov', 'integer', 1}; 
defaults(end+1,:) = {'stopOnError','boolean',false};
defaults(end+1,:) = {'eTable','freeType',[]};
defaults(end+1,:) = {'contrastLow','fraction',0};
defaults(end+1,:) = {'contrastHigh','fraction',1};
defaults(end+1,:) = {'fidGain','nonnegative',1};
defaults(end+1,:) = {'datGain','nonnegative',1}; 
defaults(end+1,:) = {'saveFig','boolean',true};
defaults(end+1,:) = {'showFidXY','boolean',true};
defaults(end+1,:) = {'showFidXZ','boolean',false};
defaults(end+1,:) = {'showDatXY','boolean',false};
defaults(end+1,:) = {'showDatXZ','boolean',false};
defaults(end+1,:) = {'rescaleTile','positive',.1};
defaults(end+1,:) = {'colormap','colormap','jet'};
pars = ParseVariableArguments(varargin,defaults,mfilename);

%% Main function

[tileLabels_fid,tileLabels_dat] = TileLabelsFromEtable(pars.eTable);

if pars.saveFig
    saveFormats = {'png','fig'};
else
    saveFormats = {'png'};
end

% Plot and save data
if pars.showFidXY
    fids = cat(4,fiducialAlignFrames{:});
    im =  pars.fidGain*squeeze(max(fids,[],3));
    imO = IncreaseContrast(Ncolor(im)); 
    overlayFig = figure(1); clf; imagesc(imO);
    overlayFig.Name = ['fov',num2str(pars.fov,'%03d'),'_fid_overlayFig'];
    SaveFigure(overlayFig,'formats',saveFormats,'overwrite',true,'saveData',pars.saveData);

    tileFig = figure(2); clf; 
    TileImageStack(imresize(im,pars.rescaleTile),'tileLabels',tileLabels_fid,'colormap',pars.colormap);
    tileFig.Name = ['fov',num2str(pars.fov,'%03d'),'_fid_tileFig'];
    SaveFigure(tileFig,'formats',saveFormats,'overwrite',true,'saveData',pars.saveData);
end

if pars.showFidXZ
    fids = cat(4,fiducialAlignFrames{:});
    zProj = pars.fidGain*squeeze(max(fids,[],2));
    zProjFig = figure(3); clf; Ncolor(permute(zProj,[2,1,3]));
    zProjFig.Name = ['fov',num2str(pars.fov,'%03d'),'_fid_zProjFig'];
    SaveFigure(zProjFig,'formats',saveFormats,'overwrite',true,'saveData',pars.saveData);

    imToTile = imresize(zProj,pars.rescaleTile);
    imToTile = permute(imToTile,[2,1,3]);
    zTileFig = figure(4); clf; TileImageStack(imToTile,'tileLabels',tileLabels_fid,'colormap',pars.colormap);
    zTileFig.Name = ['fov',num2str(pars.fov,'%03d'),'_fid_zTileFig'];
    SaveFigure(zTileFig,'formats',saveFormats,'overwrite',true,'saveData',pars.saveData);
end

if pars.showDatXY
    dats = cat(4,dataAlignFrames{:});
    im =  pars.datGain*squeeze(max(dats,[],3));
    imO = IncreaseContrast( Ncolor(im) );
    overlayFig = figure(5); clf; imagesc(imO);
    overlayFig.Name = ['fov',num2str(pars.fov,'%03d'),'_data_overlayFig'];
    SaveFigure(overlayFig,'formats',saveFormats,'overwrite',true,'saveData',pars.saveData);

    tileFig = figure(6); clf; 
    TileImageStack(imresize(im,pars.rescaleTile),'tileLabels',tileLabels_dat,'colormap',pars.colormap);
    tileFig.Name = ['fov',num2str(pars.fov,'%03d'),'_data_tileFig'];
    SaveFigure(tileFig,'formats',saveFormats,'overwrite',true,'saveData',pars.saveData);
end

if pars.showDatXZ
    dats = cat(4,dataAlignFrames{:});
    zProj = pars.datGain*squeeze(max(dats,[],2));
    zProjFig = figure(7); clf; Ncolor(permute(zProj,[2,1,3]));
    zProjFig.Name = ['fov',num2str(pars.fov,'%03d'),'_data_zProjFig'];
    SaveFigure(zProjFig,'formats',saveFormats,'overwrite',true,'saveData',pars.saveData);

    imToTile = imresize(zProj,pars.rescaleTile);
    imToTile = permute(imToTile,[2,1,3]);
    zTileFig = figure(8); clf; TileImageStack(imToTile,'tileLabels',tileLabels_dat,'colormap',pars.colormap);
    zTileFig.Name = ['fov',num2str(pars.fov,'%03d'),'_data_zTileFig'];
    SaveFigure(zTileFig,'formats',saveFormats,'overwrite',true,'saveData',pars.saveData);
end
