function [fidSpots,dataSpots,pars] = ChrTracer_CropSpots(fiducialFrames,dataFrames,varargin)
% pars = ChrTracer_CropAndPlot('parameters',pars);
% pars = ChrTracer_CropAndPlot('parameters',pars,'optionName',optionName);
% 
% -------------------------------------------------------------------------
% Required Inputs
% -------------------------------------------------------------------------
% 'parameters',pars
% A structure containing pars.fiducialFrames, and pars.dataFrames
% Both are 4D matrices, nRows x nCols x nStacks x nHybes 
% 
% -------------------------------------------------------------------------
% Optional inputs
% -------------------------------------------------------------------------
% 
% 
% -------------------------------------------------------------------------
% Outputs
% -------------------------------------------------------------------------
% 
% 
% -------------------------------------------------------------------------
% Notes
% -------------------------------------------------------------------------
% called by ChrTracer


defaults = cell(0,3);
% key parameters
defaults(end+1,:) = {'boxWidth', 'positive', 16};
defaults(end+1,:) = {'goodHybes','array',[]};
defaults(end+1,:) = {'showFolderNames', 'boolean',false};
defaults(end+1,:) = {'fidGain','positive',1};
defaults(end+1,:) = {'datGain','positive',1};
% FOV parameters
defaults(end+1,:) = {'fov', 'integer', 0};  % Field of view in experiment,  0 for not known;
defaults(end+1,:) = {'showPlots','boolean',true};
defaults(end+1,:) = {'saveData','boolean',false};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'firstSpot', 'integer', 1};
defaults(end+1,:) = {'lastSpot', 'integer', []};
defaults(end+1,:) = {'numParallel', 'integer',1};
pars = ParseVariableArguments(varargin, defaults, mfilename);


%% A little more parameter parsing
% break some variables out of 'pars', handle missing parameters 
tic;

% close open figure handles, remember where the user put them; 
allFigs = findobj('Type','Figure');
tileLocs = cell(4,1);
del_i = false;
for i=1:length(allFigs)
    isTile = strfind(allFigs(i).Name,'fid_tileXY');
    if ~isempty(isTile)
        tileLocs{1} = allFigs(i).Position;
        del_i = true;
    end
    isTile = strfind(allFigs(i).Name,'fid_tileXZ');
    if ~isempty(isTile)
        tileLocs{2} = allFigs(i).Position;
        del_i = true;
    end
    isTile = strfind(allFigs(i).Name,'dat_tileXY');
    if ~isempty(isTile)
        tileLocs{3} = allFigs(i).Position;
        del_i = true;
    end
    isTile = strfind(allFigs(i).Name,'dat_tileXZ');
    if ~isempty(isTile)
        tileLocs{4} = allFigs(i).Position;
        del_i = true;
    end
    if del_i
        delete(allFigs(i));
    end
    del_i = false;
end

[numHybes,numDataChns] = size(dataFrames); 
numSpots = size(pars.lociXY,1);

pars.goodHybes = logical(pars.goodHybes);
if isempty(pars.goodHybes)
    pars.goodHybes = true(numHybes,1);
else 
    numHybes = sum(pars.goodHybes);
end

if isempty(pars.lastSpot)
    pars.lastSpot = numSpots;
end


try
[tileLabels_fid,tileLabels_dat] = TileLabelsFromEtable(pars.eTable(pars.goodHybes,:));
catch
   error('problem'); 
end

%% Crop images

boxWidth = pars.boxWidth; % image width
spots = pars.lociXY;
fidSpots = cell(numSpots,1);
dataSpots = cell(numSpots,1);
[numRows,numCols,~] = size(fiducialFrames{1});
for s=pars.firstSpot:pars.lastSpot
    selectRows = max(1,spots(s,2)-boxWidth/2):min(spots(s,2)+boxWidth/2,numRows);
    selectCols = max(1,spots(s,1)-boxWidth/2):min(spots(s,1)+boxWidth/2,numCols);
    for h=find(pars.goodHybes)
        fidSpots{s}{h} = fiducialFrames{h}(selectRows,selectCols,:);
        for n=1:numDataChns
            dataSpots{s}{h,n} = dataFrames{h,n}(selectRows,selectCols,:);
        end
    end
end

%% Save images of individual spots as png and fig files
% a little more code to detect errors, supress unhelpful warnings. 
hasFid = ~isempty(fiducialFrames);
hasDat = ~isempty(dataFrames); 

if pars.showPlots
    vis = 'on';
else
    vis = 'off';
    warning('off','MATLAB:prnRenderer:opengl');
end

   
if pars.numParallel > 2
    vis = 'off'; % for parallel processing we don't want this on; 
    warning('off','MATLAB:prnRenderer:opengl');   
    warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');
    if isempty(gcp('nocreate'))
        parpool(pars.numParallel); 
    end
end    

% to avoid passing pars inside the parfor loop
saveData = pars.saveData;
saveFolder = pars.saveFolder;
firstSpot = pars.firstSpot;
lastSpot = pars.lastSpot;
lociXYx = pars.lociXY(:,1);
lociXYy = pars.lociXY(:,2);
fidGain = pars.fidGain;
datGain = pars.datGain;
fov = pars.fov;
numParallel = pars.numParallel;

disp('numParallel');
numParallel

if pars.showPlots
% Display images
if numParallel > 1
    parfor s=firstSpot:lastSpot  % par for
        warning('off','MATLAB:prnRenderer:opengl');
        SetFigureSavePath(saveFolder,'verbose',false); %

        if hasFid
            fSpot = cat(4,fidSpots{s}{:});
            im = squeeze(max(fSpot(:,:,:,:),[],3));
            pks_fid = squeeze(max(max(im))); 
            im = fidGain*im;
            pkFig = figure(1); clf; 
                pkFig.Visible = vis;
                subplot(1,2,1); bar(pks_fid); title('fid peaks');

            overlayFig = figure(2); clf; 
                overlayFig.Visible = vis;
                subplot(1,2,1);
                imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_fid_overlayFig'];
                overlayFig.Name = imName;
                imagesc(IncreaseContrast(Ncolor(im),'low',.1,'high',.999 )); title(imName,'interpreter','none');

            tileFigFid = figure(4); clf; 
                % tileFig.Position = [681 60 1382 531]; % (could become a parameter)
                tileFigFid.Visible = vis;
                imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_fid_tileXY'];
                tileFigFid.Name = imName;
                TileImageStack(im,'tileLabels',tileLabels_fid);  
                SaveFigure(tileFigFid,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);

            % zproj fid
            im = fidGain*max(permute(fSpot(:,:,:,:),[3,2,1,4]),[],3);
            zProjFig = figure(3); clf;
                zProjFig.Visible = vis;
                subplot(1,2,1);
                imName =['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_fid_zProjFig'];
                zProjFig.Name = imName;
                imagesc(IncreaseContrast(Ncolor(im),'low',.1,'high',.999 )); title(imName,'interpreter','none');
             zTileFigFid = figure(5);  clf;
                % zTileFigFid.Position = [472 648 1565 425];
                zTileFigFid.Visible = vis;
                imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_fid_tileXZ'];
                zTileFigFid.Name = imName;
                TileImageStack(im,'tileLabels',tileLabels_fid); 
                SaveFigure(zTileFigFid,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);
        end
        if hasDat
            dSpot = cat(4,dataSpots{s}{:});   
            % dat
            im = squeeze(max(dSpot(:,:,:,:),[],3));
            pks_dat = squeeze(max(max(im)));
            im = datGain*im;
            figure(pkFig); 
                subplot(1,2,2); bar(pks_dat); title('dat peaks');
                imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_peakFig'];
                SaveFigure(pkFig,'name',imName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);       
            figure(overlayFig);  
                subplot(1,2,2);
                imTitle = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_dat_overlayFig'];
                imagesc(IncreaseContrast(Ncolor(im),'low',.1,'high',.999 )); 
                title(imTitle,'interpreter','none');
                imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_dat_overlayFig'];
                SaveFigure(overlayFig,'name',imName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);       

            tileFigDat = figure(6);  clf;
                tileName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_dat_tileXY'];
                tileFigDat.Name = tileName;
                tileFigDat.Visible = vis;
                TileImageStack(im,'tileLabels',tileLabels_dat);
                SaveFigure(tileFigDat,'name',tileName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);

            % zproj dat
            im = datGain*max(permute(dSpot(:,:,:,:),[3,2,1,4]),[],3);
            figure(zProjFig);  
                subplot(1,2,2);
                imTitle =['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_dat_zProjFig'];
                imagesc(IncreaseContrast(Ncolor(im),'low',.1,'high',.999 )); 
                title(imTitle,'interpreter','none');
                imName =['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_dat_zProjFig'];
                SaveFigure(zProjFig,'name',imName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);
            zTileFigDat = figure(7); clf;
            ztileName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_dat_tileXZ'];    
                zTileFigDat.Name = ztileName;
                tileFigDat.Visible = vis;
                TileImageStack(im,'tileLabels',tileLabels_dat); 
                SaveFigure(zTileFigDat,'name',ztileName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);
        end
    end
    
    
else  %  non-parallel-processing
    
   
   for s=firstSpot:lastSpot  % par for
        warning('off','MATLAB:prnRenderer:opengl');
        SetFigureSavePath(saveFolder,'verbose',false); %

        if hasFid
            fSpot = cat(4,fidSpots{s}{:});
            im = squeeze(max(fSpot(:,:,:,:),[],3));
            pks_fid = squeeze(max(max(im))); 
            im = fidGain*im;
            pkFig = figure(1); clf; 
                pkFig.Visible = vis;
                subplot(1,2,1); bar(pks_fid); title('fid peaks');

            overlayFig = figure(2); clf; 
                overlayFig.Visible = vis;
                subplot(1,2,1);
                imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_fid_overlayFig'];
                overlayFig.Name = imName;
                imagesc(IncreaseContrast(Ncolor(im),'low',.1,'high',.999 )); title(imName,'interpreter','none');

            tileFigFid = figure(4); clf; 
                if ~isempty(tileLocs{1}); tileFigFid.Position = tileLocs{1}; end % location of last figure
                tileFigFid.Visible = vis;
                imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_fid_tileXY'];
                tileFigFid.Name = imName;
                TileImageStack(im,'tileLabels',tileLabels_fid);  
                SaveFigure(tileFigFid,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);

            % zproj fid
            im = fidGain*max(permute(fSpot(:,:,:,:),[3,2,1,4]),[],3);
            zProjFig = figure(3); clf;
                zProjFig.Visible = vis;
                subplot(1,2,1);
                imName =['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_fid_zProjFig'];
                zProjFig.Name = imName;
                imagesc(IncreaseContrast(Ncolor(im),'low',.1,'high',.999 )); title(imName,'interpreter','none');
             zTileFigFid = figure(5);  clf;
                if ~isempty(tileLocs{2}); zTileFigFid.Position = tileLocs{2}; end
                zTileFigFid.Visible = vis;
                imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_fid_tileXZ'];
                zTileFigFid.Name = imName;
                TileImageStack(im,'tileLabels',tileLabels_fid); 
                SaveFigure(zTileFigFid,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);
        end
        if hasDat
            dSpot = cat(4,dataSpots{s}{:});   
            % dat
            im = squeeze(max(dSpot(:,:,:,:),[],3));
            pks_dat = squeeze(max(max(im)));
            im = datGain*im;
            figure(pkFig); 
                subplot(1,2,2); bar(pks_dat); title('dat peaks');
                imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_peakFig'];
                SaveFigure(pkFig,'name',imName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);       
            figure(overlayFig);  
                subplot(1,2,2);
                imTitle = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_dat_overlayFig'];
                imagesc(IncreaseContrast(Ncolor(im),'low',.1,'high',.999 )); 
                title(imTitle,'interpreter','none');
                imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_dat_overlayFig'];
                SaveFigure(overlayFig,'name',imName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);       

            tileFigDat = figure(6);  clf;
                tileName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_dat_tileXY'];
                tileFigDat.Name = tileName;
                if ~isempty(tileLocs{3}); tileFigDat.Position = tileLocs{3}; end
                tileFigDat.Visible = vis;
                TileImageStack(im,'tileLabels',tileLabels_dat);
                SaveFigure(tileFigDat,'name',tileName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);

            % zproj dat
            im = datGain*max(permute(dSpot(:,:,:,:),[3,2,1,4]),[],3);
            figure(zProjFig);  
                subplot(1,2,2);
                imTitle =['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_dat_zProjFig'];
                imagesc(IncreaseContrast(Ncolor(im),'low',.1,'high',.999 )); 
                title(imTitle,'interpreter','none');
                imName =['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_dat_zProjFig'];
                SaveFigure(zProjFig,'name',imName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);
            zTileFigDat = figure(7); clf;
            ztileName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_dat_tileXZ'];    
                zTileFigDat.Name = ztileName;
                if ~isempty(tileLocs{4}); zTileFigDat.Position = tileLocs{4}; end;
                tileFigDat.Visible = vis;
                TileImageStack(im,'tileLabels',tileLabels_dat); 
                SaveFigure(zTileFigDat,'name',ztileName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);
        end
   end
end


t = toc;
if pars.verbose
    disp(['finished crop and plot in t=',num2str(t/60),'min.']);
end

end