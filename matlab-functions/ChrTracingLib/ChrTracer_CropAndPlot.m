function [fidSpots,dataSpots,pars] = ChrTracer_CropAndPlot(fiducialFrames,dataFrames,varargin)
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
% maybe output should be ctChains, a lighter-weight data structure where we
% drop the full multi-gig data of the aligned images. 


defaults = cell(0,3);
% key parameters
defaults(end+1,:) = {'boxWidth', 'positive', 16};
defaults(end+1,:) = {'goodHybes','array',[]};
defaults(end+1,:) = {'showFolderNames', 'boolean',false};
defaults(end+1,:) = {'fidGain','positive',1};
defaults(end+1,:) = {'datGain','positive',1};
% FOV parameters
defaults(end+1,:) = {'zstart', 'integer', 6}; 
defaults(end+1,:) = {'zstop', 'integer', 55}; 
defaults(end+1,:) = {'fov', 'integer', 0};  % Field of view in experiment,  0 for not known;
defaults(end+1,:) = {'showPlots','boolean',true};
defaults(end+1,:) = {'saveData','boolean',false};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'firstSpot', 'integer', 1};
defaults(end+1,:) = {'lastSpot', 'integer', []};
defaults(end+1,:) = {'numParallel', 'integer',1};
pars = ParseVariableArguments(varargin, defaults, mfilename);


%% Crop 
tic;
boxWidth = pars.boxWidth; % image width
spots = pars.lociXY;
fov = pars.fov;
[numRows,numCols,~,numHybes] = size(fiducialFrames);
numSpots = size(spots,1);
zs = pars.zstart:pars.zstop;
if isempty(pars.lastSpot)
    pars.lastSpot = numSpots;
end


if isfield(pars,'eTable') && pars.showFolderNames
    tileLabels = pars.eTable.FolderName;
else
    tileLabels = cellstr(num2str((1:numHybes)'));
end

hasCy3 = ~isempty(fiducialFrames);
hasCy5 = ~isempty(dataFrames); 

fidSpots = cell(numSpots,1);
dataSpots = cell(numSpots,1);
for s=pars.firstSpot:pars.lastSpot
    selectRows = max(1,spots(s,2)-boxWidth/2):min(spots(s,2)+boxWidth/2,numRows);
    selectCols = max(1,spots(s,1)-boxWidth/2):min(spots(s,1)+boxWidth/2,numCols);
    fidSpots{s} = fiducialFrames(selectRows,selectCols,:,:);
    dataSpots{s} = dataFrames(selectRows,selectCols,:,:);
end


%% Save images of individual spots as png and fig files
if pars.showPlots
    vis = 'on';
else
    vis = 'off';
    warning('off','MATLAB:prnRenderer:opengl');
end
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');

saveData = pars.saveData;


   
if pars.numParallel < 2
    lociXYx = pars.lociXY(:,1);
    lociXYy = pars.lociXY(:,2);
    for s=pars.firstSpot:pars.lastSpot
        % SetFigureSavePath(saveFolder,'verbose',false); % needs testing. last time these got dumped into the working directory
        % cy3 
        if hasCy3
            im = pars.fidGain*squeeze(max(fidSpots{s}(:,:,zs,:),[],3));
            try pks_cy3 = max(reshape(im,(boxWidth+1)^2,numHybes)); 
            catch er
                disp(er.message)
                pks_cy3 = 0;
            end
            pkFig = figure(5); 
                pkFig.Visible = vis;
                clf; subplot(1,2,1); 
                bar(pks_cy3); title('cy3 peaks');
            overlayFig = figure(1); 
                overlayFig.Visible = vis;
                clf; subplot(1,2,1);
                imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy3_overlayFig'];
                imagesc(IncreaseContrast(Ncolor(im),'low',.1,'high',.999 )); 
                title(imName,'interpreter','none');
            tileFig = figure(2); clf;   
                tileFig.Position = [681 60 1382 531]; % (could become a parameter)
                tileFig.Visible = vis;
                subplot(1,2,1);
                imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy3_tileFig'];
                TileImage(IncreaseContrast(im,'low',.1,'high',.9999),'tileLabels',tileLabels);
                title(imName,'interpreter','none');
            % zproj cy3
            im = pars.fidGain*max(permute(fidSpots{s}(:,:,zs,:),[3,2,1,4]),[],3);
            zProjFig = figure(3); 
                zProjFig.Visible = vis;
                clf;  subplot(1,2,1);
                imName =['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy3_zProjFig'];
                imagesc(IncreaseContrast(Ncolor(im),'low',.1,'high',.999 )); title(imName,'interpreter','none');
            zTileFig = figure(4); 
                zTileFig.Visible = vis;
                clf;  subplot(1,2,1); 
                zTileFig.Position = [472 648 1565 425];
                imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy3_zTileFig'];
                TileImage(IncreaseContrast(im,'low',.1,'high',.9999),'tileLabels',tileLabels); 
                title(imName,'interpreter','none');
        end
        if hasCy5
            % cy5
            im = pars.datGain*squeeze(max(dataSpots{s}(:,:,zs,:),[],3));
            try pks_cy5 = max(reshape(im,(boxWidth+1)^2,numHybes));
            catch er
                disp(er.message);
                pks_cy5 = 0; 
            end
            figure(pkFig); 
                subplot(1,2,2); bar(pks_cy5); title('cy5 peaks');
                imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_peakFig'];
                SaveFigure(pkFig,'name',imName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);       
            figure(overlayFig);  
                subplot(1,2,2);
                imTitle = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy5_overlayFig'];
                imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy3-5_overlayFig'];
                imagesc(IncreaseContrast(Ncolor(im),'low',.1,'high',.999 )); 
                title(imTitle,'interpreter','none');
                SaveFigure(overlayFig,'name',imName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);       
            figure(tileFig);  
                subplot(1,2,2);
                tileTitle = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy5_tileFig'];
                tileName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy3-5_tileFig'];
                [~,pars.labelOffsetsXY] = TileImage(IncreaseContrast(im,'low',.1,'high',.9999),'tileLabels',tileLabels,'showImage',true); 
                title(tileTitle,'interpreter','none');
                SaveFigure(tileFig,'name',tileName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);
            % zproj cy5
            im = pars.datGain*max(permute(dataSpots{s}(:,:,zs,:),[3,2,1,4]),[],3);
            figure(zProjFig);  subplot(1,2,2);
                imTitle =['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy5_zProjFig'];
                imName =['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy3-5_zProjFig'];
                imagesc(IncreaseContrast(Ncolor(im),'low',.1,'high',.999 )); 
                title(imTitle,'interpreter','none');
                SaveFigure(zProjFig,'name',imName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);
            figure(zTileFig); subplot(1,2,2);
                ztileTitle = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy5_zTileFig'];
                ztileName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy3-5_zTileFig'];
                [~,pars.labelOffsetsXZ] = TileImage(IncreaseContrast(im,'low',.1,'high',.9999),'tileLabels',tileLabels,'showImage',true); 
                title(ztileTitle,'interpreter','none');
                SaveFigure(zTileFig,'name',ztileName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);
        end
        pars.figHandles = {overlayFig,tileFig,zProjFig,zTileFig};
    end

else
    
    if isempty(gcp('nocreate'))
        parpool(pars.numParallel); 
    end
    
    saveData = pars.saveData;
    saveFolder = pars.saveFolder;
    vis = 'off'; % for parallel processing we don't want this on; 
    warning('off','MATLAB:prnRenderer:opengl');
    firstSpot = pars.firstSpot;
    lastSpot = pars.lastSpot;
    lociXYx = pars.lociXY(:,1);
    lociXYy = pars.lociXY(:,2);
    labelOffsetsXY = {};
    labelOffsetsXZ = {};
    fidGain = pars.fidGain;
    datGain = pars.datGain;
    parfor s=firstSpot:lastSpot
        warning('off','MATLAB:prnRenderer:opengl');
        SetFigureSavePath(saveFolder,'verbose',false); % needs testing. last time these got dumped into the working directory
        if hasCy3
            fSpot = fidSpots{s};
            im = fidGain*squeeze(max(fSpot(:,:,zs,:),[],3));
            try pks_cy3 = max(reshape(im,(boxWidth+1)^2,numHybes)); 
            catch er
                disp(er.message)
                pks_cy3 = 0;
            end
            pkFig = figure(5); 
                pkFig.Visible = vis;
                clf; subplot(1,2,1); bar(pks_cy3); title('cy3 peaks');
            overlayFig = figure(1); 
                overlayFig.Visible = vis;
                clf; subplot(1,2,1);
                imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy3_overlayFig'];
                imagesc(IncreaseContrast(Ncolor(im),'low',.1,'high',.999 )); title(imName,'interpreter','none');
             tileFig = figure(2); 
                tileFig.Position = [681 60 1382 531]; % (could become a parameter)
                tileFig.Visible = vis;
                clf;  subplot(1,2,1); 
                imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy3_tileFig'];
                TileImage(IncreaseContrast(im,'low',.1,'high',.9999),'tileLabels',tileLabels); title(imName,'interpreter','none');
            % zproj cy3
            im = fidGain*max(permute(fSpot(:,:,zs,:),[3,2,1,4]),[],3);
            zProjFig = figure(3); 
                zProjFig.Visible = vis;
                clf;  subplot(1,2,1);
                imName =['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy3_zProjFig'];
                imagesc(IncreaseContrast(Ncolor(im),'low',.1,'high',.999 )); title(imName,'interpreter','none');
             zTileFig = figure(4); 
                zTileFig.Position = [472 648 1565 425];
                zTileFig.Visible = vis;
                clf;  subplot(1,2,1); 
                imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy3_zTileFig'];
                TileImage(IncreaseContrast(im,'low',.1,'high',.9999),'tileLabels',tileLabels); title(imName,'interpreter','none');
        end
        if hasCy5
            dSpot = dataSpots{s};   
            % cy5
            im = datGain*squeeze(max(dSpot(:,:,zs,:),[],3));
            try pks_cy5 = max(reshape(im,(boxWidth+1)^2,numHybes));
            catch er
                pks_cy5 = 0;
                disp(er.message);
            end
            figure(pkFig); 
                subplot(1,2,2); bar(pks_cy5); title('cy5 peaks');
                imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_peakFig'];
                SaveFigure(pkFig,'name',imName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);       
            figure(overlayFig);  
                subplot(1,2,2);
                imTitle = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy5_overlayFig'];
                imagesc(IncreaseContrast(Ncolor(im),'low',.1,'high',.999 )); 
                title(imTitle,'interpreter','none');
                imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy3-5_overlayFig'];
                SaveFigure(overlayFig,'name',imName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);       
            figure(tileFig);  
                subplot(1,2,2);
                tileName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy3-5_tileFig'];
                [~,labelOffsetsXY{s}] = TileImage(IncreaseContrast(im,'low',.1,'high',.9999),'tileLabels',tileLabels,'showImage',true); title(imName,'interpreter','none');
                SaveFigure(tileFig,'name',tileName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);
            % zproj cy5
            im = datGain*max(permute(dSpot(:,:,zs,:),[3,2,1,4]),[],3);
            figure(zProjFig);  
                subplot(1,2,2);
                imTitle =['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy5_zProjFig'];
                imagesc(IncreaseContrast(Ncolor(im),'low',.1,'high',.999 )); 
                title(imTitle,'interpreter','none');
                imName =['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy3-5_zProjFig'];
                SaveFigure(zProjFig,'name',imName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);
            figure(zTileFig); 
                subplot(1,2,2);
                ztileTitle = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy5_zTileFig'];
                [~,labelOffsetsXZ{s}] = TileImage(IncreaseContrast(im,'low',.1,'high',.9999),'tileLabels',tileLabels,'showImage',true);
                title(ztileTitle,'interpreter','none');
                ztileName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy3-5_zTileFig'];
                SaveFigure(zTileFig,'name',ztileName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);
        end
    end
    if pars.numParallel > 1
       pars.labelOffsetsXY =  labelOffsetsXY{lastSpot};
       pars.labelOffsetsXZ =  labelOffsetsXZ{lastSpot};
    end
end

t = toc;
if pars.verbose
    disp(['finished crop and plot in t=',num2str(t/60),'min.']);
end
