function [spotTable,pars] = ChrTracer_FitPSFs(fidSpots,dataSpots,varargin)
% [dTable,pars] = ChrTracer_FitPSFs(fidSpots,dataSpots,varargin)
% Loop over hybes, apply 3D Gaussian fit to all spots using common
% thresholds and fit parameters for fidSpots and dataSpots separately.
% Plot fit results in Figure 2 and Figure 4 on top of plots setup by
% ChrTracer_CropAndPlot.  
% Return a table of spot locations, also recording brightness and flag to
% distinguish data  from fiducial.
% 
global scratchPath;
defaults = cell(0,3);
defaults(end+1,:) = {'zstart', 'integer', 6}; 
defaults(end+1,:) = {'zstop', 'integer', 55}; 
defaults(end+1,:) = {'fov','integer',1};
defaults(end+1,:) = {'lociXY', 'array', [0,0]};
defaults(end+1,:) = {'fidMinPeakHeight', 'positive', 600};
defaults(end+1,:) = {'fidCameraBackground', 'nonnegative', 0};
defaults(end+1,:) = {'fidPeakBlur', 'nonnegative', .5};
defaults(end+1,:) = {'fidTroubleshoot', 'boolean', false};
defaults(end+1,:) = {'fidMaxFitWidth', 'positive', 6};
defaults(end+1,:) = {'fidMinSep', 'nonnegative', 3};
defaults(end+1,:) = {'datMinPeakHeight', 'positive', 2000};
defaults(end+1,:) = {'datCameraBackground', 'nonnegative', 500};
defaults(end+1,:) = {'datPeakBlur', 'nonnegative', .5};
defaults(end+1,:) = {'datTroubleshoot', 'boolean', false};
defaults(end+1,:) = {'datMaxFitWidth', 'positive', 6};
defaults(end+1,:) = {'datMinSep', 'nonnegative', 3};
defaults(end+1,:) = {'verbose', 'boolean', true};
defaults(end+1,:) = {'showPlots', 'boolean', true};
defaults(end+1,:) = {'saveData', 'boolean', false};
defaults(end+1,:) = {'saveFolder', 'string', scratchPath};
defaults(end+1,:) = {'labelOffsetsXY', 'array', []};
defaults(end+1,:) = {'labelOffsetsXZ', 'array', []};
defaults(end+1,:) = {'firstSpot', 'integer', 1};
defaults(end+1,:) = {'lastSpot', 'integer', []};
defaults(end+1,:) = {'figHandles','cell',{}};
defaults(end+1,:) = {'reloadFigs','boolean',false};
defaults(end+1,:) = {'numParallel','integer',1};
defaults(end+1,:) = {'nmXYpix','positive',154};
defaults(end+1,:) = {'nmZpix','positive',100};
defaults(end+1,:) = {'fidKeepBrightest','integer',2};
defaults(end+1,:) = {'fidRelativeHeight','fraction',0};
defaults(end+1,:) = {'fidMinHBratio','nonnegative',1.2}; % peak value over background value
defaults(end+1,:) = {'fidMinAHratio','nonnegative',.25}; % fitted height over background vs peak value
defaults(end+1,:) = {'fidMaxUncert','nonnegative',2}; % pixels
defaults(end+1,:) = {'datKeepBrightest','integer',2};
defaults(end+1,:) = {'datRelativeHeight','fraction',0};
defaults(end+1,:) = {'datMinHBratio','nonnegative',1.2}; % peak value over background value
defaults(end+1,:) = {'datMinAHratio','nonnegative',.25}; % fitted height over background vs peak value
defaults(end+1,:) = {'datMaxUncert','nonnegative',2}; % pixels

pars = ParseVariableArguments(varargin, defaults, mfilename);

if isempty(pars.labelOffsetsXZ) && pars.showPlots 
    warning('missing labelOffsetsXY needed for plotting');
    warning('fits will be computed and saved but not plotted');
    pars.showPlots = false; 
end

% turn off some annoying and non-relevant warnings
warning('off','MATLAB:prnRenderer:opengl');
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');
warning('off','MATLAB:singularMatrix');

%% Fit 3D Gaussians in matlab using FitPsf3D
[~,~,~,numHybes] = size(fidSpots{pars.firstSpot});
numSpots = length(fidSpots);
zs = pars.zstart:pars.zstop;
if isempty(pars.lastSpot)
    pars.lastSpot = numSpots;
end

if isempty(pars.figHandles)
    pars.figHandles = {1,2,3,4};
    pars.figHandles{2} = figure(pars.figHandles{2});
    pars.figHandles{4} = figure(pars.figHandles{4});
    pars.reloadFigs = true;
end



cellTable = cell(numSpots*numHybes,1);
%%
if pars.numParallel < 2
    lociXYx = pars.lociXY(:,1);
    lociXYy = pars.lociXY(:,2);
    tic
    k=0;
    for s=pars.firstSpot:pars.lastSpot  % another good place for a parfor loop
        for h=1:numHybes % h=3
            k=k+1;
            if pars.fidTroubleshoot; figure(11); clf; end
            dTable = FitPsf3D(squeeze(fidSpots{s}(:,:,zs,h)),...
                'minPeakHeight',pars.fidMinPeakHeight,...
                'peakBlur',pars.fidPeakBlur,...
                'cameraBackground',pars.fidCameraBackground,...
                'maxFitWidth',pars.fidMaxFitWidth,...
                'minSep',pars.fidMinSep,...
                'keepBrightest',pars.fidKeepBrightest,...
                'relativeHeight',pars.fidRelativeHeight,...
                'minHBratio',pars.fidMinHBratio,...
                'minAHratio',pars.fidMinAHratio,...
                'maxUncert',pars.fidMaxUncert,...
                'troubleshoot',pars.fidTroubleshoot);
            
            % convert units
            dTable.x = dTable.x*pars.nmXYpix;
            dTable.y = dTable.y*pars.nmXYpix;
            dTable.z = dTable.z*pars.nmZpix;
            dTable.xL = dTable.xL*pars.nmXYpix;
            dTable.xU = dTable.xU*pars.nmXYpix;
            dTable.yL = dTable.yL*pars.nmXYpix;
            dTable.yU = dTable.yU*pars.nmXYpix;
            dTable.zL = dTable.zL*pars.nmZpix;
            dTable.zU = dTable.zU*pars.nmZpix;
            dTable.wx = dTable.wx*pars.nmXYpix;
            dTable.wy = dTable.wy*pars.nmXYpix;
            dTable.wz = dTable.wz*pars.nmZpix;
            % keep track of hybe, spot number, and is fiduial. 
            dTable.hybe = h*ones(length(dTable.x),1);
            dTable.s = s*ones(length(dTable.x),1);
            dTable.isfid = true(length(dTable.x),1);
            dTable.fov = pars.fov*ones(length(dTable.x),1);
            dTable.locusX = pars.lociXY(s,1)*ones(length(dTable.x),1);
            dTable.locusY = pars.lociXY(s,2)*ones(length(dTable.x),1);
            cellTable{k} = dTable;
            if pars.fidTroubleshoot; pause; end
            if pars.showPlots
                % xy fig
                tileName = ['fov',num2str(pars.fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy3-5_tileFig'];
                tileSaveName = ['fov',num2str(pars.fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy3-5_tileFig_Fit'];
                if h==1
                    try
                        pars.figHandles{2}.Position = [581 60 1382 531];
                        allLines = findobj('Type','Line');
                        delete(allLines); 
                    catch 
                        pars.reloadFigs = true;  
                        pars.figHandles{2} = figure(2); 
                    end
                    if pars.reloadFigs
                        figure(pars.figHandles{2}); close; 
                        uiopen([pars.saveFolder,filesep,tileName,'.fig'],1);  
                        pars.figHandles{2} = gcf;
                    end
                end
                figure(pars.figHandles{2}); subplot(1,2,1); hold on; 
                plot(pars.labelOffsetsXY(h,1)+dTable.x/pars.nmXYpix,...
                    pars.labelOffsetsXY(h,2)+dTable.y/pars.nmXYpix,'yo','MarkerSize',20);

                % xz fig
                ztileName = ['fov',num2str(pars.fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy3-5_zTileFig'];
                ztileSaveName = ['fov',num2str(pars.fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy3-5_zTileFig_Fit'];
                if h==1 
                    try
                        pars.figHandles{4}.Position = [372 548 1565 425];
                        allLines = findobj('Type','Line');
                        delete(allLines); 
                    catch
                        pars.reloadFigs = true;
                        pars.figHandles{4} = figure(4); 
                    end
                    if pars.reloadFigs
                        figure(pars.figHandles{4}); close; 
                        uiopen([pars.saveFolder,filesep,ztileName,'.fig'],1);  
                        pars.figHandles{4} = gcf;
                    end
                end
                figure(pars.figHandles{4}); subplot(1,2,1); hold on; 
                plot(pars.labelOffsetsXZ(h,1)+dTable.x/pars.nmXYpix,...
                    pars.labelOffsetsXZ(h,2)+dTable.z/pars.nmZpix,'yo','MarkerSize',20);
            end

            k=k+1;
            if pars.datTroubleshoot; figure(11); clf; end
            dTable = FitPsf3D(squeeze(dataSpots{s}(:,:,zs,h)),...
                'minPeakHeight',pars.datMinPeakHeight,...
                'peakBlur',pars.datPeakBlur,...
                'cameraBackground',pars.datCameraBackground,...
                'maxFitWidth',pars.datMaxFitWidth,...
                'minSep',pars.datMinSep,...
                'keepBrightest',pars.datKeepBrightest,...
                'relativeHeight',pars.datRelativeHeight,...
                'minHBratio',pars.datMinHBratio,...
                'minAHratio',pars.datMinAHratio,...
                'maxUncert',pars.datMaxUncert,...
                'troubleshoot',pars.datTroubleshoot);
            % convert units
            dTable.x = dTable.x*pars.nmXYpix;
            dTable.y = dTable.y*pars.nmXYpix;
            dTable.z = dTable.z*pars.nmZpix;
            dTable.xL = dTable.xL*pars.nmXYpix;
            dTable.xU = dTable.xU*pars.nmXYpix;
            dTable.yL = dTable.yL*pars.nmXYpix;
            dTable.yU = dTable.yU*pars.nmXYpix;
            dTable.zL = dTable.zL*pars.nmZpix;
            dTable.zU = dTable.zU*pars.nmZpix;
            dTable.wx = dTable.wx*pars.nmXYpix;
            dTable.wy = dTable.wy*pars.nmXYpix;
            dTable.wz = dTable.wz*pars.nmZpix;
            % keep track of hybe, spot number, and is fiduial. 
            dTable.hybe = h*ones(length(dTable.x),1);
            dTable.s = s*ones(length(dTable.x),1);
            dTable.isfid = false(length(dTable.x),1);  
            dTable.fov = pars.fov*ones(length(dTable.x),1);
            dTable.locusX = pars.lociXY(s,1)*ones(length(dTable.x),1);
            dTable.locusY = pars.lociXY(s,2)*ones(length(dTable.x),1);
            cellTable{k} = dTable;
            if pars.showPlots
                figure(pars.figHandles{2}); subplot(1,2,2); hold on;    
                plot(pars.labelOffsetsXY(h,1)+dTable.x/pars.nmXYpix,...
                    pars.labelOffsetsXY(h,2)+dTable.y/pars.nmXYpix,'yo','MarkerSize',20);
                figure(pars.figHandles{4}); subplot(1,2,2); hold on; 
                plot(pars.labelOffsetsXZ(h,1)+dTable.x/pars.nmXYpix,...
                    pars.labelOffsetsXZ(h,2)+dTable.z/pars.nmZpix,'yo','MarkerSize',20);
            end
        end
        if pars.saveData
            SaveFigure(pars.figHandles{2},'name',tileSaveName,'formats',{'png','fig'},'overwrite',true,'saveData',pars.saveData);
            SaveFigure(pars.figHandles{4},'name',ztileSaveName,'formats',{'png','fig'},'overwrite',true,'saveData',pars.saveData);
        end
        if pars.verbose
            disp(['Fitting 3D gaussian spots: ' num2str(100*s/numSpots,3),' % complete...']);
        end

    end
    t= toc;
    spotTable = cat(1,cellTable{:});
    if pars.verbose
       disp(['FitPsf3D fit ',num2str(height(spotTable)),' spots in ',num2str(t/60),' minutes']);
    end
    
else % run parallel;
    if isempty(gcp('nocreate'))
        p = parpool(pars.numParallel); 
    else
        p = gcp;
    end
    
    % parfor can't parse structures
    firstSpot = pars.firstSpot;
    lastSpot = pars.lastSpot;
    fidMinPeakHeight = pars.fidMinPeakHeight;
    fidPeakBlur = pars.fidPeakBlur;
    fidCameraBackground = pars.fidCameraBackground;
    fidMaxFitWidth = pars.fidMaxFitWidth;
    fidMinSep = pars.fidMinSep;
    fidKeepBrightest = pars.fidKeepBrightest;
    fidRelativeHeight = pars.fidRelativeHeight;
    fidMinHBratio = pars.fidMinHBratio;
    fidMinAHratio = pars.fidMinAHratio;
    fidMaxUncert = pars.fidMaxUncert;
    fov = pars.fov; 
    lociXYx = pars.lociXY(:,1);
    lociXYy = pars.lociXY(:,2);
    showPlots = pars.showPlots;
    saveFolder = pars.saveFolder; 
    saveData = pars.saveData;
    labelOffsetsXY = pars.labelOffsetsXY;
    labelOffsetsXZ = pars.labelOffsetsXZ;
    datMinPeakHeight = pars.datMinPeakHeight;
    datPeakBlur = pars.datPeakBlur;
    datCameraBackground = pars.datCameraBackground;
    datMaxFitWidth = pars.datMaxFitWidth;
    datMinSep = pars.datMinSep;
    datKeepBrightest = pars.datKeepBrightest;
    datRelativeHeight = pars.datRelativeHeight;
    datMinHBratio = pars.datMinHBratio;
    datMinAHratio = pars.datMinAHratio;
    datMaxUncert = pars.datMaxUncert;
    verbose = pars.verbose;
    nmXYpix = pars.nmXYpix;
    nmZpix = pars.nmZpix;
    
    tic
    cellTable = {};
    parfor s=firstSpot:lastSpot  % another good place for a parfor loop
        fidTable = {}; 
        datTable = {};
        f2 = figure(2); 
        f4 = figure(4); 
        for h=1:numHybes % h=3
            dTable = FitPsf3D(squeeze(fidSpots{s}(:,:,zs,h)),...
                'minPeakHeight',fidMinPeakHeight,...
                'peakBlur',fidPeakBlur,...
                'cameraBackground',fidCameraBackground,...
                'maxFitWidth',fidMaxFitWidth,...
                'minSep',fidMinSep,...
                'keepBrightest',fidKeepBrightest,...
                'relativeHeight',fidRelativeHeight,...
                'minHBratio',fidMinHBratio,...
                'minAHratio',fidMinAHratio,...
                'maxUncert',fidMaxUncert,...
                'troubleshoot',false);
            % convert units
            dTable.x = dTable.x*nmXYpix;
            dTable.y = dTable.y*nmXYpix;
            dTable.z = dTable.z*nmZpix;
            dTable.xL = dTable.xL*nmXYpix;
            dTable.xU = dTable.xU*nmXYpix;
            dTable.yL = dTable.yL*nmXYpix;
            dTable.yU = dTable.yU*nmXYpix;
            dTable.zL = dTable.zL*nmZpix;
            dTable.zU = dTable.zU*nmZpix;
            dTable.wx = dTable.wx*nmXYpix;
            dTable.wy = dTable.wy*nmXYpix;
            dTable.wz = dTable.wz*nmZpix;
            % keep track of hybe, spot number, and is fiduial. 
            dTable.hybe = h*ones(length(dTable.x),1);
            dTable.s = s*ones(length(dTable.x),1);
            dTable.isfid = true(length(dTable.x),1);
            dTable.fov = fov*ones(length(dTable.x),1);
            dTable.locusX = lociXYx(s)*ones(length(dTable.x),1);
            dTable.locusY = lociXYy(s)*ones(length(dTable.x),1);
            fidTable{h} = dTable;
            if showPlots
                % xy fig
                tileName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy3-5_tileFig'];
                tileSaveName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy3-5_tileFig_Fit'];
                if h==1
                    figure(f2); close; 
                    uiopen([saveFolder,filesep,tileName,'.fig'],1);  
                    f2 = gcf;
                    f2.Position = [581 60 1382 531];
                    f2.Visible = 'off';
                end
                figure(f2); subplot(1,2,1); hold on; 
                plot(labelOffsetsXY(h,1)+dTable.x/nmXYpix,...
                    labelOffsetsXY(h,2)+dTable.y/nmXYpix,'yo','MarkerSize',20);

                % xz fig
                ztileName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy3-5_zTileFig'];
                ztileSaveName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_cy3-5_zTileFig_Fit'];
                if h==1 
                    figure(f4); close;
                    uiopen([saveFolder,filesep,ztileName,'.fig'],1);  
                    f4 = gcf;
                    f4.Position = [372 548 1565 425];
                    f4.Visible = 'off';
                end
                figure(f4); subplot(1,2,1); hold on; 
                plot(labelOffsetsXZ(h,1)+dTable.x/nmXYpix,...
                    labelOffsetsXZ(h,2)+dTable.z/nmZpix,'yo','MarkerSize',20);
            end

            dTable = FitPsf3D(squeeze(dataSpots{s}(:,:,zs,h)),...
                'minPeakHeight',datMinPeakHeight,...
                'peakBlur',datPeakBlur,...
                'cameraBackground',datCameraBackground,...
                'maxFitWidth',datMaxFitWidth,...
                'minSep',datMinSep,...
                'keepBrightest',datKeepBrightest,...
                'relativeHeight',datRelativeHeight,...
                'minHBratio',datMinHBratio,...
                'minAHratio',datMinAHratio,...
                'maxUncert',datMaxUncert,...
                'troubleshoot',false);
            % convert units
            dTable.x = dTable.x*nmXYpix;
            dTable.y = dTable.y*nmXYpix;
            dTable.z = dTable.z*nmZpix;
            dTable.xL = dTable.xL*nmXYpix;
            dTable.xU = dTable.xU*nmXYpix;
            dTable.yL = dTable.yL*nmXYpix;
            dTable.yU = dTable.yU*nmXYpix;
            dTable.zL = dTable.zL*nmZpix;
            dTable.zU = dTable.zU*nmZpix;
            dTable.wx = dTable.wx*nmXYpix;
            dTable.wy = dTable.wy*nmXYpix;
            dTable.wz = dTable.wz*nmZpix;
            % keep track of hybe, spot number, and is fiduial. 
            dTable.hybe = h*ones(length(dTable.x),1);
            dTable.s = s*ones(length(dTable.x),1);
            dTable.isfid = false(length(dTable.x),1);  
            dTable.fov = fov*ones(length(dTable.x),1);
            dTable.locusX = lociXYx(s)*ones(length(dTable.x),1);
            dTable.locusY = lociXYy(s)*ones(length(dTable.x),1);
            datTable{h} = dTable;
            if showPlots
                figure(f2); subplot(1,2,2); hold on;    
                plot(labelOffsetsXY(h,1)+dTable.x/nmXYpix,...
                    labelOffsetsXY(h,2)+dTable.y,'yo','MarkerSize',20);
                figure(f4); subplot(1,2,2); hold on; 
                plot(labelOffsetsXZ(h,1)+dTable.x/nmXYpix,...
                    labelOffsetsXZ(h,2)+dTable.z/nmZpix,'yo','MarkerSize',20);
            end
        end
        if saveData
            warning('off','MATLAB:prnRenderer:opengl');
            SaveFigure(f2,'name',tileSaveName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);
            SaveFigure(f4,'name',ztileSaveName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);
        end
        cellTable{s} = cat(1,datTable{:},fidTable{:}); 
    end
    t= toc;
    spotTable = cat(1,cellTable{:});
    if pars.verbose
       disp(['FitPsf3D fit ',num2str(height(spotTable)),' spots in ',num2str(t/60),' minutes']);
    end
end





