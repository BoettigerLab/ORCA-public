function [spotTable,pars] = ChrTracer_FitSpots(fidSpots,dataSpots,varargin)
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
defaults(end+1,:) = {'fov','integer',1};
defaults(end+1,:) = {'lociXY', 'array', [0,0]};
defaults(end+1,:) = {'verbose', 'boolean', true};
defaults(end+1,:) = {'showPlots', 'boolean', true};
defaults(end+1,:) = {'saveData', 'boolean', false};
defaults(end+1,:) = {'saveFolder', 'string', scratchPath};
defaults(end+1,:) = {'firstSpot', 'integer', 1};
defaults(end+1,:) = {'lastSpot', 'integer', []};
defaults(end+1,:) = {'numParallel','integer',1};
defaults(end+1,:) = {'nmXYpix','positive',154};
defaults(end+1,:) = {'nmZpix','positive',100};
defaults(end+1,:) = {'figHandles','freeType',[]}; 
defaults(end+1,:) = {'goodHybes','freeType',[]};
defaults(end+1,:) = {'eTable','freeType',[]};

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

defaults(end+1,:) = {'datZframes','positive',40};
pars = ParseVariableArguments(varargin, defaults, mfilename);


% turn off some annoying and non-relevant warnings
warning('off','MATLAB:prnRenderer:opengl');
warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');
warning('off','MATLAB:singularMatrix');

%% Fit 3D Gaussians in matlab using FitPsf3D
numSpots = length(fidSpots)- sum(cellfun(@isempty,fidSpots));

[numHybes,numDataBits] = size(dataSpots{pars.firstSpot});

if isempty(pars.lastSpot) || pars.lastSpot > numSpots
    pars.lastSpot = numSpots;
end


pars.goodHybes = logical(pars.goodHybes);
if isempty(pars.goodHybes)
    pars.goodHybes = true(numHybes,1);
end
    
numGoodHybes = sum(pars.goodHybes);
hybNum = find(pars.goodHybes);
if ~isempty(pars.eTable)
    dataType = pars.eTable.DataType;
    if ~iscell(pars.eTable.Bit(1))
        bitNum = pars.eTable.Bit;
    else
        bitNum = cellfun(@(x) str2double(strsplit(x,',')),pars.eTable.Bit,'UniformOutput',false); % added st2double
        bitNum = cat(1,bitNum{:});
        bitNum = bitNum(pars.goodHybes,:);
    end
else
    bitNum = 1:numHybes;
end



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

%%
% parfor loops don't like structures
%    "Valid indices for "pars" are restricted"
firstSpot = pars.firstSpot;
lastSpot = pars.lastSpot;
fov = pars.fov; 
lociXYx = pars.lociXY(:,1);
lociXYy = pars.lociXY(:,2);
showPlots = pars.showPlots;
saveFolder = pars.saveFolder; 
saveData = pars.saveData;

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
fidTroubleshoot = pars.fidTroubleshoot;

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
datTroubleshoot = pars.datTroubleshoot;

datZframes = pars.datZframes;

nmXYpix = pars.nmXYpix;
nmZpix = pars.nmZpix;
verbose = pars.verbose;
numParallel = pars.numParallel;

if pars.numParallel > 1
   fidTroubleshoot = false;
   datTroubleshoot = false;
   vis = 'off';
   if isempty(gcp('nocreate'))
        parpool(pars.numParallel); 
   end
   if pars.verbose
      disp('running in parallel mode, setting figure visibility to:');
      disp(vis);
   end
else
    vis = 'on';
end
%% Main loop
cellTable = cell(numSpots*numGoodHybes,1);
tic
if verbose
    disp(['fitting spots ',num2str(firstSpot),' to ',num2str(lastSpot) '.']);
    disp(['using ',num2str(numParallel),' processors']);
end




%%

% currently blows up gradually in memory. Need to fix. 
if numParallel < 2
    for s=firstSpot:lastSpot % another good place for a parfor loop
        fidTable = {}; 
        datTable = {};
        for h=1:numGoodHybes % h=3
            if fidTroubleshoot; figure(11); clf; end
            dTable = FitPsf3D(squeeze(fidSpots{s}{hybNum(h)}(:,:,:)),...
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
                'xyUnitConvert',nmXYpix,...
                'zUnitConvert',nmZpix,...
                'troubleshoot',fidTroubleshoot); 

            % keep track of hybe, spot number, and is fiduial. 
            dTable.panel = h*ones(length(dTable.x),1);
            dTable.hybe = hybNum(h)*ones(length(dTable.x),1);
            dTable.dataType = repmat(string(dataType{hybNum(h)}),length(dTable.x),1); 
            dTable.bit = repmat(bitNum(hybNum(h),1),length(dTable.x),1); 
            dTable.s = s*ones(length(dTable.x),1);
            dTable.isfid = true(length(dTable.x),1);
            dTable.fov = fov*ones(length(dTable.x),1);
            dTable.locusX = lociXYx(s)*ones(length(dTable.x),1);
            dTable.locusY = lociXYy(s)*ones(length(dTable.x),1);
            fidTable{h} = dTable;  %#ok<AGROW> % Save

            % Plot spots
            if fidTroubleshoot; pause; end
            if showPlots
                % xy fig
                tileName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_fid_tileXY'];
                tileSaveName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_fid_tileXY_Fit'];
                if h==1
                    uiopen([saveFolder,filesep,tileName,'.fig'],1); 
                    tileFidxy = gcf;
                    tileFidxy.Visible = vis;
                    tileFidxy.Name = tileSaveName;
                    panelsFidxy = flipud(tileFidxy.Children);  % (need to store these because calling axes() set the current axes as the first child, changing the order). 
                    if ~isempty(tileLocs{1})
                       tileFidxy.Position =  tileLocs{1};
                    end
                end
                axes(panelsFidxy(h));    %#ok<*LAXES>  
                plot(dTable.x/nmXYpix,dTable.y/nmXYpix,'yo','MarkerSize',15);

                % xz fig
                ztileName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_fid_tileXZ'];
                ztileSaveName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_fid_tileXZ_Fit'];
                if h==1 
                    uiopen([saveFolder,filesep,ztileName,'.fig'],1);  
                    tileFidxz = gcf;
                    tileFidxz.Visible = vis;
                    tileFidxz.Name = ztileSaveName;
                    panelsFidxz = flipud(tileFidxz.Children);
                    if ~isempty(tileLocs{2})
                       tileFidxz.Position =  tileLocs{2};
                    end
                end
                axes(panelsFidxz(h)); hold on; 
                plot(dTable.x/nmXYpix,dTable.z/nmZpix,'yo','MarkerSize',15);
               % text(dTable.x/pars.nmXYpix,dTable.z/pars.nmZpix,num2str(h),'FontSize',15,'color','w');
            end 
        end
        numZframes = size(squeeze(fidSpots{s}{hybNum(1)}(:,:,:)),3);
        meanZ = cellfun(@(x) x.z,fidTable,'UniformOutput',false); 
        meanZ = floor(nanmean(cat(1,meanZ{:}))/nmZpix);
        datZstart = meanZ - floor(datZframes/2);
        datZstop = meanZ + ceil(datZframes/2)-1;
        datZstart(datZstart<1) = 1;
        datZstop(datZstop>numZframes) = numZframes;
        datZstart(isnan(datZstart)) =1;
        datZstop(isnan(datZstop)) = numZframes;
        p = 0; % subplot index
        for n=1:numDataBits
            for h=1:numGoodHybes 
                p=p+1;
                if datTroubleshoot; figure(11); clf; end
                dTable = FitPsf3D(squeeze(dataSpots{s}{hybNum(h),n}(:,:,datZstart:datZstop)),...
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
                    'xyUnitConvert',nmXYpix,...
                    'zUnitConvert',nmZpix,...
                    'troubleshoot',datTroubleshoot);

                % keep track of hybe, spot number, and is fiduial. 
                dTable.z = dTable.z + datZstart*nmZpix;
                dTable.panel = p*ones(length(dTable.x),1);
                dTable.hybe = hybNum(h)*ones(length(dTable.x),1);
                dTable.dataType = repmat(string(dataType{hybNum(h)}),length(dTable.x),1);  % needs to be string not char
                dTable.bit = repmat(bitNum(hybNum(h),n),length(dTable.x),1);
                dTable.s = s*ones(length(dTable.x),1);
                dTable.isfid = false(length(dTable.x),1);  
                dTable.fov = fov*ones(length(dTable.x),1);
                dTable.locusX = lociXYx(s)*ones(length(dTable.x),1);
                dTable.locusY = lociXYy(s)*ones(length(dTable.x),1);
                datTable{p} = dTable; %#ok<AGROW>
                if showPlots
                     % xy fig
                    tileName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_dat_tileXY'];
                    tileSaveName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_dat_tileXY_Fit'];
                    if p==1
                        uiopen([saveFolder,filesep,tileName,'.fig'],1);  
                        tileDatxy = gcf;
                        tileDatxy.Visible = vis;
                        tileDatxy.Name = tileSaveName;
                        panelsDatxy = flipud(tileDatxy.Children);
                        if ~isempty(tileLocs{3})
                           tileDatxy.Position =  tileLocs{3};
                        end
                    end
                    axes(panelsDatxy(p)); hold on;                 
                    plot(dTable.x/nmXYpix,dTable.y/nmXYpix,'yo','MarkerSize',15);


                    % xz fig
                    ztileName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_dat_tileXZ'];
                    ztileSaveName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_dat_tileXZ_Fit'];
                    if p==1 
                        uiopen([saveFolder,filesep,ztileName,'.fig'],1);  
                        tileDatxz = gcf;
                        tileDatxz.Visible = vis;
                        tileDatxz.Name = ztileSaveName;
                        panelsDatxz = flipud(tileDatxz.Children);
                        if ~isempty(tileLocs{4})
                           tileDatxz.Position =  tileLocs{4};
                        end
                    end
                    axes(panelsDatxz(p)); hold on; 
                    plot(dTable.x/nmXYpix,dTable.z/nmZpix,'yo','MarkerSize',15);
                    plot([0,30],[datZstart,datZstart],'-','color',[.8 .8 .8]);
                    plot([0,30],[datZstop,datZstop],'-','color',[.8 .8 .8]);
                end       
            end
        end

        cellTable{s} = cat(1,datTable{:},fidTable{:}); 
        if saveData && showPlots
            warning('off','MATLAB:prnRenderer:opengl');
            SetFigureSavePath(saveFolder,'verbose',false); 
            SaveFigure(tileFidxy,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);
            SaveFigure(tileFidxz,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);
            SaveFigure(tileDatxy,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);
            SaveFigure(tileDatxz,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);
        end
        if saveData
            tableName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_spottable.csv'];
            writetable(cellTable{s},[saveFolder,tableName]); 
            if verbose
               disp(['wrote ',tableName]); 
            end
        end
        
        if numParallel < 2 && verbose
            disp(['Fitting 3D gaussian spots: ' num2str(100*s/numSpots,3),' % complete...']);
        end
    end

 
else 
    % =================================================================== % 
    %  Identical to above but in parfor                                   %
    % =================================================================== % 
    parfor s=firstSpot:lastSpot % parfor %  another good place for a parfor loop
        fidTable = {}; 
        datTable = {};
        for h=1:numGoodHybes % h=3
            if fidTroubleshoot; figure(11); clf; end
            dTable = FitPsf3D(squeeze(fidSpots{s}{hybNum(h)}(:,:,:)),...
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
                'xyUnitConvert',nmXYpix,...
                'zUnitConvert',nmZpix,...
                'troubleshoot',fidTroubleshoot); %#ok<PFBNS>

            % keep track of hybe, spot number, and is fiduial. 
            dTable.panel = h*ones(length(dTable.x),1);
            dTable.hybe = hybNum(h)*ones(length(dTable.x),1);
            dTable.dataType = repmat(string(dataType{hybNum(h)}),length(dTable.x),1); %#ok<PFBNS>
            dTable.bit = repmat(bitNum(hybNum(h),1),length(dTable.x),1); %#ok<PFBNS>
            dTable.s = s*ones(length(dTable.x),1);
            dTable.isfid = true(length(dTable.x),1);
            dTable.fov = fov*ones(length(dTable.x),1);
            dTable.locusX = lociXYx(s)*ones(length(dTable.x),1);
            dTable.locusY = lociXYy(s)*ones(length(dTable.x),1);
            fidTable{h} = dTable; 

            % Plot spots
            if fidTroubleshoot; pause; end
            if showPlots
                % xy fig
                tileName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_fid_tileXY'];
                tileSaveName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_fid_tileXY_Fit'];
                if exist([tileSaveName,'.fig'],'file')>0
                    alreadyWritten = true;
                else
                    alreadyWritten = false;
                end
                if ~alreadyWritten    
                    if h==1
                        uiopen([saveFolder,filesep,tileName,'.fig'],1); 
                        tileFidxy = gcf;
                        tileFidxy.Visible = vis;
                        tileFidxy.Name = tileSaveName;
                        panelsFidxy = flipud(tileFidxy.Children);  % (need to store these because calling axes() set the current axes as the first child, changing the order). 
                        if ~isempty(tileLocs{1}) %#ok<PFBNS>
                           tileFidxy.Position =  tileLocs{1};
                        end
                    end
                    axes(panelsFidxy(h));    %#ok<*LAXES>  
                    plot(dTable.x/nmXYpix,dTable.y/nmXYpix,'yo','MarkerSize',15);

                    % xz fig
                    ztileName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_fid_tileXZ'];
                    ztileSaveName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_fid_tileXZ_Fit'];
                    if h==1 
                        uiopen([saveFolder,filesep,ztileName,'.fig'],1);  
                        tileFidxz = gcf;
                        tileFidxz.Visible = vis;
                        tileFidxz.Name = ztileSaveName;
                        panelsFidxz = flipud(tileFidxz.Children);
                        if ~isempty(tileLocs{2})
                           tileFidxz.Position =  tileLocs{2};
                        end
                    end
                    axes(panelsFidxz(h)); hold on; 
                    plot(dTable.x/nmXYpix,dTable.z/nmZpix,'yo','MarkerSize',15);
                   % text(dTable.x/pars.nmXYpix,dTable.z/pars.nmZpix,num2str(h),'FontSize',15,'color','w');
                end
            end
        end
        numZframes = size(squeeze(fidSpots{s}{hybNum(1)}(:,:,:)),3);
        meanZ = cellfun(@(x) x.z,fidTable,'UniformOutput',false); 
        meanZ = floor(nanmean(cat(1,meanZ{:}))/nmZpix);
        datZstart = meanZ - floor(datZframes/2);
        datZstop = meanZ + ceil(datZframes/2)-1;
        datZstart(datZstart<1) = 1;
        datZstop(datZstop>numZframes) = numZframes;
        datZstart(isnan(datZstart)) =1;
        datZstop(isnan(datZstop)) = numZframes;
        p = 0; % subplot index
        for n=1:numDataBits
            for h=1:numGoodHybes 
                p=p+1;
                if datTroubleshoot; figure(11); clf; end
                dTable = FitPsf3D(squeeze(dataSpots{s}{hybNum(h),n}(:,:,datZstart:datZstop)),...
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
                    'xyUnitConvert',nmXYpix,...
                    'zUnitConvert',nmZpix,...
                    'troubleshoot',datTroubleshoot);

                % keep track of hybe, spot number, and is fiduial. 
                dTable.z = dTable.z + datZstart*nmZpix;
                dTable.panel = p*ones(length(dTable.x),1);
                dTable.hybe = hybNum(h)*ones(length(dTable.x),1);
                dTable.dataType = repmat(string(dataType{hybNum(h)}),length(dTable.x),1);  % needs to be string not char
                dTable.bit = repmat(bitNum(hybNum(h),n),length(dTable.x),1);
                dTable.s = s*ones(length(dTable.x),1);
                dTable.isfid = false(length(dTable.x),1);  
                dTable.fov = fov*ones(length(dTable.x),1);
                dTable.locusX = lociXYx(s)*ones(length(dTable.x),1);
                dTable.locusY = lociXYy(s)*ones(length(dTable.x),1);
                datTable{p} = dTable; 
                
                if showPlots
                    tileName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_dat_tileXY'];
                    tileSaveName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_dat_tileXY_Fit'];
                    if exist([tileSaveName,'.fig'],'file') > 0
                        alreadyWritten = true;
                    else
                        alreadyWritten = false;
                    end
                    if ~alreadyWritten
                         % xy fig
                        if p==1
                            uiopen([saveFolder,filesep,tileName,'.fig'],1);  
                            tileDatxy = gcf;
                            tileDatxy.Visible = vis;
                            tileDatxy.Name = tileSaveName;
                            panelsDatxy = flipud(tileDatxy.Children);
                            if ~isempty(tileLocs{3})
                               tileDatxy.Position =  tileLocs{3};
                            end
                        end
                        axes(panelsDatxy(p)); hold on;                 
                        plot(dTable.x/nmXYpix,dTable.y/nmXYpix,'yo','MarkerSize',15);


                        % xz fig
                        ztileName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_dat_tileXZ'];
                        ztileSaveName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_dat_tileXZ_Fit'];
                        if p==1 
                            uiopen([saveFolder,filesep,ztileName,'.fig'],1);  
                            tileDatxz = gcf;
                            tileDatxz.Visible = vis;
                            tileDatxz.Name = ztileSaveName;
                            panelsDatxz = flipud(tileDatxz.Children);
                            if ~isempty(tileLocs{4})
                               tileDatxz.Position =  tileLocs{4};
                            end
                        end
                        axes(panelsDatxz(p)); hold on; 
                        plot(dTable.x/nmXYpix,dTable.z/nmZpix,'yo','MarkerSize',15);
                        plot([0,30],[datZstart,datZstart],'-','color',[.8 .8 .8]);
                        plot([0,30],[datZstop,datZstop],'-','color',[.8 .8 .8]);
                    end
                end
            end
        end

        cellTable{s} = cat(1,datTable{:},fidTable{:}); 
        if saveData && showPlots
            warning('off','MATLAB:prnRenderer:opengl');
            SetFigureSavePath(saveFolder,'verbose',false); 
            SaveFigure(tileFidxy,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);
            SaveFigure(tileFidxz,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);
            SaveFigure(tileDatxy,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);
            SaveFigure(tileDatxz,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);
            tileFidxy = []; % attempt to clear memory, prevent RAM fillup.
            tileFidxz = [];
            tileDatxy = [];
            tileDatxz = [];
            dTable = []; %#ok<*NASGU>
            datTable = [];
            fidTable = [];
        end
        if saveData
            tableName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_L',num2str(lociXYx(s)),'-',num2str(lociXYy(s)),'_spottable.csv'];
            writetable(cellTable{s},[saveFolder,tableName]); 
            if verbose
               disp(['wrote ',tableName]); 
            end
        end

        if numParallel < 2 && verbose
            disp(['Fitting 3D gaussian spots: ' num2str(100*s/numSpots,3),' % complete...']);
        end
        
    end
end







t= toc;
spotTable = cat(1,cellTable{:});
if pars.verbose
   disp(['FitPsf3D fit ',num2str(height(spotTable)),' spots in ',num2str(t/60),' minutes']);
end
 



%% A little more elegant than copy paste
% 
% if numParallel < 2
%     for s=firstSpot:lastSpot
%         cellTable{s} = FitSpots();
%     end
% else
%     parfor s=firstSpot:lastSpot
%         cellTable{s} = FitSpots();
%     end
% end
% 
% function cellTable = FitSpots(fov,lociXYx,lociXYy,showPlots,saveFolder,saveData,...
%     fidMinPeakHeight,fidPeakBlur,fidCameraBackground,fidMaxFitWidth,...
%     fidMinSep,fidKeepBrightest,fidRelativeHeight,fidMinHBratio)
% 
