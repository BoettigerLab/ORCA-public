function [dat_table,distMap,allShifts,allShifts_h,distMap2D] = UniqueSpotTableToMap(spotTable,varargin)
%%  single point matcher
% Matching assumptions
%  Max of 1 localization per hybe in ROI (by prior filtering)
%  All localizations in ROI are on target

global scratchPath 

% Parse Variable Input
defaults = cell(0,3);
defaults(end+1,:) = {'fov','integer',1};
defaults(end+1,:) = {'firstSpot','integer',1};
defaults(end+1,:) = {'lastSpot','integer',[]};
defaults(end+1,:) = {'saveData','boolean',false};
defaults(end+1,:) = {'showplots','boolean',true};
defaults(end+1,:) = {'stopOnError','boolean',false};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'numParallel','integer',1};
defaults(end+1,:) = {'saveFolder', 'string', scratchPath};
defaults(end+1,:) = {'renderInterpMethod',{'pchip','spline','linear','nearest','next','previous','cubic'},'pchip'};
defaults(end+1,:) = {'renderTubeRadius','positive',25};
defaults(end+1,:) = {'renderBallRadius','positive',30};
defaults(end+1,:) = {'minFracFidPoints','fraction',.5};
defaults(end+1,:) = {'alignToPrevious','boolean',false};
defaults(end+1,:) = {'refHybe','integer',1};
defaults(end+1,:) = {'fixUnique','boolean',false};
pars = ParseVariableArguments(varargin,defaults,mfilename);

SetFigureSavePath(pars.saveFolder);

fov = pars.fov;
numHybes = max(spotTable.hybe);
numSpots = max(spotTable.s); 
numBits = max(spotTable.bit);

if isempty(pars.lastSpot)
    pars.lastSpot = max([numSpots,pars.firstSpot]);
end

pars.numHybes = numHybes;
pars.numSpots = numSpots;
pars.numBits = numBits;

allShifts = nan(numHybes,3,numSpots); % cumulative shifts
allShifts_h = nan(numHybes,3,numSpots); % per hybe shift
distMap = nan(numHybes,numHybes,numSpots);
distMap2D = nan(numHybes,numHybes,numSpots);
dat_tables = cell(numSpots,1);
table_fidC = cell(numSpots,1);
table_datC = cell(numSpots,1);

if pars.numParallel > 1 && isempty(gcp('nocreate'))
    parpool(pars.numParallel); 
end

% pre-process the data to create cell arrays for parallel processing
for s = pars.firstSpot:pars.lastSpot
    try 
        %% Extracting data from table
        fid_s = spotTable.s == s & spotTable.isfid;
        table_fid = spotTable(fid_s,:); 
        dat_s = spotTable.s == s & ~spotTable.isfid;
        table_dat = spotTable(dat_s,:); 
        try
            table_fid = table_fid(table_fid.locusX==table_fid.locusX(1),:);
            if pars.fixUnique
                table_fid = unique(table_fid,'rows');
            end
        catch
            table_fid = table();
        end
        try
            table_dat = table_dat(table_dat.locusX==table_dat.locusX(1),:);
            if pars.fixUnique
                table_dat = unique(table_dat,'rows');
            end
        catch
            table_dat = table();
        end
    catch
    end
    table_fidC{s} = table_fid;
    table_datC{s} = table_dat;
end

if pars.numParallel < 2
    vis = 'on';
    for s=pars.firstSpot:pars.lastSpot    
        table_fid = table_fidC{s};
        table_dat = table_datC{s};
        if height(table_fid) < pars.minFracFidPoints*numHybes || height(table_dat) < 5
           continue 
        end
        [dat_table_s,distMap_s,distMap2D_s,shifts,shifts_h] = MainLoop(table_fid,table_dat,pars,s,vis);
        distMap(:,:,s) = distMap_s; % squareform(pdist(dat_xyz_aligned));
        distMap2D(:,:,s) = distMap2D_s;% squareform(pdist(dat_xyz_aligned(:,1:2)));
        dat_tables{s} = dat_table_s;
        allShifts(:,:,s) = shifts;
        allShifts_h(:,:,s) = shifts_h;
    end
else
    vis = 'off';
    warning('off','MATLAB:prnRenderer:opengl');   
    warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');
    parfor s=pars.firstSpot:pars.lastSpot    
        table_fid = table_fidC{s};
        table_dat = table_datC{s};
        if height(table_fid) < pars.minFracFidPoints*numHybes || height(table_dat) < 5
           continue 
        end
        try
            [dat_table_s,distMap_s,distMap2D_s,shifts,shifts_h] = MainLoop(table_fid,table_dat,pars,s,vis);
            distMap(:,:,s) = distMap_s; % squareform(pdist(dat_xyz_aligned));
            distMap2D(:,:,s) = distMap2D_s;% squareform(pdist(dat_xyz_aligned(:,1:2)));
            dat_tables{s} = dat_table_s;
            allShifts(:,:,s) = shifts;
            allShifts_h(:,:,s) = shifts_h;
        catch
            warning(['error on fov ',num2str(pars.fov),' spot ',num2str(s)]);
        end
    end
end

dat_table = cat(1,dat_tables{:});

function  [dat_table_s,distMap_s,distMap2D_s,shifts,shifts_h] = MainLoop(table_fid,table_dat,pars,s,vis)
    try
        SetFigureSavePath(pars.saveFolder,'verbose',false);
        numHybes = pars.numHybes;
        fov = pars.fov;
        numBits = pars.numBits;
        
        fid_xyz = nan(numHybes,3);
        fid_xyz(table_fid.hybe,:) = table_fid{:,1:3}; % sort in hybe order

        dat_xyz = nan(numHybes,3);
        dat_xyz(table_dat.hybe,:) = table_dat{:,1:3};
        try
            dat_bitName = cell(numHybes,1);
            dat_bitName(table_dat.hybe) = cellstr(num2str(table_dat.bit));
            dat_bit = nan(numBits,1);
            dat_bit(table_dat.bit,:) = table_dat.bit(:);
        catch
            dat_bitName = cellstr(num2str( (1:numHybes)'  ));
        end

        % allDat(:,:,s) = dat_xyz;
        %% compute shifts
        shifts = zeros(numHybes,3);
        shifts_h = zeros(numHybes,3);
        for h=2:numHybes %   h =2;        
            if pars.alignToPrevious
                % deal with missing data
                alignHybe = h-1;
                good = ~isnan(fid_xyz(alignHybe,1));
                while ~good && alignHybe > 1
                    alignHybe = alignHybe -1;
                    good = ~isnan(fid_xyz(alignHybe,1));
                end
                if ~isnan(fid_xyz(h,1))
                    shift_h = fid_xyz(h,:) - fid_xyz(alignHybe,:);
                else
                    shift_h = [0,0,0];
                end
                shifts(h,:) = shifts(h-1,:) + shift_h;
                shifts_h(h,:) = shift_h;
            else
                alignHybe = pars.refHybe;
                shifts(h,:) = fid_xyz(h,:) - fid_xyz(alignHybe,:);
                shifts_h(h,:) = shifts(h,:) - shifts(h-1,:);
            end
        end

        % allShifts(:,:,s) = shifts;
        % allShifts_h(:,:,s) = shifts_h;
        %% apply shifts
        fid_xyz_aligned = fid_xyz;
        dat_xyz_aligned = dat_xyz;
        for h=1:numHybes
            fid_xyz_aligned(h,:) = fid_xyz_aligned(h,:) - shifts(h,:);
            dat_xyz_aligned(h,:) = dat_xyz_aligned(h,:) - shifts(h,:);
        end

        distMap_s = squareform(pdist(dat_xyz_aligned));
        distMap2D_s = squareform(pdist(dat_xyz_aligned(:,1:2)));

        xc = dat_xyz_aligned(table_dat.hybe,1);
        yc = dat_xyz_aligned(table_dat.hybe,2);
        zc = dat_xyz_aligned(table_dat.hybe,3);

        [~,fidID] = intersect(table_fid.hybe,table_dat.hybe);
        [~,ia,ib] = setxor(table_fid.hybe,table_dat.hybe);
        fid_table =table_fid(fidID,1:16);

        if ~isempty(ib) || ~isempty(ia) % if there are data without matching fiducial
            % % Option 1: remove these data points 
            % table_dat(ib,:) = []; 
            % xc(ib) = []; yc(ib) = []; zc(ib) = [];

            % Option 2: insert NaNs to matching fid table
            fid_table = table();
            for t=1:height(table_dat)
                m = table_fid.hybe == table_dat.hybe(t);
                if sum(m) == 0
                    blankTable = cell2table(repmat({NaN},1,16));
                    blankTable.Properties.VariableNames = table_fid.Properties.VariableNames(1:16);
                    fid_table = cat(1,fid_table,blankTable);
                end
                fid_table = cat(1,fid_table,table_fid(m,1:16));
            end
        end

        varNames = fid_table.Properties.VariableNames;
        varNames = cellfun(@(x) ['fid_',x],varNames,'UniformOutput',false);
        fid_table.Properties.VariableNames = varNames;


        try
            dat_table_s = cat(2,table_dat,table(xc,yc,zc),fid_table);
        catch
            error('failed cat');
        end

        %% just plotting 
        if pars.showplots
            scatterFig = figure(10); clf; 
            scatterFig.Visible = vis;
            subplot(1,2,1);
            scatter(dat_xyz(:,1),dat_xyz(:,2),2,jet(numHybes)); hold on;
            scatter(dat_xyz_aligned(:,1),dat_xyz_aligned(:,2),25,jet(numHybes));
            scatter(fid_xyz(:,1),fid_xyz(:,2),5,jet(numHybes),'d');
            scatter(fid_xyz_aligned(:,1),fid_xyz_aligned(:,2),25,jet(numHybes),'d');
            subplot(1,2,2);
            scatter(dat_xyz(:,3),dat_xyz(:,2),2,jet(numHybes)); hold on;
            scatter(dat_xyz_aligned(:,3),dat_xyz_aligned(:,2),25,jet(numHybes));
            scatter(fid_xyz(:,3),fid_xyz(:,2),5,jet(numHybes),'d');
            scatter(fid_xyz_aligned(:,3),fid_xyz_aligned(:,2),25,jet(numHybes),'d');
            imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_scatterFig'];
            try
                SaveFigure(scatterFig,'name',imName,'overwrite',true,'formats',{'png','fig'},'saveData',pars.saveData);
            catch er
               warning(er.getReport)
            end

            tubeFig = figure(11); clf;
            tubeFig.Visible = vis;
            noDat = isnan(dat_xyz_aligned(:,1));
            plotDat = dat_xyz_aligned(~noDat,:);
            datHybes = find(~noDat);
            PlotTube(plotDat,'r',pars.renderTubeRadius,'interpPts',10,'method',pars.renderInterpMethod,'colormap',jet(numHybes)); 
            set(gca,'color','k'); hold on;
            alpha .5;
            freezeColors; 
            PlotSpheres(plotDat,'r',pars.renderBallRadius,'color',jet(numHybes),'alpha',.55);
            text(plotDat(:,1)+pars.renderBallRadius,plotDat(:,2),plotDat(:,3),dat_bitName(datHybes),'color','w'); % cellstr(num2str(datHybes))
            axis equal;
            imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_tubeFig'];
            try 
                SaveFigure(tubeFig,'name',imName,'overwrite',true,'formats',{'png','fig'},'saveData',pars.saveData);
            catch er
                warning(er.getReport)
            end

            mapFig = figure(12); clf;
            mapFig.Visible = vis;
            imagesc(distMap_s); colorbar; caxis([0,700]);
            GetColorMap('redToWhiteSat');
            imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_mapFig'];
            try
                SaveFigure(mapFig,'name',imName,'overwrite',true,'formats',{'png','fig'},'saveData',pars.saveData);       
            catch er
                warning(er.getReport)
            end
        end
    catch er
        disp(er.message);
        warning(['error on spot ',num2str(s)]);
    end


