function [cy5LinkTables,shifts,pars] = ChrTracer_SpotTableToMap(spotTable,varargin)
% [cy5Chains,cy3Chains,M] = ChrTracer_SpotTableToMap(spotTable)
% 

global scratchPath;

% Parse Variable Input
defaults = cell(0,3);
defaults(end+1,:) = {'fov','integer',1};
defaults(end+1,:) = {'eTable','array',[]};
defaults(end+1,:) = {'firstSpot','integer',1};
defaults(end+1,:) = {'lastSpot','integer',[]};
defaults(end+1,:) = {'saveData','boolean',false};
defaults(end+1,:) = {'stopOnError','boolean',false};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'veryverbose','boolean',false};
defaults(end+1,:) = {'numParallel','integer',1};
defaults(end+1,:) = {'saveFolder', 'string', scratchPath};
defaults(end+1,:) = {'maxJump','positive',1000}; % in nm 
defaults(end+1,:) = {'maxDistFromStart','positive',2000}; % in nm 
defaults(end+1,:) = {'renderInterpMethod',{'pchip','spline','linear','nearest','next','previous','cubic'},'pchip'};
defaults(end+1,:) = {'renderTubeRadius','positive',25};
defaults(end+1,:) = {'renderBallRadius','positive',30};
defaults(end+1,:) = {'minFracFidPoints','fraction',.5};
pars = ParseVariableArguments(varargin,defaults,mfilename);

%% Main Function
numSpots = max(spotTable.s);
numHybes = max(spotTable.hybe);
fov = pars.fov;
maxJump = pars.maxJump;
maxDistFromStart = pars.maxDistFromStart;  
renderTubeRadius = pars.renderTubeRadius;
renderBallRadius = pars.renderBallRadius;
renderInterpMethod = pars.renderInterpMethod;
stopOnError = pars.stopOnError;
saveData = pars.saveData;
verbose = pars.verbose;
saveFolder = pars.saveFolder;
minFracFidPoints = pars.minFracFidPoints;

if isempty(pars.lastSpot)
    pars.lastSpot = max([numSpots,pars.firstSpot]);
end
firstSpot = pars.firstSpot;
lastSpot = pars.lastSpot;


cy5LinkTables = cell(numSpots,1);
cy3Chains = cell(numSpots,1);
cy5Chains = cell(numSpots,1); 
shifts = cell(numSpots,1); 

cy3Tables = cell(numSpots,1);
cy5Tables = cell(numSpots,1);
for s=pars.firstSpot:pars.lastSpot
    cy3_s = spotTable.s == s & spotTable.isfid;
    cy3Tables{s}  = spotTable(cy3_s,:); 
    cy5_s = spotTable.s == s & ~spotTable.isfid;
    cy5Tables{s}  = spotTable(cy5_s,:); 
end


if pars.numParallel < 2
    for s=firstSpot:lastSpot 
            shifts{s} = nan(numHybes,3); 
            table_cy3 = cy3Tables{s};
            xyz_cy3 = table_cy3{:,1:3};
            hybe_cy3 = table_cy3.hybe ;

            if isempty(table_cy3)
                if verbose
                    cprintf([1 .5 0], ['warning: no data in fiducial table for spot ',num2str(s),'. Skipping spot.']);
                end
                continue;
            end
            
            if height(table_cy3) < minFracFidPoints*numHybes
                if verbose
                    cprintf([1 .5 0], ['warning: less than ',num2str(minFracFidPoints*numHybes),' fiducial points for spot ',num2str(s),'. Skipping spot.']); 
                end
                continue
            end
            
            % Link molecules
            [chainPts_cy3,chainQuality_cy3,chainIdx_cy3] = LinkMolecules(...
                xyz_cy3,hybe_cy3,...
                'maxJump',maxJump,...
                'maxDistFromStart',maxDistFromStart,...
                'brightness',table_cy3.h);

            % compute chain lengths
            cy3ChainLengths = cellfun(@(x) sum(~isnan(sum(x,2))), chainPts_cy3);

            % discard super-short chains, be sure to keep at least 1 chain
            chainPts_cy3 = chainPts_cy3(cy3ChainLengths > min([10,max(cy3ChainLengths)-1]) );
            chainIdx_cy3 = chainIdx_cy3(cy3ChainLengths > min([10,max(cy3ChainLengths)-1]) );

            % compute chain Centroids (average positions)
            cy3ChainPos = cell2mat(cellfun(@nanmean,chainPts_cy3,'UniformOutput',false));


            %% Link cy5 chains
            % For Future: Maybe this should be done after cy3 based drift correction is complete? 


            % extract the data from our spot Tables
            cy5_s = spotTable.s == s & ~spotTable.isfid;
            table_cy5 = spotTable(cy5_s,:);
            xyz_cy5 = table_cy5{:,1:3};
            hybe_cy5 = table_cy5.hybe;

            
            if isempty(table_cy5)
                if verbose
                    cprintf([1 .5 0], ['warning: no data in cy5 table for spot ',num2str(s)]);
                end
                continue;
            end
            
            % Link molecules cy5
            [chainPts_cy5,chainQuality_cy5, chainIdx_cy5] = LinkMolecules(...
                xyz_cy5,hybe_cy5,...
                'maxJump',maxJump,...
                'maxDistFromStart',maxDistFromStart,...
                'startPos',cy3ChainPos,....
                'brightness',table_cy5.h);

            % compute chain lengths cy5
            cy5ChainLengths = cellfun(@(x) sum(~isnan(sum(x,2))), chainPts_cy5);

            % discard short chains
            chainPts_cy5 = chainPts_cy5(cy5ChainLengths > min([3,max(cy5ChainLengths)-1]) ); 
            chainIdx_cy5 = chainIdx_cy5(cy5ChainLengths > min([3,max(cy5ChainLengths)-1]) );

            % compute centroids
            cy5ChainPos = cell2mat(cellfun(@nanmean,chainPts_cy5,'UniformOutput',false));


            %% Match 
            % sort so we have 2 cell arrays, equal length, matched order
            % two cy5 chains may use the same cy3 chain

            [ids,ds] = knnsearch(cy3ChainPos(:,:),cy5ChainPos(:,:));
            % eliminate cy5 without cy3 pair
            ids(ds>maxDistFromStart) = [];  % might want this to be a unique parameter 
            chainPts_cy5(ds>maxDistFromStart) = [];
            chainPts_cy3 = chainPts_cy3(ids);
            chainIdx_cy3 = chainIdx_cy3(ids); 

            numChains = length(chainPts_cy3);
            numCy3Pts = max(hybe_cy3); % length(chainPts_cy3{1}); 
            numCy5Pts = max(hybe_cy5); % length(chainPts_cy5{1});
            %% Update tables
            % now these tables match (which they may have already done). 
            linkTable_cy3 = cell(1,numChains);
            linkTable_cy5 = cell(1,numChains); 
            for n=1:numChains
                linkTable_cy3{n} = table_cy3(nonzeros(chainIdx_cy3{n}),:);
                chainID = n*ones(height(linkTable_cy3{n}),1);
                linkTable_cy3{n} = cat(2,linkTable_cy3{n},table(chainID));
                linkTable_cy5{n} = table_cy5(nonzeros(chainIdx_cy5{n}),:);
                chainID = n*ones(height(linkTable_cy5{n}),1);
                linkTable_cy5{n} = cat(2,linkTable_cy5{n},table(chainID));
            end
            figure(1); clf; 
            subplot(1,2,1); RenderSpotTable(linkTable_cy3{1},'gain',2);
            subplot(1,2,2); RenderSpotTable(linkTable_cy5{1},'gain',3); 




            %% correct via cy3
            chainPts_cy5_align = chainPts_cy5;
            fidChain = 1;
            cntrlChain0 = chainPts_cy3{1};
            cntrlChain = chainPts_cy3{end};

            % compute shifts
            numDataHybes = min([numHybes,numCy3Pts]);
            xshift = 0;
            yshift = 0;
            zshift = 0;
            xshifts = zeros(numDataHybes,1);
            yshifts = zeros(numDataHybes,1);
            zshifts = zeros(numDataHybes,1); 
            for h=2:numDataHybes %   h =2;
                % deal with missing data
                alignHybe = h-1;
                good = ~isnan(chainPts_cy3{fidChain}(alignHybe,1));
                while ~good && alignHybe > 1
                    alignHybe = alignHybe -1;
                    good = ~isnan(chainPts_cy3{fidChain}(alignHybe,1));
                end
                if ~isnan(chainPts_cy3{fidChain}(h,1))
                    xshift_h = chainPts_cy3{fidChain}(h,1) - chainPts_cy3{fidChain}(alignHybe,1);
                    yshift_h = chainPts_cy3{fidChain}(h,2) - chainPts_cy3{fidChain}(alignHybe,2);
                    zshift_h = chainPts_cy3{fidChain}(h,3) - chainPts_cy3{fidChain}(alignHybe,3);
                else
                    xshift_h = 0;
                    yshift_h = 0;
                    zshift_h = 0;
                end
                xshift = xshift + xshift_h;
                yshift = yshift + yshift_h;
                zshift = zshift + zshift_h;
                xshifts(h) = xshift;
                yshifts(h) = yshift;
                zshifts(h) = zshift;
            end

            % apply shifts
            shifts{s}(1:numDataHybes,:) = [xshifts,yshifts,zshifts];
            for h=1:numDataHybes %
                cntrlChain0(h,1) = cntrlChain0(h,1) - xshifts(h);
                cntrlChain0(h,2) = cntrlChain0(h,2) - yshifts(h);
                cntrlChain0(h,3) = cntrlChain0(h,3) - zshifts(h);

                cntrlChain(h,1) = cntrlChain(h,1) - xshifts(h);
                cntrlChain(h,2) = cntrlChain(h,2) - yshifts(h);
                cntrlChain(h,3) = cntrlChain(h,3) - zshifts(h);

                hIds = linkTable_cy3{n}.hybe == h;
                xIds = StringFind(linkTable_cy3{n}.Properties.VariableNames,{'x','xL','xU'},'exactly',true);
                yIds = StringFind(linkTable_cy3{n}.Properties.VariableNames,{'y','yL','yU'},'exactly',true);
                zIds = StringFind(linkTable_cy3{n}.Properties.VariableNames,{'z','zL','zU'},'exactly',true);
                linkTable_cy3{n}{hIds,xIds} = linkTable_cy3{n}{hIds,xIds} - xshifts(h);
                linkTable_cy3{n}{hIds,yIds} = linkTable_cy3{n}{hIds,yIds} - yshifts(h);
                linkTable_cy3{n}{hIds,zIds} = linkTable_cy3{n}{hIds,zIds} - zshifts(h); 
            end

            % apply to all cy5 chains

            for h=1:min([numHybes,numCy5Pts,numCy3Pts])
                for n=1:length(chainPts_cy5_align)
                    chainPts_cy5_align{n}(h,1) = chainPts_cy5_align{n}(h,1) - xshifts(h);
                    chainPts_cy5_align{n}(h,2) = chainPts_cy5_align{n}(h,2) - yshifts(h);
                    chainPts_cy5_align{n}(h,3) = chainPts_cy5_align{n}(h,3) - zshifts(h);

                    hIds = linkTable_cy5{n}.hybe == h;
                    xIds = StringFind(linkTable_cy5{n}.Properties.VariableNames,{'x','xL','xU'},'exactly',true);
                    yIds = StringFind(linkTable_cy5{n}.Properties.VariableNames,{'y','yL','yU'},'exactly',true);
                    zIds = StringFind(linkTable_cy5{n}.Properties.VariableNames,{'z','zL','zU'},'exactly',true);
                    linkTable_cy5{n}{hIds,xIds} = linkTable_cy5{n}{hIds,xIds} - xshifts(h);
                    linkTable_cy5{n}{hIds,yIds} = linkTable_cy5{n}{hIds,yIds} - yshifts(h);
                    linkTable_cy5{n}{hIds,zIds} = linkTable_cy5{n}{hIds,zIds} - zshifts(h);    
                end
            end

            try
            spotFig = figure(1);  
            subplot(1,2,1); title('fiducial before correction');
            subplot(1,2,2); cla; RenderSpotTable(linkTable_cy5{1},'gain',3,'err',0);
            title('corrected data');
            imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_spotFig'];
            SaveFigure(spotFig,'name',imName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);        
            catch
                disp('error on ');
                disp(['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_spotFig']);
            end

            tubeFig = figure(100); clf;
            PlotTube(linkTable_cy5{1}{:,1:3},'r',renderTubeRadius,'interpPts',10,'method',renderInterpMethod,'colormap',jet(numCy5Pts)); 
            set(gca,'color','k'); hold on;
            alpha .5;
            freezeColors; 
            PlotSpheres(chainPts_cy5_align{1},'r',renderBallRadius,'color',jet(numCy5Pts),'alpha',.55);
            text(linkTable_cy5{1}{:,1}+renderBallRadius,linkTable_cy5{1}{:,2},linkTable_cy5{1}{:,3},num2str(linkTable_cy5{1}.hybe),'color','w');
            axis equal;
            imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_tubeFig'];
            SaveFigure(tubeFig,'name',imName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);

            %% create maps and save data
            numCy5Chains = 1;%  length(chainPts_cy5_align);
            hic5Fig = figure(101); clf; 
            for n = 1:numCy5Chains
                im = squareform(pdist(chainPts_cy5_align{n}(:,:)));
                cmap = GetColorMap('redToWhite',100); colormap(cmap);
                subplot(1,numCy5Chains,n); imagesc(im); colorbar; % caxis([0,3]);
            end
            imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_hic5Fig'];
            title(imName,'interpreter','none');
            SaveFigure(hic5Fig,'name',imName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);

            cy5Chains{s} = chainPts_cy5_align;
            cy3Chains{s} = cntrlChain;
            cy5LinkTables{s} = linkTable_cy5;
    end
    
    
else
%% ====================== Parallel Processing  =============================  

    if isempty(gcp('nocreate'))
        parpool(pars.numParallel); 
    end    

    parfor s=firstSpot:lastSpot  % parfor
        warning('off','MATLAB:prnRenderer:opengl');
        warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');
        SetFigureSavePath(saveFolder,'verbose',false);
        try
            %% link 
            % extract the data from our spot Tables
            shifts{s} = nan(numHybes,3);
            table_cy3 = cy3Tables{s};
            xyz_cy3 = table_cy3{:,1:3};
            hybe_cy3 = table_cy3.hybe ;

            if isempty(table_cy3)
                if verbose
                    cprintf([1 .5 0], ['warning: no data in cy3 table for spot ',num2str(s),'. Skipping spot.']);
                end
                continue;
            end
                        
            if height(table_cy3) < minFracFidPoints*numHybes
                if verbose
                    cprintf([1 .5 0], ['warning: less than ',num2str(minFracFidPoints*numHybes),' fiducial points for spot ',num2str(s),'. Skipping spot.']); 
                end
                continue
            end
            
            % Link molecules
            [chainPts_cy3,chainQuality_cy3,chainIdx_cy3] = LinkMolecules(...
                xyz_cy3,hybe_cy3,...
                'maxJump',maxJump,...
                'maxDistFromStart',maxDistFromStart,...
                'brightness',table_cy3.h);
           
            % compute chain lengths
            cy3ChainLengths = cellfun(@(x) sum(~isnan(sum(x,2))), chainPts_cy3);

            % discard super-short chains, be sure to keep at least 1 chain
            chainPts_cy3 = chainPts_cy3(cy3ChainLengths > min([10,max(cy3ChainLengths)-1]) );
            chainIdx_cy3 = chainIdx_cy3(cy3ChainLengths > min([10,max(cy3ChainLengths)-1]) );

            % compute chain Centroids (average positions)
            cy3ChainPos = cell2mat(cellfun(@nanmean,chainPts_cy3,'UniformOutput',false));


            %% Link cy5 chains
            % For Future: Maybe this should be done after cy3 based drift correction is complete? 


            % extract the data from our spot Tables
            table_cy5 = cy5Tables{s};
            xyz_cy5 = table_cy5{:,1:3};
            hybe_cy5 = table_cy5.hybe;

            if isempty(table_cy5)
                if verbose
                    cprintf([1 .5 0], ['warning: no data in cy5 table for spot ',num2str(s)]);
                end
                continue;
            end
            
            % Link molecules cy5
            [chainPts_cy5,chainQuality_cy5, chainIdx_cy5] = LinkMolecules(...
                xyz_cy5,hybe_cy5,...
                'maxJump',maxJump,...
                'maxDistFromStart',maxDistFromStart,...
                'startPos',cy3ChainPos,....
                'brightness',table_cy5.h);
            
            % compute chain lengths cy5
            cy5ChainLengths = cellfun(@(x) sum(~isnan(sum(x,2))), chainPts_cy5);

            % discard short chains
            chainPts_cy5 = chainPts_cy5(cy5ChainLengths > min([3,max(cy5ChainLengths)-1]) ); 
            chainIdx_cy5 = chainIdx_cy5(cy5ChainLengths > min([3,max(cy5ChainLengths)-1]) );

            % compute centroids
            cy5ChainPos = cell2mat(cellfun(@nanmean,chainPts_cy5,'UniformOutput',false));


            %% Match 
            % sort so we have 2 cell arrays, equal length, matched order
            % two cy5 chains may use the same cy3 chain

            [ids,ds] = knnsearch(cy3ChainPos(:,:),cy5ChainPos(:,:));
            % eliminate cy5 without cy3 pair
            ids(ds>maxDistFromStart) = [];  % might want this to be a unique parameter 
            chainPts_cy5(ds>maxDistFromStart) = [];
            chainPts_cy3 = chainPts_cy3(ids);
            chainIdx_cy3 = chainIdx_cy3(ids); 

            numChains = length(chainPts_cy3);
            numCy3Pts = max(hybe_cy3); % length(chainPts_cy3{1}); 
            numCy5Pts = max(hybe_cy5); % length(chainPts_cy5{1});
            
            %% Update tables
            % now these tables match (which they may have already done). 
            linkTable_cy3 = cell(1,numChains);
            linkTable_cy5 = cell(1,numChains); 
            for n=1:numChains
                linkTable_cy3{n} = table_cy3(nonzeros(chainIdx_cy3{n}),:);
                chainID = n*ones(height(linkTable_cy3{n}),1);
                linkTable_cy3{n} = cat(2,linkTable_cy3{n},table(chainID));
                linkTable_cy5{n} = table_cy5(nonzeros(chainIdx_cy5{n}),:);
                chainID = n*ones(height(linkTable_cy5{n}),1);
                linkTable_cy5{n} = cat(2,linkTable_cy5{n},table(chainID));
            end
%             figure(1); clf; 
%             subplot(1,2,1); RenderSpotTable(linkTable_cy3{1},'gain',2);
%             subplot(1,2,2); RenderSpotTable(linkTable_cy5{1},'gain',3); 




            %% correct via cy3
            chainPts_cy5_align = chainPts_cy5;
            fidChain = 1;
            cntrlChain0 = chainPts_cy3{1};
            cntrlChain = chainPts_cy3{end};

            % compute shifts
            numDataHybes = min([numHybes,numCy3Pts]);
            xshift = 0;
            yshift = 0;
            zshift = 0;
            xshifts = zeros(numDataHybes,1);
            yshifts = zeros(numDataHybes,1);
            zshifts = zeros(numDataHybes,1); 
            for h=2:numDataHybes %   h =2;
                % deal with missing data
                alignHybe = h-1;
                good = ~isnan(chainPts_cy3{fidChain}(alignHybe,1));
                while ~good && alignHybe > 1
                    alignHybe = alignHybe -1;
                    good = ~isnan(chainPts_cy3{fidChain}(alignHybe,1));
                end
                if ~isnan(chainPts_cy3{fidChain}(h,1))
                    xshift_h = chainPts_cy3{fidChain}(h,1) - chainPts_cy3{fidChain}(alignHybe,1);
                    yshift_h = chainPts_cy3{fidChain}(h,2) - chainPts_cy3{fidChain}(alignHybe,2);
                    zshift_h = chainPts_cy3{fidChain}(h,3) - chainPts_cy3{fidChain}(alignHybe,3);
                else
                    xshift_h = 0;
                    yshift_h = 0;
                    zshift_h = 0;
                end
                xshift = xshift + xshift_h;
                yshift = yshift + yshift_h;
                zshift = zshift + zshift_h;
                xshifts(h) = xshift;
                yshifts(h) = yshift;
                zshifts(h) = zshift;
            end

            % apply shifts
            shifts{s}(1:numDataHybes,:) = [xshifts,yshifts,zshifts];
            for h=1:numDataHybes %
                cntrlChain0(h,1) = cntrlChain0(h,1) - xshifts(h);
                cntrlChain0(h,2) = cntrlChain0(h,2) - yshifts(h);
                cntrlChain0(h,3) = cntrlChain0(h,3) - zshifts(h);

                cntrlChain(h,1) = cntrlChain(h,1) - xshifts(h);
                cntrlChain(h,2) = cntrlChain(h,2) - yshifts(h);
                cntrlChain(h,3) = cntrlChain(h,3) - zshifts(h);

                hIds = linkTable_cy3{n}.hybe == h;
                xIds = StringFind(linkTable_cy3{n}.Properties.VariableNames,{'x','xL','xU'},'exactly',true);
                yIds = StringFind(linkTable_cy3{n}.Properties.VariableNames,{'y','yL','yU'},'exactly',true);
                zIds = StringFind(linkTable_cy3{n}.Properties.VariableNames,{'z','zL','zU'},'exactly',true);
                linkTable_cy3{n}{hIds,xIds} = linkTable_cy3{n}{hIds,xIds} - xshifts(h);
                linkTable_cy3{n}{hIds,yIds} = linkTable_cy3{n}{hIds,yIds} - yshifts(h);
                linkTable_cy3{n}{hIds,zIds} = linkTable_cy3{n}{hIds,zIds} - zshifts(h); 
            end

            % apply to all cy5 chains
            for h=1:min([numHybes,numCy5Pts,numCy3Pts])
                for n=1:length(chainPts_cy5_align)
                    chainPts_cy5_align{n}(h,1) = chainPts_cy5_align{n}(h,1) - xshifts(h);
                    chainPts_cy5_align{n}(h,2) = chainPts_cy5_align{n}(h,2) - yshifts(h);
                    chainPts_cy5_align{n}(h,3) = chainPts_cy5_align{n}(h,3) - zshifts(h);

                    hIds = linkTable_cy5{n}.hybe == h;
                    xIds = StringFind(linkTable_cy5{n}.Properties.VariableNames,{'x','xL','xU'},'exactly',true);
                    yIds = StringFind(linkTable_cy5{n}.Properties.VariableNames,{'y','yL','yU'},'exactly',true);
                    zIds = StringFind(linkTable_cy5{n}.Properties.VariableNames,{'z','zL','zU'},'exactly',true);
                    linkTable_cy5{n}{hIds,xIds} = linkTable_cy5{n}{hIds,xIds} - xshifts(h);
                    linkTable_cy5{n}{hIds,yIds} = linkTable_cy5{n}{hIds,yIds} - yshifts(h);
                    linkTable_cy5{n}{hIds,zIds} = linkTable_cy5{n}{hIds,zIds} - zshifts(h); 
                end
            end

            disp(['Finished fiducial drift correction for spot ',num2str(s)]);
            disp(['saving data for spot ',num2str(s)]);
            
            try
                if height(linkTable_cy5{1}) > 1
                    tubeFig = figure(100); clf;
                    tubeFig.Visible = 'off';
                    PlotTube(linkTable_cy5{1}{:,1:3},'r',renderTubeRadius,'interpPts',10,'method',renderInterpMethod,'colormap',jet(numCy5Pts)); 
                    set(gca,'color','k'); hold on;
                    alpha .5;
                    freezeColors; 
                    PlotSpheres(chainPts_cy5_align{1},'r',renderBallRadius,'color',jet(numCy5Pts),'alpha',.55);
                    text(linkTable_cy5{1}{:,1}+renderBallRadius,linkTable_cy5{1}{:,2},linkTable_cy5{1}{:,3},num2str(linkTable_cy5{1}.hybe),'color','w');
                    axis equal;
                    imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_tubeFig'];
                    SaveFigure(tubeFig,'name',imName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);
                elseif verbose
                    warning(['spot ',num2str(s),' contains 1 of fewer cy5 points!']);
                end
             catch er
                error(['error rendering TubePlot for spot ',num2str(s)]); 
             end
            %% create maps and save data
            numCy5Chains = 1;%  length(chainPts_cy5_align);
            hic5Fig = figure(101); clf; 
            hic5Fig.Visible = 'off';
            for n = 1:numCy5Chains
                im = squareform(pdist(chainPts_cy5_align{n}(:,:)));
                cmap = GetColorMap('redToWhite',100); colormap(cmap);
                subplot(1,numCy5Chains,n); imagesc(im); colorbar; % caxis([0,3]);
            end
            imName = ['fov',num2str(fov,'%03d'),'_spot',num2str(s,'%03d'),'_hic5Fig'];
            SaveFigure(hic5Fig,'name',imName,'formats',{'png','fig'},'overwrite',true,'saveData',saveData);

            cy5Chains{s} = chainPts_cy5_align;
            cy3Chains{s} = cntrlChain;
            cy5LinkTables{s} = linkTable_cy5;
        catch er
            if verbose
                disp(er.message);
                disp(er.getReport);          
            end
            if stopOnError
                error(['error on line ',num2str(er.stack.line),': ',er.message]);
            else
                disp(['error on spot ',num2str(s),'!']);
            end
        end
    end
end  
    
%%
if numSpots > 10
    k = 0;
    M = NaN(numHybes,numHybes,1000);
    for s=1:numSpots
        for t=1:length(cy5Chains{s})
           m = squareform(pdist(cy5Chains{s}{t}(:,1:2)));
           [n,n] = size(m);
           k=k+1;
           M(1:n,1:n,k) = m;
        end
    end
    figure(6); clf; imagesc(nanmedian(M,3)); colorbar;
end
