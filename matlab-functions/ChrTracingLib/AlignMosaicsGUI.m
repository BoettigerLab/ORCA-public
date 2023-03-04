classdef AlignMosaicsGUI < StepGUIclass

% MosaicBuilder specific dynamic properties
    properties % editable parameters
        % GUI properties
        currStep = 'Load Mosaic1';
        
        % --- data properties
        chn1Data = []; % storing data in structures makes it easier to add more global properties to the  
        chn2Data = [];
        comboData =[];
%         % channel 1        
%         analysisFolder1 = '';
%         chnNames1 = {};
%         imageTiles1 = {};
%         uls1 = [];
%         mosaic1 =[];
%         % channel 2
%         analysisFolder2 ='';
%         chnNames2 = {};
%         imageTiles2 ={};
%         imTilesRot2 = {};
%         uls2 = [];
%         coarseAlign = struct(); 
%         comboMosaic = []; % 
%         uls2to1 = []; % coordinates of mosaic 2 panels in mosaic 1 
%         refULs = []; % aligned ULs, in mosaic1 coords. 
    end
    
    % MosaciBuilder specific static properties
    %   though all stepGUIs will need to define these properties
    properties % (Access = private) % 
       allStepNames = {'Load Mosaic1',...
                       'Load Mosaic2',...
                       'Coarse Align Mosaics',...
                       'Fine Align Mosaics',...
                       'Validate Alignment',...
                       'Merge Mosaics',...
                       'Save Mosaic Data',...
                       'Explore Mosaic'};
       stepDirections = {'Step 1: load mosaic 1 (typically DNA/ORCA data).  Validate that FOVs array properly into the mosaic.';
                        'Step 2: load mosaic 2 (typically RNA)';
                        'Step 3: coarse align mosaics';
                        'Step 4: fine align mosiacs';
                        'Step 5: validate fine alignment';
                        'Step 6: assemble mosaic';
                        'Step 7: save alignment data'}
      allFxnHandles = {@(inputs,pars) LoadMosaic1(inputs,pars),...
                @(inputs,pars) LoadMosaic2(inputs,pars),...
                @(inputs,pars) CoarseAlign(inputs,pars),...
                @(inputs,pars) FineAlign(inputs,pars),...
                @(inputs,pars) ValidateAlign(inputs,pars),...
                @(inputs,pars) AssembleMosaic(inputs,pars),...
                @(inputs,pars) SaveAlignment(inputs,pars),...
                @(inputs,pars) ExploreMosaic(inputs,pars)};
       stepParameters = {};
    end

    
    methods            
        function self = AlignMosaicsGUI
            % inherets from MosaicBuilderClass
            % GUI interface is inhereted 
            % button properties are inhereted 
            % Load Mosaic 1 (DNA, default fid only, 1 hyb only)
             self.gui_h.figure1.Name='Align Mosaics GUI'; 
            
            defaults = cell(0,3);
            defaults(end+1,:) = {'chnTableFile','string','Enter Path To DNA Experimental Table: C:\'};
            defaults(end+1,:) = {'analysisFolder','string','.\Analysis\'};
            defaults(end+1,:) = {'backgroundCorrect','freeType','none'}; % options: 'none', 'flatten', bkdImage (actual image data, not a string);
            defaults(end+1,:) = {'method',{'mean','sum','edgeBlur','last','first'},'edgeBlur'}; 
            defaults(end+1,:) = {'mosaicContrastLow','fraction',.1};
            defaults(end+1,:) = {'mosaicContrastHigh','fraction',.9999};
            defaults(end+1,:) = {'fovs','freeType',inf};
            defaults(end+1,:) = {'hybs','freeType',1};
            defaults(end+1,:) = {'dataType',{'fiducial','data','all'},'fiducial'};
            defaults(end+1,:) = {'scope',{'autoDetect','scope1','scope2','scope3','other'},'autoDetect'};  %
%             defaults(end+1,:) = {'transpose','boolean',true}; % listed in order of op
%             defaults(end+1,:) = {'fliplr','boolean',true}; % listed in order of op
%             defaults(end+1,:) = {'flipud','boolean',false}; % listed in order of op
            defaults(end+1,:) = {'trimBorder','integer',10};
            parsLoadMosaic1  = ParseVariableArguments([],defaults,'MosaicAnalyzer LoadMosaic1');
            % Load Mosaic 2 (RNA, default 1 copy fid, all hybs of data)
            defaults = cell(0,3);
            defaults(end+1,:) = {'chnTableFile','string','Enter Path To RNA Experimental Table: C:\'};
            defaults(end+1,:) = {'analysisFolder','string','.\Analysis\'};
            defaults(end+1,:) = {'backgroundCorrect','freeType','none'}; % options: 'none', 'flatten', bkdImage (actual image data, not a string);
            defaults(end+1,:) = {'method',{'mean','sum','edgeBlur','last','first'},'edgeBlur'}; 
            defaults(end+1,:) = {'mosaicContrastLow','fraction',.1};
            defaults(end+1,:) = {'mosaicContrastHigh','fraction',.9999};
            defaults(end+1,:) = {'fovs','freeType',inf};
            defaults(end+1,:) = {'hybs','freeType',inf};
            defaults(end+1,:) = {'dataType',{'fiducial','data','all'},'all'};
            defaults(end+1,:) = {'scope',{'autoDetect','scope1','scope2','scope3','other'},'autoDetect'};  %
%             defaults(end+1,:) = {'transpose','boolean',true}; % listed in order of op
%             defaults(end+1,:) = {'fliplr','boolean',true}; % listed in order of op
%             defaults(end+1,:) = {'flipud','boolean',false}; % listed in order of op
            defaults(end+1,:) = {'trimBorder','integer',10};
            parsLoadMosaic2  = ParseVariableArguments([],defaults,'MosaicAnalyzer LoadMosaic2');
            % CoarseAlign
            defaults = cell(0,3);
            defaults(end+1,:) = {'align',{'matchFirstUpperLeft','matchFirstLowerRight','matchFirstCenter'},'matchFirstUpperLeft'};
            defaults(end+1,:) = {'angles','freeType',-12:.5:12};
            defaults(end+1,:) = {'downsample','positive',10};
            defaults(end+1,:) = {'maxSize','positive',400};
            defaults(end+1,:) = {'alignContrastHigh','fraction',.999};
            defaults(end+1,:) = {'alignContrastLow','fraction',.9};
            defaults(end+1,:) = {'maxShift','integer',300};
            defaults(end+1,:) = {'showPlot','boolean',true};
            defaults(end+1,:) = {'verbose','boolean',true};
            defaults(end+1,:) = {'troubleshoot','boolean',false};
            parsCoarseAlign =  ParseVariableArguments([],defaults,'MosaicAnalyzer CoarseAlign');
            % FineAlign
            defaults = cell(0,3);          
            defaults(end+1,:) = {'interactive','boolean',true};
            defaults(end+1,:) = {'corrAlign','boolean',true};
            defaults(end+1,:) = {'padAlign','integer',150}; % pixels to add on either side of reference image
            defaults(end+1,:) = {'corrAlignHigh','boolean',.998};
            defaults(end+1,:) = {'corrAlignLow','boolean',.65};
            defaults(end+1,:) = {'angles','float',0}; % -10:1:10
            defaults(end+1,:) = {'scales','float',1}; % 0.9:0.01:1.10
            defaults(end+1,:) = {'showplot', 'boolean', true};
            defaults(end+1,:) = {'maxSize', 'positive', 400};
            defaults(end+1,:) = {'showExtraPlot','boolean',false};
            defaults(end+1,:) = {'fineBox', 'freeType', 250}; 
            defaults(end+1,:) = {'minFineImprovement', 'float', .1};
            parsFineAlign =  ParseVariableArguments([],defaults,'MosaicAnalyzer FineAlign');
            % ValidateAlign
            defaults = cell(0,3);  
            defaults(end+1,:) = {'redContrastLow','fraction',.1};
            defaults(end+1,:) = {'redContrastHigh','fraction',.9999};
            defaults(end+1,:) = {'cyanContrastLow','fraction',.1};
            defaults(end+1,:) = {'cyanContrastHigh','fraction',.9999};
            parsValidateAlign =  ParseVariableArguments([],defaults,'MosaicAnalyzer ValidateAlign');
            % AssembleMosaic
            defaults = cell(0,3);
            defaults(end+1,:) = {'verbose','boolean',true};
            parsMerge =  ParseVariableArguments([],defaults,'MosaicAnalyzer Merge');
            % SaveAlign
            defaults = cell(0,3);
            defaults(end+1,:) = {'saveName','string',''};
            defaults(end+1,:) = {'saveFolder','string',''};
            parsSaveAlign =  ParseVariableArguments([],defaults,'MosaicAnalyzer SaveAlign');
            % ExploreMosaic
            defaults = cell(0,3);
            defaults(end+1,:) = {'verbose','boolean',true};
            parsExplore =  ParseVariableArguments([],defaults,'MosaicAnalyzer Explorer');
            % stack parameters
            self.stepParameters = {parsLoadMosaic1,...
                            parsLoadMosaic2,...
                             parsCoarseAlign,...
                             parsFineAlign,...
                             parsValidateAlign,...
                             parsMerge,...
                             parsSaveAlign,...
                             parsExplore};
        end
        


        function self = LoadMosaic1(self,pars)
            % --- read experiment table
            eTable1_xls = pars.chnTableFile;
            eTable1 = readtable(eTable1_xls);
            % --- extract analysis folder [make this a sub-function)
            if strcmp(pars.analysisFolder(1),'.') % current directory 
                saveName = fileparts(pars.chnTableFile);
                saveName = [saveName,pars.analysisFolder(2:end)];
            else
                saveName = pars.analysisFolder;
            end
            if ~strcmp(saveName(end),filesep)
                saveName = [saveName,filesep];
            end
            self.chn1Data.analysisFolder = saveName;
            % --- interpet channel names
            if pars.hybs == 1 && strcmp(pars.dataType,'fiducial')
                self.chn1Data.chnNames = {'fiducial1'}; % default load condition
            else
                [~,self.chn1Data.chnNames] = TileLabelsFromEtable(eTable1,'style','names');      % Record Names
                % otherwise take names from table
            end 
            % load names of all data requested
            [~,~,~,daxNames1] =  LoadDaxFromEtable(eTable1_xls,...
                            'parameters',pars,'saveFolder',self.chn1Data.analysisFolder,...
                            'hybNumber',pars.hybs,'fov',pars.fovs);
            % update pars.fovs
            if isinf(pars.fovs)
                self.stepParameters{1}.fovs = GetFOVnums(daxNames1);
            end
            % create hybe-registered image tiles of all data
            [self.chn1Data.imageTiles,stageXY1,mosPars] =LoadRegDaxAsMosaic(daxNames1,...        % Image tiles in memory
                'parameters',pars,'saveFolder',self.chn1Data.analysisFolder);
            self.stepParameters{1}.scope = mosPars.scope; % save autodetected scope name
            xyShifts1 = readtable([self.chn1Data.analysisFolder,'regFOV.csv']);     % 
            xyShifts1 = [xyShifts1.xShift,xyShifts1.yShift]; % convert table to matrix  
            self.chn1Data.uls = stageXY1 + xyShifts1;                                % Corrected image positions in memory
            mos1 = MosaicViewerRender(self.chn1Data.imageTiles,self.chn1Data.uls,...
                'parameters',pars,'downsample',4,...
                'fliplr',false,'flipud',false,'transpose',false);
            figure(1); clf; 
            imagesc(IncreaseContrast(mos1,'high',.9999,'low',.3)); 
            colormap(gray); 
            set(self.gui_h.TextDir,'String','Step 1: Load Mosaic 1 complete. Select Next Step.');
        end
        
        function self = LoadMosaic2(self,pars)
            % load experiment table
            eTable2_xls = pars.chnTableFile;
            eTable2 = readtable(eTable2_xls);
            % parse out analysis folder (this is clutzy)
            if strcmp(pars.analysisFolder(1),'.') % current directory 
                saveName = fileparts(pars.chnTableFile);
                saveName = [saveName,pars.analysisFolder(2:end)];
            else
                saveName = pars.analysisFolder;
            end
            if ~strcmp(saveName(end),filesep)
                saveName = [saveName,filesep];
            end
            self.chn2Data.analysisFolder = saveName;
            % get channel names from experiment table
            [~,chnNames] = TileLabelsFromEtable(eTable2,'style','names');
            self.chn2Data.chnNames = cat(1,chnNames,'fiducial2');
            % load names of all data requested
            [~,~,~,daxNames2] =  LoadDaxFromEtable(eTable2_xls,...
                            'parameters',pars,'saveFolder',self.chn2Data.analysisFolder,...
                            'hybNumber',pars.hybs,'fov',pars.fovs);
            % create hybe-registered image tiles of all data
            [self.chn2Data.imageTiles,stageXY2] =LoadRegDaxAsMosaic(daxNames2,...
                'parameters',pars,'saveFolder',self.chn2Data.analysisFolder);
            % fetch xy-shifts from CSV if not in memory
            xyShifts2 = readtable([self.chn2Data.analysisFolder,'regFOV.csv']);
            xyShifts2 = [xyShifts2.xShift,xyShifts2.yShift]; % convert table to matrix  NO LONGER NEEDED 
            self.chn2Data.uls = stageXY2 + xyShifts2;
            % plot the mosaic
              mos2 = MosaicViewerRender(self.chn2Data.imageTiles(:,end),self.chn2Data.uls,...
                'parameters',pars,'downsample',4,...
                'fliplr',false,'flipud',false,'transpose',false);            
             figure(3); clf; imagesc(mos2);  
             imagesc(IncreaseContrast(mos2,'high',.9999,'low',.3)); 
            colormap(gray); 
            set(self.gui_h.TextDir,'String','Step 2: Load Mosaic 2 complete. Select Next Step.');
        end
    
        function CoarseAlign(self,pars)
            % -- check for existing coarse alignment folder
            coarseAlignFile = [self.chn2Data.analysisFolder,'coarseAlign.csv'];
            answer = 'Overwrite';
            if exist(coarseAlignFile,'file')
                  answer = questdlg(['found existing coarse alignment file'], ... % question
                    'Prompt ', ...  % pop-up label
                    'Load','Overwrite','Cancel','Load'); % op1 op2 default
            end 
            switch answer
                case 'Cancel'
                    disp('run step canceled by user');
                case 'Load'
                    coarseAlignTable = readtable(coarseAlignFile);                
                    self.comboData.coarseAlign = table2struct(coarseAlignTable);
                    if isfield(self.comboData.coarseAlign,'center_1') % convert center_1,center_2 back to array
                        self.comboData.coarseAlign.center = [self.comboData.coarseAlign.center_1,self.comboData.coarseAlign.center_2];
                    end
                case 'Overwrite'
                    % ---- Create downsampled mosaics
                    mos1 = TilesToMosaic(self.chn1Data.imageTiles,self.chn1Data.uls,'padMosaic',1,'parameters',pars);
                    mos2 = TilesToMosaic(self.chn2Data.imageTiles(:,end),self.chn2Data.uls,'padMosaic',1,'parameters',pars);                    
                    % --- Merge and contrast
                    imO = MergeImages(mos1,mos2,'align',pars.align); % pad appropriately to make same size.
                    [h1,w1] = size(mos1);
                    [h2,w2] = size(mos2);
                    uls_Overlay = self.chn2Data.uls;
                    if strcmp(pars.align,'matchFirstCenter')
                         uls_Overlay(:,1) = uls_Overlay(:,1) - ceil((w2-w1)/2*pars.downsample);
                         uls_Overlay(:,2) = uls_Overlay(:,2) - ceil((h2-h1)/2*pars.downsample);
                    elseif strcmp(pars.align,'matchFirstLowerRight')
                         uls_Overlay(:,1) =  uls_Overlay(:,1) - (w2-w1)*pars.downsample;
                         uls_Overlay(:,2) = uls_Overlay(:,2) - (h2-h1)*pars.downsample;
                    end
                    self.chn2Data.uls_Overlay = uls_Overlay;
                    imD = IncreaseContrast(imO(:,:,1),'high',pars.alignContrastHigh,'low',pars.alignContrastLow );
                    imR = IncreaseContrast(imO(:,:,2),'high',pars.alignContrastHigh,'low',pars.alignContrastLow );
                    if pars.troubleshoot
                        figure(100); clf; % for troubleshooting
                        subplot(1,3,1); imagesc(imD); 
                        subplot(1,3,2); imagesc(imR);  
                        subplot(1,3,3); Ncolor(cat(3,imD,imR)); 
                        colormap(gray);
                        title('can you see how these match');
                        disp('paused, press any key to continue, or ctrl-c to abort');
                        pause();
                    end
                    % ----- compute and show the alignment
                    figure(4); clf; pause(.1); 
                    disp('computing coarse alignment...');
                    alignValues = CorrAlignFast(imD,imR,...
                                            'parameters',pars);
                    imR2 = ApplyReg(imR,alignValues);
                    figure(3); clf; Ncolor(cat(3,imD,imR2));
        %             imD2 = ApplyReg(imD,alignValues,'invert',true);
        %             figure(4); clf; Ncolor(cat(3,imD2,imR));
                    % save variables
                    alignValues.xshift = pars.downsample*alignValues.xshift; % (inconsistent capitalization, to clean up sometime)  
                    alignValues.yshift = pars.downsample*alignValues.yshift;
                    alignValues.theta = alignValues.theta;
                    alignValues.center = pars.downsample*round(size(imD)/2);
                    self.comboData.coarseAlign = alignValues;
                    disp('coarse alignment complete'); % move this
                    %--- write table
                    coarseAlignTable = struct2table(alignValues);
                    writetable(coarseAlignTable,coarseAlignFile,'delimiter','\t');                   
            end
            set(self.gui_h.TextDir,'String','Step 3: coarse alignment complete. Select Next Step.');
        end
    
        function FineAlign(self,pars)
            % (used to be done after crop. Trying more general version,
            % lets do this on the mosaic as whole).
            [numFOV2, numChns2] = size(self.chn2Data.imageTiles);           
            fovToDock = 1:numFOV2;
            tableName = [self.chn2Data.analysisFolder,'uls2to1.csv'];
            if exist(tableName,'file')==2  % matlabs number system here is annoying, we gave it a flag, it should return a boolean 
                uls2to1_table = readtable(tableName);
                hasData =  ~isnan(uls2to1_table.x)';
                % prompt, found data
                if sum(hasData) > 0
                    disp('found existing docking positions for FOVs:'); disp(find(hasData));
                    answer = questdlg(['found ',num2str(sum(hasData)),' existing docking positions'], ... % question
                        'Prompt ', ...  % pop-up label
                        'Skip','Overwrite','Cancel','Skip'); % op1 op2 default
                else
                    answer = 'New';
                end
                switch answer
                    case 'Skip'
                        fovToDock(hasData) = [];
                        disp('using existing docking positions for FOVs:');
                        disp(find(hasData));
                    case 'Overwrite'
                        disp('will now overwrite docking positions for FOVs:');
                        disp(find(hasData));
                    case 'Cancel'
                        warning('canceled by user');
                        return
                    case 'New'
                        disp(['processing data for ',num2str(numFOV2),' FOVs']);
                end               
            else % write a blank table
                x = nan(numFOV2,1);
                y = nan(numFOV2,1);
                theta = nan(numFOV2,1); % 
                scale = ones(numFOV2,1);
                theta2 = zeros(numFOV2,1);
                uls2to1_table = table(x,y,theta,scale,theta2);
                writetable(uls2to1_table,tableName);
                disp(['wrote ',tableName]);
            end
            % apply rotation translation to RNA data. 
            theta = self.comboData.coarseAlign.theta;
            disp('applying rotation to imTiles for Mosaic2');
            [uls2_align,imTiles2] = RotateTranslatePoints(...
                self.chn2Data.uls_Overlay,self.comboData.coarseAlign,...
                                'center',self.comboData.coarseAlign.center,...
                                'imageTiles',self.chn2Data.imageTiles,...
                                'invert',false,...
                                'opOrder','TR');
            disp('creating full resolution image for Mosaic1');
            [self.chn1Data.mosaic,self.comboData.refULs] = TilesToMosaic(self.chn1Data.imageTiles,self.chn1Data.uls,'padMosaic',1);
            [numFOV2,~] = size(imTiles2);
            alignValues = cell(numFOV2,1);
            disp('docking tiles of Mosaic2 into Mosaic1 at full resolution');
            for f = fovToDock
                % loop through FOV and dock mosaic2 panels (RNA) into
                % fully assembled mosaic1 (DNA).
                % Note, currently does not allow skipping. I think this is
                % okay, we're not laying the files onto one another, just
                % onto the DNA mosaic.
                % figure(12); clf; imagesc(IncreaseContrast(imTiles2{f,end},'high',.999));
                [alignValues{f},~,valid] = InteractiveAlign(self.chn1Data.mosaic,imTiles2{f,end},... 
                            'parameters',pars,...
                            'fineCenter',pars.padAlign+size(imTiles2{f,end})/2,...
                            'startPosition',uls2_align(f,:),...
                            'imName',num2str(f));        
                ul = ScaleRotateTranslatePoints(uls2_align(f,:),alignValues{f},'imSize',size(imTiles2{f,end}));       
                if (valid == 1) || (valid == 4)
                    uls2to1_table{f,:} = [ul(1),ul(2),theta,alignValues{f}.rescale,alignValues{f}.theta];
                    writetable(uls2to1_table,tableName); % will overwrite existing table
                else
                    disp(['flagged FOV ',num2str(f) ' to align later']);
                end  
            end 
            self.chn2Data.imTilesRot = imTiles2;
            set(self.gui_h.TextDir,'String','Step 4: Fine alignment complete. Select Next Step.');
        end
        
        function ValidateAlign(self,pars)
            try
            [numFOV2, numChns2] = size(self.chn2Data.imageTiles);
            [h_m,w_m] = size(self.chn1Data.mosaic);
            tableName = [self.chn2Data.analysisFolder,'uls2to1.csv'];
            uls2to1_table = readtable(tableName); 
            % check table is complete
            missing = isnan(uls2to1_table.x);
            if any(missing)
                wasSkipped = find(missing);
                cprintf([1 0 0],['Error, unable to find data for FOVs ',num2str(wasSkipped')]);
                cprintf([1 0 0],['Return to previous step and correct these FOVs']);
                error('Not all tiles have been placed yet!');
            end
            uls2_align2 = [uls2to1_table.x,uls2to1_table.y]; 
            % update tiles; 
            disp('updating tile rotation and scaling, please wait...');
            imTiles2 = cell(numFOV2,1); % a new variable outside of 'self', so that 'self' does not get parfor broadcast
            for n=1:numFOV2 
                im = self.chn2Data.imageTiles{n,numChns2};
                alignValues.rescale = uls2to1_table.scale(n);
                alignValues.theta = uls2to1_table.theta(n);
                alignValues.rescale2 = 1;
                alignValues.theta2 = uls2to1_table.theta2(n);
                im = ScaleRotateIm(im,alignValues);
                imTiles2{n} = im;
            end
            mosaicRNA = TilesToMosaic(imTiles2,uls2_align2,'mosaicSize',[h_m,w_m]);
            mosaicRNA = IncreaseContrast(mosaicRNA,'high',pars.cyanContrastHigh,'low',pars.cyanContrastLow);
            mosaicDNA = IncreaseContrast(self.chn1Data.mosaic,'high',pars.redContrastHigh,'low',pars.redContrastLow);
           overlayMosaic = cat(3,mosaicDNA,mosaicRNA); 
           figure(); clf; Ncolor(overlayMosaic);  title('mosaic1-DNA-red mosaic2-RNA-cyan');
           % NcolorApp(overlayMosaic,'names',{'mosaic1-DNA','mosaic2-RNA'});
           set(self.gui_h.TextDir,'String','Zoom in and validate alignment');
            catch er
                cprintf([1 0 0],er.getReport);
                disp('place debug here');
                error(er.message);
            end
        end
        
        
        function AssembleMosaic(self,pars)
            % Use the assembled mosaic1 and the computed rotation +
            % translation 
            % this should probably call a stand-alone function, as this can
            % be used to load and align also the 3D dax, since it works
            % directly on imageTiles. Of course the tile blending will have
            % to be computed anew, which is fitting once it's not
            % max-project.
            % 
            [numFOV2, numChns2] = size(self.chn2Data.imageTiles);
            [h_m,w_m] = size(self.chn1Data.mosaic);
            tableName = [self.chn2Data.analysisFolder,'uls2to1.csv'];
            uls2to1_table = readtable(tableName); 
            uls2_align2 = [uls2to1_table.x,uls2to1_table.y]; 
            % update tiles; 
            disp('updating tile rotation and scaling, please wait...');
            imTiles2 = cell(numFOV2,numChns2); % a new variable outside of 'self', so that 'self' does not get parfor broadcast
            for n=1:numFOV2
                for c=1:numChns2
                    im = self.chn2Data.imageTiles{n,c};
                    alignValues.rescale = uls2to1_table.scale(n);
                    alignValues.theta = uls2to1_table.theta(n);
                    alignValues.rescale2 = 1;
                    alignValues.theta2 = uls2to1_table.theta2(n);
                    im = ScaleRotateIm(im,alignValues);
                    imTiles2{n,c} = im;
                end
            end
            mosaicStack = cell(numChns2,1); % enables use of parfor 
            disp('assembling combined mosaic, this will take time...');
            parfor c=1:numChns2  % c = 6
                disp(['processing channel ',num2str(c), ' of ',num2str(numChns2)]);
                mosaicStack{c} = TilesToMosaic(imTiles2(:,c),uls2_align2,'mosaicSize',[h_m,w_m]);
            end
           self.comboData.comboMosaic = cat(3,self.chn1Data.mosaic,mosaicStack{:}); 
           chnName = cat(1,self.chn1Data.chnNames,self.chn2Data.chnNames);
           NcolorApp(self.comboData.comboMosaic,'names',chnName);
           set(self.gui_h.TextDir,'String','Step 5: Mosaic assembly complete. Select Next Step.');
        end

        function ExploreMosaic(self,pars)
            chnName = cat(1,self.chn1Data.chnNames,self.chn2Data.chnNames);
            NcolorApp(self.comboData.comboMosaic,'names',chnName);
        end
            
        function SaveAlignment(self,pars) % SaveAlignment
            % default save root
            if isempty(pars.saveName)
                saveRoot = 'ComboMosaic';
            else
                saveRoot = pars.saveName;
            end
            % default save folder
            if isempty(pars.saveFolder)
                saveFolder = self.chn2Data.analysisFolder;
            else
                saveFolder = pars.saveFolder;
            end
            tifName = [saveFolder,saveRoot,'.tif'];
            % check if files exist
            [~,~,nChns] = size(self.comboData.comboMosaic);
            skip = false;
            if exist(tifName,'file')~=0
                answer = questdlg([tifName ' exists, overwrite?'], ... % question
                    'Prompt ', ...  % pop-up label
                    'Yes','No','No'); % op1 op2 default
                switch answer
                    case 'Yes'
                        delete(tifName);
                        skip = false;
                    case 'No'
                        skip = true;
                end         
            end
            if ~skip
                % write table of channel names
                tifLabels = [saveFolder,saveRoot,'Labels.txt'];
                chnName = cat(1,self.chn1Data.chnNames,self.chn2Data.chnNames);
                writetable(table(chnName),tifLabels);
                % write table of uls
                ulTable = [saveFolder,saveRoot,'tileULs.csv'];
                x = self.comboData.refULs(:,1);
                y = self.comboData.refULs(:,2);
                [w,h] = size(self.chn1Data.imageTiles{1});
                w= w*ones(length(x),1);
                h= h*ones(length(x),1);
                trim = self.stepParameters{1}.trimBorder*ones(length(x),1);
%                 fliplr= self.stepParameters{1}.fliplr*ones(length(x),1);
%                 flipud= self.stepParameters{1}.flipud*ones(length(x),1);
%                 transpose= self.stepParameters{1}.transpose*ones(length(x),1);
                scope= repmat(self.stepParameters{1}.scope,length(x),1);
                fov = Column(self.stepParameters{1}.fovs);
                writetable(table(x,y,w,h,trim,scope,fov),ulTable);
                try
                    for c=1:nChns
                        mosaicOut= self.comboData.comboMosaic(:,:,c); 
                        tifSaveName = regexprep(tifName,'.tif',['_c',num2str(c,'%03d'),'.tif']);
                        imwrite(mosaicOut,tifSaveName);
                        % imwrite(mosaicOut,tifName,'WriteMode','append');
                    end
                catch er
                    warning(er.getReport);
                    disp('something wrong');
                end
                disp(['finished writing ',tifName]);                     
            end
        end
        
    end
end