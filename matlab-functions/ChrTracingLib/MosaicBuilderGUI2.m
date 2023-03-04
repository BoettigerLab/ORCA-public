classdef MosaicBuilderGUI2 < StepGUIclass
% Update Log
% An update of MosaicBuilderClass, to build as a child of StepGUIclass
% instead of being a direct decendent of handles. 
    
% MosaicBuilder specific dynamic properties
    properties % editable parameters
        currStep = 'Load Data' % overload
        saveFolder = '';
        daxInHybFov = {};
        maxNameInHybFov = {};
        imageTiles = {};
        stageXY = [];
        mosaicGray = [];
        mosaicImage = []; % all channels
        xyShifts =[];
        ulOut = [];
        chnNames = {};
    end
    
    % MosaciBuilder specific static properties
    %   though all stepGUIs will need to define these properties
    properties % (Access = private) % 
       allStepNames = {'Load Data',...
                       'Align Hybs',...
                       'Validate Mosaic',...
                       'Align FOV',...
                       'Explore Mosaic',...
                       'Save Tif'};
       stepDirections = {['Step 1: Create max projections of all FOVs. For ORCA DNA data, set parameter "dataType"="fiducial" and "hybs=1" to accelerate.',...
                        'For data such as multiplex RNA labeling, set parameter "dataType"="all" and "hybs"=inf.'];
                        'Step 2: Correct hyb-to-hyb drift across all requested hybs.';
                        'Step 3: Load mosaic tiles for inspection and prepare for fov-to-fov drift correction';
                        'Step 4: Correct fov-to-fov drift across all FOV. Requires always running step 3 first.';
                        'Step 5: Assemble and Explore Mosaic (Optional).';
                        'Step 6: Save Tif (Optional). Requires step 5 to have run.';
                        'Step 7: Browse assembled mosaic (requires step 5 to have run)'}
       stepSuccess = {'Max projections and experiment table loaded.',...
                      'Drift correction complete.',...
                      'Validation complete.',...
                      'FOV-to-FOV drift corrected.'};
       allFxnHandles = {@(inputs,pars) LoadData(inputs,pars),...
                        @(inputs,pars) AlignHybs(inputs,pars),...
                        @(inputs,pars) ValidateMosaic(inputs,pars),...
                        @(inputs,pars) AlignFOV(inputs,pars),...
                        @(inputs,pars) ExploreMosaic(inputs,pars),...
                        @(inputs,pars) SaveTif(inputs,pars),...
                        @(inputs,pars) BrowseMosaic(inputs,pars)};
       stepParameters;
    end

    
    methods            
        function self = BuildMosaicsGUI
            % The following are inhereted from StepGUIclass
            % inherets self.gui_h = guihandles(MultistepAnalysisGUI);
            % inherets: "ButtonRun"
            %           "ButtonNext"
            %           "ButtonBack"
            %           "ButtonPars"
            % inherets: Close_fcn
          
            self.gui_h.figure1.Name='Build Mosaics GUI'; % GUI name;
            
            % setup default parameters
            % Load Data         % S:\2018-05-15_WT_2-4hrs(HoxRNA)_(3kbDNA)\RNA_Expt\rnaTableFile.xlsx
            defaults = cell(0,3); % M:\2019-06-04_Rad21_TM3TwiGFP_18-20hrs_(HoxRNA)_(BXC3kbEven)\RNA_Expt\rnaTable_3color.xlsx
            defaults(end+1,:) = {'chnTableFile','string','M:\2019-06-04_Rad21_TM3TwiGFP_18-20hrs_(HoxRNA)_(BXC3kbEven)\RNA_Expt\rnaTable_3color.xlsx'};
            defaults(end+1,:) = {'analysisFolder','string','.\Analysis\'};
            defaults(end+1,:) = {'dataType',{'fiducial','data','all'},'all'};
            defaults(end+1,:) = {'hybNumber', 'freeType', inf}; 
            defaults(end+1,:) = {'fov','freeType',inf};
            parsLoadData  = ParseVariableArguments([],defaults,'MosaicBuilder LoadData');
            % Align Hybs
            defaults = cell(0,3);
            defaults(end+1,:) = {'alignContrastLow', 'fraction', .7}; % low image threshold for contrast balance prior to coarse alignment
            defaults(end+1,:) = {'alignContrastHigh', 'fraction', .9995}; % high threshold  for contrast balance prior to coarse alignment
            defaults(end+1,:) = {'refHybe','integer',1}; % hybe to use to start alignment
            defaults(end+1,:) = {'maxSize', 'positive', 400}; % rescale all images to this size for alignment
            defaults(end+1,:) = {'fineBox', 'freeType', 200};  % perform fine scale alignment using a box of this size around the brightest point.
            defaults(end+1,:) = {'fineUpsample', 'positive', 1};  
            defaults(end+1,:) = {'maxShift', 'freeType', inf};
            defaults(end+1,:) = {'gradMax', 'boolean', true};
            defaults(end+1,:) = {'minGrad', 'float', -inf};
            defaults(end+1,:) = {'fineMaxShift', 'nonnegative', 10};
            defaults(end+1,:) = {'fineCenter','array',[0,0]};
            defaults(end+1,:) = {'showplot', 'boolean', true};
            defaults(end+1,:) = {'fastDisplay', 'boolean', true};
            defaults(end+1,:) = {'displayWidth', 'integer', 500};
            defaults(end+1,:) = {'showExtraPlot', 'boolean', false};
            defaults(end+1,:) = {'minFineImprovement', 'float', 0}; 
            defaults(end+1,:) = {'showCorrAlign', 'boolean', true};
            parsAlignHybs  = ParseVariableArguments([],defaults,'MosaicBuilder AlignHybs');
            % Validate Mosaic
            defaults = cell(0,3);
            defaults(end+1,:) = {'backgroundCorrect','freeType','none'}; % options: 'none', 'flatten', bkdImage (actual image data, not a string);
            defaults(end+1,:) = {'method',{'mean','sum','edgeBlur','last','first'},'edgeBlur'}; 
            defaults(end+1,:) = {'mosaicContrastLow','fraction',.1};
            defaults(end+1,:) = {'mosaicContrastHigh','fraction',.9999};    % for MosaicViewerApp
            defaults(end+1,:) = {'fovs','freeType',inf};
            defaults(end+1,:) = {'hybs','freeType',1};
            defaults(end+1,:) = {'dataType',{'fiducial','data','all'},'fiducial'}; 
            defaults(end+1,:) = {'scope',{'autoDetect','scope1','scope2','scope3','other'},'autoDetect'};  %
            defaults(end+1,:) = {'positionsFile','string',''}; 
            defaults(end+1,:) = {'trimBorder','integer',10};
%             defaults(end+1,:) = {'transpose','boolean',true}; % listed in order of op
%             defaults(end+1,:) = {'fliplr','boolean',true}; % listed in order of op
%             defaults(end+1,:) = {'flipud','boolean',false}; % listed in order of op. Tr, Lr, Ud: Scope1=T T F. Scope2=T F T. 
            parsValidateMosaic  = ParseVariableArguments([],defaults,'MosaicBuilder ValidateMosaic');
            % AlignFOV
            defaults = cell(0,3);
            defaults(end+1,:) = {'interactive','boolean',true};
            defaults(end+1,:) = {'corrAlign','boolean',true};
            defaults(end+1,:) = {'gradMax', 'boolean', true};
            defaults(end+1,:) = {'showplot','boolean',false};
            defaults(end+1,:) = {'method',{'mean','sum','edgeBlur','last','first'},'edgeBlur'}; 
            defaults(end+1,:) = {'corrAlignLow','fraction',.8};
            defaults(end+1,:) = {'corrAlignHigh','fraction',.999}; % for corrAlign
            defaults(end+1,:) = {'mosaicContrastLow','fraction',.1};
            defaults(end+1,:) = {'mosaicContrastHigh','fraction',.9999};    % for MosaicViewerApp
            defaults(end+1,:) = {'fovs','freeType',inf};
            defaults(end+1,:) = {'maxShift','freeType',.2};
            parsAlignFOV  = ParseVariableArguments([],defaults,'MosaicBuilder AlignFOV');
            % Explore Mosaic
            defaults = cell(0,3);
            defaults(end+1,:) = {'backgroundCorrect','freeType','none'}; % options: 'none', 'flatten', bkdImage (actual image data, not a string);
            defaults(end+1,:) = {'method',{'mean','sum','edgeBlur','last','first'},'edgeBlur'}; 
            defaults(end+1,:) = {'mosaicContrastLow','fraction',.1};
            defaults(end+1,:) = {'mosaicContrastHigh','fraction',.9999};
            defaults(end+1,:) = {'fovs','freeType',inf};
            defaults(end+1,:) = {'hybs','freeType',inf};
            defaults(end+1,:) = {'dataType',{'fiducial','data','all'},'all'};
            defaults(end+1,:) = {'scope',{'autoDetect','scope1','scope2','scope3','other'},'autoDetect'};  %
            defaults(end+1,:) = {'positionsFile','string',''}; 
            defaults(end+1,:) = {'trimBorder','integer',10};
%             defaults(end+1,:) = {'transpose','boolean',true};
%             defaults(end+1,:) = {'fliplr','boolean',true};
%             defaults(end+1,:) = {'flipud','boolean',false};
            parsExploreMosaic  = ParseVariableArguments([],defaults,'MosaicBuilder ValidateMosaic');
            % Save Tif
            defaults = cell(0,3);
            defaults(end+1,:) = {'saveName','string','Mosaic'};
            defaults(end+1,:) = {'saveFolder','string',''};
            parsSaveTif  = ParseVariableArguments([],defaults,'MosaicBuilder ValidateMosaic');
            
            % stack
            self.stepParameters = {parsLoadData,...
                             parsAlignHybs,...
                             parsValidateMosaic,...
                             parsAlignFOV,...
                             parsExploreMosaic,...
                             parsSaveTif};
                         
            set(self.gui_h.TextDir,'String',self.stepDirections{1});
          
        end     
    end
    % ==== 
    methods   
% ========================== STEP 1:  Load Data ================================ % 
        function self = LoadData(self,pars)
            % get analysis folder
            if strcmp(pars.analysisFolder(1),'.') % current directory 
                saveName = fileparts(pars.chnTableFile);
                saveName = [saveName,pars.analysisFolder(2:end)];
            else
                saveName = pars.analysisFolder;
            end
            if ~strcmp(saveName(end),filesep)
                saveName = [saveName,filesep];
            end
            SetFigureSavePath(saveName,'makeDir',true);
            self.saveFolder = saveName;
           [self.daxInHybFov,~,~,...
            self.maxNameInHybFov] =  LoadDaxFromEtable(pars.chnTableFile,...
                                        'parameters',pars,...
                                        'saveFolder',self.saveFolder);
            if ~isempty(self.maxNameInHybFov)
               set(self.gui_h.TextDir,'String','data loaded, select next step.'); 
            end
            % update step 5 parameters to match
            self.stepParameters{5}.hybs = pars.hybNumber;
            self.stepParameters{5}.dataType = pars.dataType;
            self.stepParameters{5}.fovs = pars.fov;
        end
% ========================== STEP 2:  Align Hybs ================================ % 
        % ---- Align Hybs 
        function self = AlignHybs(self,pars)  
            fidChn = 1;
            [nHybs,nFov] = size(self.daxInHybFov(:,:,fidChn)); %#ok<ASGLU>
            if nHybs > 1
                ChrTracer3p2_FixGlobalDrift(self.daxInHybFov(:,:,fidChn),...
                    'saveFolder',self.saveFolder,'parameters',pars);
                 set(self.gui_h.TextDir,'String','Align data complete.');
            else
                textOut = 'no need to align single hyb, advance to next step';
                disp(textOut);
                set(self.gui_h.TextDir,'String',textOut); 
            end
        end
% ========================== STEP 3: Validate Mosaic ================================ % 
        function self = ValidateMosaic(self,pars)
            [self.imageTiles,self.stageXY,scopePars] = LoadRegDaxAsMosaic(self.maxNameInHybFov,'parameters',pars);
            % display results
            % figure(1); clf; imagesc(IncreaseContrast(self.imageTiles{1},'high',pars.mosaicContrastHigh));
            mva = MosaicViewerApp(self.imageTiles,self.stageXY,...
                'parameters',pars,...   
                'displaySize',inf,...
                'displayScale',.5); %   overrides displaySize
            disp('select "Done" when finished to record values and resume MosaicBuilder'); 
            uiwait(); % wait until mva is closed
            % in future, we will allow the user to save these manual shifts  
            manualShifts = evalin('base','manualShifts'); % pull from base workspace
            % in future, we will allow the user to save these contrast options   
            adjustedTiles = evalin('base','adjustedTiles');    
            % propegate parameters
            self.stepParameters{3}.scope = scopePars.scope;
            self.stepParameters{5}.scope = scopePars.scope;
            self.stepParameters{5}.backgroundCorrect = pars.backgroundCorrect;
            % save parameters
            savePars = struct2table(pars,'AsArray',true);
            writetable(savePars,[self.saveFolder,'ParsValidateMosaic.csv']);
        end

        % === AlignFOV step  4
        function self = AlignFOV(self,pars) 
            disp('if you wish to skip this step, abort (Ctrl+C), and set parameters CorrAlign and Interactive to "false"');
            % parse fov request (must do here, not passed to
            % MosaicViewerRender, which is also the aligner. 
            if isinf(pars.fovs)
                fovs = 1:length(self.imageTiles);
            else
                fovs=pars.fovs;
            end
            [self.mosaicGray,~,self.xyShifts] = ... 
                        MosaicViewerRender(self.imageTiles(fovs),...
                                           self.stageXY(fovs,:),...
                                           'parameters',pars,...
                                           'saveFolder',self.saveFolder);
           textDir = 'Optional: select individual tiles to validate alignment';
           set(self.gui_h.TextDir,'String',textDir);
           disp(textDir);
            MosaicViewerApp(self.imageTiles,self.stageXY+self.xyShifts,...
                'parameters',pars,...
                'displaySize',inf,...
                'displayScale',.5); % These parameters are already taken care of in DaxToImageTiles, and we don't want to apply these transforms twice  
        end
        % === Explore Mosaic  Step 5
        function self = ExploreMosaic(self,pars)
            % get eTable from step 1 pars
            eTableXls = self.stepParameters{1}.chnTableFile;
            [~,tempNames] = TileLabelsFromEtable(readtable(eTableXls),'style','names');
            if strcmp(pars.dataType,'all')
                self.chnNames = cat(1,tempNames,'fiducial');
            elseif strcmp(pars.dataType,'data')
                 self.chnNames = tempNames;
            elseif strcmp(pars.dataType,'fiducial')
                 self.chnNames ={'fiducial'};
            else % fiduical
                error([pars.dataType,' not recognized dataType']);
            end
            % load names of all data requested
            [~,~,~,self.maxNameInHybFov] =  LoadDaxFromEtable(eTableXls,...
                            'parameters',pars,'saveFolder',self.saveFolder,...
                            'hybNumber',pars.hybs,'fov',pars.fovs,...
                            'readDax',false,'verbose',false);
                        
            % update pars.fovs with fov nums for export
            if isinf(pars.fovs)  % 
                % self.stepParameters{5}.fovs = GetFOVnums(daxNames1);
                self.stepParameters{5}.fovs = GetFOVnums(self.maxNameInHybFov(1,:,1)); 
            end
            % create hybe-registered image tiles of all data
            [self.imageTiles,self.stageXY] =LoadRegDaxAsMosaic(self.maxNameInHybFov,...
                'parameters',pars,'saveFolder',self.saveFolder);
            % fetch xy-shifts from CSV if not in memory
            if isempty(self.xyShifts)
                if exist([self.saveFolder,'regFOV.csv'],'file')
                    self.xyShifts = LoadRegFov([self.saveFolder,'regFOV.csv']);
                    if ~isinf(pars.fovs)
                        self.xyShifts = self.xyShifts(pars.fovs,:);
                    end
                else
                    self.xyShifts = zeros(size(self.stageXY));
                end
            end
            % save display pars
            % 'eTableXls','pars','chnNames',
            savePars = struct2table(pars,'AsArray',true);
            writetable(savePars,[self.saveFolder,'ParsExploreMosaic.csv']);        
            [~,nChns] = size(self.imageTiles);
            ul = self.stageXY + self.xyShifts(pars.fovs,:); % stage positions + fiducial drift alignment;
            mosaicIm = cell(nChns,1);
            disp('converting tiles to mosaic...');
            tic
            for c=1:nChns
                imTiles = self.imageTiles(:,c);
%                 if ~strcmp(pars.backgroundCorrect,'none')  % LoadRegDaxAsMosaic already does this!
%                     imTiles = FlattenBackground(imTiles,'parameters',pars);
%                 end
                [mosaicIm{c},self.ulOut]= TilesToMosaic(imTiles,ul,'padMosaic',1);
            end
            toc
            mosaicIm = cat(3,mosaicIm{:});
            self.mosaicImage = mosaicIm;
            % figure(2); clf; imagesc(mosaicIm(:,:,end)); colormap(gray)
            % should have some flatten this dataset functions 
            nca = NcolorApp(mosaicIm,'names',self.chnNames);  
            
%            % troubleshoot tiles
%             nTiles = size(imTiles,1);
%             [h,w] = size(imTiles{1});
%             fovBoxes  = [ul, repmat([h,w],nTiles,1)];
%             fovBoxes(:,1:2) = fovBoxes(:,1:2) - repmat(min(ul,[],1),nTiles,1)+repmat([h,w],nTiles,1);
%             boxPostions = fovBoxes*(nca.scale)^nca.currScale;    
%             for f=1:nTiles % plot all the boxes 
%                 rectangle('Position',boxPostions(f,:),'EdgeColor','y');
%             end
%             
%             figure(1); clf; imagesc(IncreaseContrast(mosaicIm,'high',.999)); hold on;
%             for f=1:nTiles % plot all the boxes 
%                 rectangle('Position',fovBoxes(f,:),'EdgeColor','y');
%             end
            disp('mosaic rendered');
        end
        
        function self = SaveTif(self,pars)
            % default save root
            if isempty(pars.saveName)
                saveRoot = 'BuiltMosaic';
            else
                saveRoot = pars.saveName;
            end
            % default save folder
            if isempty(pars.saveFolder)
                outputFolder = self.saveFolder;
            else
                outputFolder = pars.saveFolder;
            end
            tifLabels = [outputFolder,saveRoot,'Labels.txt'];
            tifName = [outputFolder,saveRoot,'.tif'];
            % check if files exist
            [~,~,nChns] = size(self.mosaicImage);
            skip = false;
            if exist(tifLabels,'file')~=0 
                answer = questdlg([tifName ' exists, overwrite?'], ... % question
                    'Prompt ', ...  % pop-up label
                    'Yes','No','No'); % op1 op2 default
                switch answer
                    case 'Yes'
                        % delete(tifName);
                        skip = false;
                    case 'No'
                        skip = true;
                end         
            end
            if ~skip
                % write table of channel names               
                writetable(table(self.chnNames),tifLabels);
                % write table of uls
                ulTable = [outputFolder,saveRoot,'tileULsMB.csv'];
                uls = self.ulOut;% self.stageXY + self.xyShifts;
                x = uls(:,1);
                y = uls(:,2);
                [w,h] = size(self.imageTiles{1});
                w= w*ones(length(x),1);
                h= h*ones(length(x),1);
                trim = self.stepParameters{5}.trimBorder*ones(length(x),1);
%                 fliplr= self.stepParameters{5}.fliplr*ones(length(x),1);
%                 flipud= self.stepParameters{5}.flipud*ones(length(x),1);
%                 transpose= self.stepParameters{5}.transpose*ones(length(x),1);
                scope= repmat(self.stepParameters{5}.scope,length(x),1);
                fov = Column(self.stepParameters{5}.fovs);
                ulTableData = table(x,y,w,h,trim,scope,fov);
                writetable(ulTableData,ulTable);
                for c=1:nChns
                  	mosaicOut= self.mosaicImage(:,:,c); 
                    tifSaveName = regexprep(tifName,'.tif',['_c',num2str(c,'%03d'),'.tif']);
                    try
                        imwrite(mosaicOut,tifSaveName);
                    catch er
                        warning(er.getReport);
                        error(['error writing ',tifSaveName, '.  Verify you have sufficient disk space and that save directory exists.']);
                    end
                    % imwrite(mosaicOut,tifName,'WriteMode','append');
                end
                disp(['finished writing ',tifName]);                     
            end
        end
        
        function self  = BrowseMosaic(self,pars)
            NcolorApp(self.mosaicImage,'names',self.chnNames)
        end
        
    end
end