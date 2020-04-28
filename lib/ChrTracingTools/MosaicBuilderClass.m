classdef MosaicBuilderClass < handle

% MosaicBuilder specific dynamic properties
    properties % editable parameters
        gui_h;
        currStep = 'Load Data'
        saveFolder = '';
        daxInHybFov = {};
        maxNameInHybFov = {};
        imageTiles = {};
        stageXY = [];
        mosaicGray = [];
        xyShifts =[];
    end
    
    % MosaciBuilder specific static properties
    %   though all stepGUIs will need to define these properties
    properties % (Access = private) % 
       allStepNames = {'Load Data',...
                       'Align Hybs',...
                       'Validate Mosaic',...
                       'Align FOV'}
       stepDirections = {['Step 1: Create max projections of all FOVs. Optionally, also project all Hybs (all barcodes).',...
                        'Projecting all Hybs is generally unnecessary for ORCA DNA data files, just set hybs=1.'];
                        'Step 2: Correct hyb-to-hyb drift across all requested hybs.';
                        'Step 3: Validate Mosaic';
                        'Step 4: Correct fov-to-fov drift across all FOV.'}
       allFxnHandles = {@(inputs,pars) LoadData(inputs,pars),...
                        @(inputs,pars) AlignHybs(inputs,pars),...
                        @(inputs,pars) ValidateMosaic(inputs,pars),...
                        @(inputs,pars) AlignFOV(inputs,pars)};
       stepParameters;
    end

    
    methods            
        function self = MosaicBuilderClass
            % make GUI handle to the fig.m file and store it locally
            self.gui_h = guihandles(MultistepAnalysisGUI);
        
            % set the callback function for the "Run Step" button.
            set(self.gui_h.ButtonRun,'callback',@(src,event) ButtonRun_callback(self,src,event));
            
            % set the callback function for the "Next Step" button.
            set(self.gui_h.ButtonNext,'callback',@(src,event) ButtonNext_callback(self,src,event));
            
            % set the callback function for the "Back Step" button.
            set(self.gui_h.ButtonBack,'callback',@(src,event) ButtonBack_callback(self,src,event));
            
            % set the callback function for the "Parameters" button.
            set(self.gui_h.ButtonPars,'callback',@(src,event) ButtonPars_callback(self,src,event));
            
            % setup default parameters
            % Load Data
            defaults = cell(0,3);
            defaults(end+1,:) = {'chnTableFile','string','S:\2018-05-15_WT_2-4hrs(HoxRNA)_(3kbDNA)\RNA_Expt\rnaTableFile.xlsx'};
            defaults(end+1,:) = {'dataType',{'fiducial','data','all'},'all'};
            defaults(end+1,:) = {'hybNumber', 'freeType', inf}; 
            defaults(end+1,:) = {'fov','freeType',inf};
            parsLoadData  = ParseVariableArguments([],defaults,'MosaicBuilder LoadData');
            % Align Hybs
            defaults = cell(0,3);
            defaults(end+1,:) = {'refHyb', 'integer', 1}; 
            parsAlignHybs  = ParseVariableArguments([],defaults,'MosaicBuilder AlignHybs');
            % Validate Mosaic
            defaults = cell(0,3);
            defaults(end+1,:) = {'flatten','boolean',true};
            defaults(end+1,:) = {'method',{'mean','sum','edgeBlur','last','first'},'edgeBlur'}; 
            defaults(end+1,:) = {'mosaicContrastLow','fraction',.1};
            defaults(end+1,:) = {'mosaicContrastHigh','fraction',.9999};
            defaults(end+1,:) = {'fovs','freeType',inf};
            defaults(end+1,:) = {'hybs','freeType',inf};
            defaults(end+1,:) = {'dataType',{'fiducial','data','all'},'data'};
            parsValidateMosaic  = ParseVariableArguments([],defaults,'MosaicBuilder ValidateMosaic');
            % AlignFOV
            defaults = cell(0,3);
            defaults(end+1,:) = {'interactive','boolean',true};
            defaults(end+1,:) = {'flatten','boolean',true};
            defaults(end+1,:) = {'method',{'mean','sum','edgeBlur','last','first'},'edgeBlur'}; 
            defaults(end+1,:) = {'mosaicContrastLow','fraction',.1};
            defaults(end+1,:) = {'mosaicContrastHigh','fraction',.9999};
            defaults(end+1,:) = {'fovs','freeType',inf};
            parsAlignFOV  = ParseVariableArguments([],defaults,'MosaicBuilder AlignFOV');
            % stack
            self.stepParameters = {parsLoadData,...
                             parsAlignHybs,...
                             parsValidateMosaic,...
                             parsAlignFOV};

            %sets the figure close function. This lets the class know that
            %the figure wants to close and thus the class should cleanup in
            %memory as well
            set(self.gui_h.figure1,  'closerequestfcn', @(src,event) Close_fcn(self, src, event));
          
        end
        
    end
    % ==== These methods could be part of a MultistepAnalysisClass =======
    methods 
        % ---- Executes when "Run Step" button is clicked
        function ButtonRun_callback(self,src,event) %#ok<*INUSD>
            % calls the current step behavior using the current step
            % parameters
            %    see protected functions below
            for s=1:length(self.allStepNames)
                if strcmp(self.currStep,self.allStepNames(s))
                    pars = self.stepParameters{s};
                    self.allFxnHandles{s}(self,pars);
                end
            end
        end
        
        % ---- Executes when "Parameters" button is clicked
        function ButtonPars_callback(self,src,event)
            % Opens a GUI to edit the current step parameters
            for s=1:length(self.allStepNames)
                if strcmp(self.currStep,self.allStepNames(s))
                    self.stepParameters{s} = SimpleParameterGUI(self.stepParameters{s});
                end
            end
        end
        
        % ---- Executes when "Next Step" button is clicked
        function ButtonNext_callback(self,src,event)           
            step = find(strcmp(self.allStepNames, self.currStep));
            numSteps = length(self.allStepNames);
            if step < numSteps
                 step = step + 1;
                 self.currStep = self.allStepNames{step};
                 set(self.gui_h.TextDir,'String',self.stepDirections{step});
            else
                warning('Already reached last step'); 
            end
        end
         
        % ---- Executes when "Back Step" button is clicked
        function ButtonBack_callback(self,src,event)
            step = find(strcmp(self.allStepNames, self.currStep));
            if step > 1
                 step = step - 1;
                 self.currStep = self.allStepNames{step};
                 set(self.gui_h.TextDir,'String',self.stepDirections{step});
            else
                warning('Currently on step 1, can not go back.'); 
            end
        end
%     end
%     
%     
%     methods % (Access = private)
        function self = LoadData(self,pars) 
            saveName = fileparts(pars.chnTableFile);
            self.saveFolder = [saveName,filesep];
            pars.hybNumber = inf;
            % pars.fov = 1:5;
            pars.dataType = 'all';
            pars.veryverbose = false;
           [self.daxInHybFov,~,~,...
            self.maxNameInHybFov] =  LoadDaxFromEtable(pars.chnTableFile,...
                                        'parameters',pars);
%         figure(1); clf; 
%         h = 1; f=1; d=1;
%         imagesc(IncreaseContrast(self.daxInHybFov{h,f,d},'high',.999));
        end
        % ---- Align Hybs 
        function self = AlignHybs(self,pars)
            self = ChrTracer3p2_FixGlobalDrift(self.daxInHybFov(:,:,1),...
                'saveFolder',self.saveFolder,'parameters',pars);
        end
        % ---- Validate Mosaic
        function self = ValidateMosaic(self,pars)
            [numHybs,numFovs,numDatas] = size(self.maxNameInHybFov);
            if isinf(pars.hybs)
                h=1:numHybs;
            else
                h = pars.hybs;
            end
            if isinf(pars.fovs)
                f = 1:numFovs;
            else
                f = pars.fovs;
            end
            if strcmp(pars.dataType,'all')
                d = 1:numDatas;
                imTiles = self.maxNameInHyb(h,f,d);
                imTiles = reshape(permute(imTiles,[1,3,2]),numHybs*(numDatas),numFovs);   
            elseif strcmp(pars.dataType,'data')
                d = 2:numDatas;
                imTiles = self.maxNameInHyb(h,f,d);
                imTiles = reshape(permute(imTiles,[1,3,2]),numHybs*(numDatas-1),numFovs);
            elseif strcmp(pars.dataType,'fiducial')
                d = 1;
                imTiles = self.maxNameInHyb(h,f,d);
            end
            [self.imageTiles,self.stageXY] = DaxToImageTiles(imTiles,...
            'fliplr',false,'flipud',true,'loadInfoOnly',false,...
            'parameters',pars);
            MosaicViewerGUI(self.imageTiles,self.stageXY,'flatten',true,...
                'parameters',pars); % ,'method','mean'
        end

        % AlignFOV step
        function self = AlignFOV(self,pars)
            % enforce saveFolder filesep end;
            if ~isempty(self.saveFolder)
                if ~strcmp(self.saveFolder(end),filesep)
                    self.saveFolder = [self.saveFolder,filesep]; 
                end
            end
            % parse fov request
            if isinf(pars.fovs)
                fovs = 1:length(self.imageTiles);
            else
                fovs=pars.fovs;
            end
            pars.interactive = true; 
            pars.corrAlign = true;
            [self.mosaicGray,~,self.xyShifts] = ... consider supresssing all these outputs
                        MosaicViewerRender(self.imageTiles(fovs),...
                                           self.stageXY(fovs,:),...
                                           'parameters',pars,...
                                           'saveFolder',self.saveFolder);
           
                                       % these two should agree if all is well:                            
           mosaicCrop = CropEmptyPixels(mosaicGray);
           dispMosaic = IncreaseContrast(mosaicCrop,'high',.9999,'low',.5);
           figure(1); clf; imagesc(dispMosaic); colormap(gray); colorbar;

           % MosaicViewerGUI(self.imageTiles(fovs),self.stageXY(fovs,:)+xyShifts);
        end
        
        % Close Figure
        %  This (intentionally?) overloads the delete function
        %class deconstructor - handles the cleaning up of the class &
        %figure. Either the class or the figure can initiate the closing
        %condition, this function makes sure both are cleaned up
        function delete(self)
            %remove the closerequestfcn from the figure, this prevents an
            %infitie loop with the following delete command
            set(self.gui_h.figure1,  'closerequestfcn', '');
            %delete the figure
            delete(self.gui_h.figure1);
            %clear out the pointer to the figure - prevents memory leaks
            self.gui_h = [];
        end
        %function - Close_fcn
        %
        %this is the closerequestfcn of the figure. All it does here is
        %call the class delete function (presented above)
        function self = Close_fcn(self, src, event)
            disp('exiting MosaicBuilderClass');
            delete(self);
        end
        

    end
end