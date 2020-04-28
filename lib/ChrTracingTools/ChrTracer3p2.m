function varargout = ChrTracer3p2(varargin)
% CHRTRACER3P2 MATLAB code for ChrTracer3p2.fig
%      CHRTRACER3P2, by itself, creates a new CHRTRACER3P2 or raises the existing
%      singleton*.
%
%      H = CHRTRACER3P2 returns the handle to a new CHRTRACER3P2 or the handle to
%      the existing singleton*.
%
%      CHRTRACER3P2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHRTRACER3P2.M with the given input arguments.
%
%      CHRTRACER3P2('Property','Value',...) creates a new CHRTRACER3P2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ChrTracer3p2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ChrTracer3p2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
%
%% For improvement
% keep only brightest spot + spots within a fraction (say 25%) of brightest
%   done
% adaptive contrast for spot fitting, maybe with a fixed base threshold
%   done
% toss out 'non-circles' (based on 95% CI limits for XY, low a-value (a-h ratio or b-h ratio)). 
%   done
% starting spot based on clustering, instead of first frame?
%   done
% parallel processing for LinkSpots / fiducial drift correct
% mark selected link spots in fig2 and 4 (change color?)
% Parameter GUIs: more parameter types (drop-down lists).  
% auto-flag poor fits
%       - too many spots per hybe or lots of missing 'good' hybes
%       - 
% Troubleshoot load FOV data

% Edit the above text to modify the response to help ChrTracer3p2

% Last Modified by GUIDE v2.5 13-Apr-2019 09:22:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ChrTracer3p2_OpeningFcn, ...
                   'gui_OutputFcn',  @ChrTracer3p2_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ChrTracer3p2 is made visible.
function ChrTracer3p2_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ChrTracer3p2 (see VARARGIN)

% all data for the chromosome tracer is kept inside a structure array;
% so there can be multiple version of ChrTracer3p2 open, each version has its
% own version ID.  we create a cell array 
global CT
if isempty(CT)
    CT = cell(1,1);
else
    CT = [CT;cell(1,1)];
end
id = length(CT);
set(handles.CTinstance,'String',['inst id',num2str(id)]);
handles.id = id;
set(handles.ButtonRunStep,'String','Load Ex. Table');

% Choose default command line output for ChrTracer3p2
handles.output = hObject;

% setup other default parameters
DefaultParmeters(hObject,eventdata,handles);

% set current FOV save folder
% EditTextFOV_Callback(hObject, eventdata, handles)

% Update handles structure
guidata(hObject, handles);

error('This is a development branch, not for analysis. Please use ChrTracer3!');

% UIWAIT makes ChrTracer3p2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ChrTracer3p2_OutputFcn(hObject, eventdata, handles) 
global CT; %#ok<NUSED>
varargout{1} = handles.output;



function DefaultParmeters(hObject,eventdata,handles)
% These default parameters are passed forth to the indicated function
% They can also be edited in the popup GUI
% 
% 
global CT scratchPath;

% initialize some variables that get checked as empty or not;
CT{handles.id}.regData = {};
CT{handles.id}.fiducialAlignFrames = [];
CT{handles.id}.allLociXY = {};
CT{handles.id}.dataAllFrames = {};

% Load Experiment Table ------------------------------------------------% 
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true}; 
defaults(end+1,:) = {'dataFolder', 'string', ''}; 
defaults(end+1,:) = {'dataTypes', 'cell', {'any'}}; 
defaults(end+1,:) = {'daxRootDefault', 'string', 'ConvZscan*.dax'}; 
defaults(end+1,:) = {'overwrite', 'boolean', false}; 
defaults(end+1,:) = {'maxProject', 'boolean', true}; 
defaults(end+1,:) = {'saveProject', 'boolean', true}; 
defaults(end+1,:) = {'selectFOVs', 'freeType', inf}; 
CT{handles.id}.parsLoadExpTable = ParseVariableArguments([],defaults,'ChrTracer3_LoadExpTable');

% Fix Drift ------------------------------------------------% 
defaults = cell(0,3);
% TO DO:  add here other parameters from ChrTracer3_FixGlobal Drift
defaults(end+1,:) = {'selectFOVs','freeType',inf};
defaults(end+1,:) = {'saveFits', 'boolean', true}; 
% key alignment parameters
defaults(end+1,:) = {'alignContrastLow', 'fraction', .7}; % low image threshold for contrast balance prior to coarse alignment
defaults(end+1,:) = {'alignContrastHigh', 'fraction', .9995}; % high threshold  for contrast balance prior to coarse alignment
defaults(end+1,:) = {'refHybe','integer',1}; % hybe to use to start alignment
% other CorrAlignFast Parameters
defaults(end+1,:) = {'maxSize', 'positive', 400}; % rescale all images to this size for alignment
defaults(end+1,:) = {'fineBox', 'freeType', 100};  % perform fine scale alignment using a box of this size around the brightest point.
defaults(end+1,:) = {'fineUpsample', 'positive', 4};  
defaults(end+1,:) = {'maxShift', 'nonnegative', inf};
defaults(end+1,:) = {'gradMax', 'boolean', true};
defaults(end+1,:) = {'minGrad', 'float', -inf};
defaults(end+1,:) = {'angles','float',0}; % -10:1:10
defaults(end+1,:) = {'scales','float',1}; % -10:1:10
defaults(end+1,:) = {'fineMaxShift', 'nonnegative', 10};
defaults(end+1,:) = {'fineAngles','float',0}; % -1:.1:1
defaults(end+1,:) = {'fineScales','float',1}; % 0.95:0.01:.1.05
defaults(end+1,:) = {'fineCenter','array',[0,0]};
defaults(end+1,:) = {'showplot', 'boolean', true};
defaults(end+1,:) = {'fastDisplay', 'boolean', true};
defaults(end+1,:) = {'displayWidth', 'integer', 500};
defaults(end+1,:) = {'showExtraPlot', 'boolean', false};
defaults(end+1,:) = {'minFineImprovement', 'float', 0}; 
defaults(end+1,:) = {'showCorrAlign', 'boolean', true};
CT{handles.id}.parsFixDrift = ParseVariableArguments([],defaults,'ChrTracer3_FixGlobalDrift');


% FOV summary ----------------------------------------------%
defaults = cell(0,3);
defaults(end+1,:) = {'contrastLow','fraction',0};
defaults(end+1,:) = {'contrastHigh','fraction',.999};
defaults(end+1,:) = {'showFidXY','boolean',true};
defaults(end+1,:) = {'rescaleTile','positive',.1};
defaults(end+1,:) = {'colormap','colormap','hsv'};
defaults(end+1,:) = {'saveFigFile', 'boolean', true}; 
defaults(end+1,:) = {'selectFOVs','freeType',inf};
CT{handles.id}.parsFovOverlay  = ParseVariableArguments([],defaults,'ChrTracer3_FOVsummaryPlots');


% Select Spots ----------------------------------------------%
defaults = cell(0,3);
defaults(end+1,:) = {'fov','integer',1};
defaults(end+1,:) = {'autoSelectThreshold','fraction',.997};
defaults(end+1,:) = {'autoSelectDownsample','positive',3};
defaults(end+1,:) = {'border','nonnegative',2};
defaults(end+1,:) = {'coarseBackgroundSize','nonnegative',50}; % 50
defaults(end+1,:) = {'showExtraPlots','boolean',false};
defaults(end+1,:) = {'laplaceFilter','boolean',false};
defaults(end+1,:) = {'edgeThresh','nonnegative',0}; %.0001
defaults(end+1,:) = {'edgeSigma','nonnegative',1.2};
defaults(end+1,:) = {'removeEdgeStack','nonnegative',20};
defaults(end+1,:) = {'saveData','boolean',false};
defaults(end+1,:) = {'useDataChns','freeType',0};
defaults(end+1,:) = {'imContrastHigh','fraction',0.9995};
defaults(end+1,:) = {'imContrastLow','fraction',0.5};
CT{handles.id}.parsSelectSpots = ParseVariableArguments([],defaults,'ChrTracer_SelectSpots');


% Crop Spots ----------------------------------------------%
defaults = cell(0,3);
defaults(end+1,:) = {'selectSpot', 'positive', 1};
defaults(end+1,:) = {'boxWidth', 'positive', 16};
defaults(end+1,:) = {'showFolderNames', 'boolean',false};
defaults(end+1,:) = {'showPlots','boolean',true};
defaults(end+1,:) = {'fov','integer',1};
CT{handles.id}.parsCrop  = ParseVariableArguments([],defaults,'ChrTracer_CropAndPlot');

% Fit Spots  ----------------------------------------------%
defaults = cell(0,3);
% Common Data Fitting parameters
defaults(end+1,:) = {'keepBrightest','integer',1};  % Max number of peaks to allow
defaults(end+1,:) = {'maxXYstep','positive',5}; % box radius in pixels
defaults(end+1,:) = {'maxZstep','positive',7}; % box radius in pixels
defaults(end+1,:) = {'datMinSep', 'nonnegative', 3};
defaults(end+1,:) = {'datMinPeakHeight', 'positive', 500};
% Common Fid Fitting pars
defaults(end+1,:) = {'fidMinPeakHeight', 'positive', 200};
defaults(end+1,:) = {'maxXYdrift','positive',4};
defaults(end+1,:) = {'maxZdrift','positive',6};
defaults(end+1,:) = {'upsample','positive',4};  % 8 for accuracy 2 for speed
% Less common Fiducial fitting defaults
defaults(end+1,:) = {'fidMinSep', 'nonnegative', 5};  % Min separation between peaks in pixels.  Closer than this will be averaged
defaults(end+1,:) = {'fidKeepBrightest','integer',1};  % Max number of peaks to allow
defaults(end+1,:) = {'fidMaxFitWidth', 'positive', 8}; % 6
defaults(end+1,:) = {'fidMaxFitZdepth', 'positive', 12}; % 6
% Less common Data fitting defaults
defaults(end+1,:) = {'datMaxFitWidth', 'positive', 8};
defaults(end+1,:) = {'datMaxFitZdepth', 'positive', 12};
defaults(end+1,:) = {'datMinHBratio','nonnegative',1.2}; % peak value over background value
defaults(end+1,:) = {'datMinAHratio','nonnegative',.25}; % fitted height over background vs peak value
defaults(end+1,:) = {'datMaxUncert','nonnegative',2}; % pixels
% less commmon general pars
defaults(end+1,:) = {'nmXYpix','positive',154};
defaults(end+1,:) = {'nmZpix','positive',100};
defaults(end+1,:) = {'showPlots','boolean',true}; 
defaults(end+1,:) = {'saveData','boolean',false}; 
defaults(end+1,:) = {'overwrite','boolean',false}; % override FOV pars overwrite
CT{handles.id}.parsFit = ParseVariableArguments([],defaults,'ChrTracer_FitSpots');

% Select All Spots -----------------------
defaults = cell(0,3);
defaults(end+1,:) = {'selectFOVs','freeType',inf};
CT{handles.id}.parsSelectAllSpots = ParseVariableArguments([], defaults, 'AutoSelectSpots');


% Fit All Spots ------------------
defaults = cell(0,3);
defaults(end+1,:) = {'numParallel','integer',12};
defaults(end+1,:) = {'overwrite',{'Yes','No','Ask'},'No'};
CT{handles.id}.parsFitAll = ParseVariableArguments([],defaults,'ChrTracer_FitAll');

% default Steps
CT{handles.id}.stepNames = {'Load Ex. Table';
                            'Fix Global Drift';  % Save FOV Overlay
                            'Validate Drift Fix';
                            'Select Spots';
                            'Crop Spots';
                            'Fit Spots';
                            'Select All Spots';
                            'Fit All'};
CT{handles.id}.currStep = 'Load Ex. Table';

% step directions
CT{handles.id}.stepDirs = {...
'Select an experiment table to load and specify a directory in which to save the data.';
'Click "Fix Global Drift" to compute corrections x,y stage drift. Select "Step Pars" to specify a subset of FOVs if desired.';
'Click "Validate Drift Fix" to see results of x,y drift correction.  If image looks white instead of color shifted, proceed to "Next Step".';
'Click "Select Spots" to auto-ID spots.  Click "Step Pars" to change fitting parameters for selection.';
'Click "Step Pars" to select the number of a spot to plot, and adjust crop-box size (see figure from "Select Spots" step). Then click "Crop Spots".';
'Click "Step Pars" to adjust fit thresholds and fit parameters. Then click "Fit Spots".';
'Click "Step Pars" to specify FOV to process. Click "Select All Spots" to process first FOV. Click "Next Step" to record data for current FOV and advance to next.';
'Click "Fit All" to apply current fitting parameters to all FOV.'
};

% --- Executes on button press in ButtonRunStep.
function ButtonRunStep_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to ButtonRunStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CT

currStepName = CT{handles.id}.currStep;


%---------------------------- LOAD DATA ----------------------------------%    
if strcmp(currStepName,'Load Ex. Table')
        CT{handles.id}.expTableXLS = get(handles.EditTextFolder,'String'); 
        CT{handles.id}.saveFolder = get(handles.EditSaveFolder,'String');
        if isempty(CT{handles.id}.expTableXLS)
            error('Please enter an experiment table');
        end
        set(handles.TextDir,'String','Projecting data to allow global drift correction...'); guidata(hObject, handles);
        pars = CT{handles.id}.parsLoadExpTable;
        [CT{handles.id}.fiducialFrames,CT{handles.id}.eTable,...
         CT{handles.id}.rawDataNames] = ...
            ChrTracer3p2_LoadExpTable(CT{handles.id}.expTableXLS,'parameters',pars);
        % setup next step
        CT{handles.id}.parsFixDrift.selectFOVs = pars.selectFOVs; % carry forward 
        
        
%---------------------------- Fix Global Drift --------------------------%    
elseif strcmp(currStepName,'Fix Global Drift')
    pars = CT{handles.id}.parsFixDrift;
    if isempty(CT{handles.id}.fiducialFrames)
       error('no Fiducial Frames found, please rerun Load Ex. Table'); 
    end
    [CT{handles.id}.fiducialAlignFrames,CT{handles.id}.regData] = ...
        ChrTracer3p2_FixGlobalDrift(CT{handles.id}.fiducialFrames,'parameters',pars,'saveFolder',CT{handles.id}.saveFolder);
    % setup next step
    CT{handles.id}.currStep = 'Validate Drift Fix';
    
    
%---------------------------- Validate Drift Fix --------------------------%    
elseif strcmp(currStepName,'Validate Drift Fix')    
    validation = ChrTracer3p3_ValidateDriftFix(CT{handles.id}.saveFolder);
    if any(validation)
        CT{handles.id}.currStep = 'Fix Global Drift';
    else
        CT{handles.id}.currStep = 'Select Spots';
    end

%---------------------------- SELECT SPOTS ---------------------------------%
elseif strcmp(currStepName,'Select Spots') % SELECT ROI
    pars = CT{handles.id}.parsSelectSpots;
    pars.showPlots = true; % required 
    pars.numberSpots = true; % required
    pars.refHyb = CT{handles.id}.parsFixDrift.refHybe;
    pars.eTableXLS = CT{handles.id}.expTableXLS;    
    CT{handles.id}.currFOV = pars.fov;
    % ------- Manually select a spot--------------
    im = LoadDaxFromEtable(pars.eTableXLS,...
                    'dataType','fiducial',...
                    'fov',pars.fov,...
                    'hybNumber',pars.refHyb);
    CT{handles.id}.lociXY = ChrTracer3p2_ManualSelectSpot(im);
    %--------------------------------------------------
    % % Auotmated selection of spots
    % SelectSpots(handles,pars); % will populate CT{h}.lociXY
    % CT{handles.id}.allLociXY{pars.fov} = CT{handles.id}.lociXY;

    
%----------------------------- CROP SPOTS --------------------------------%    
elseif strcmp(currStepName,'Crop Spots') % PLOT SPOTS
     pars = CT{handles.id}.parsCrop;
     f = CT{handles.id}.currFOV; % pars.fov;
     pars.currentSpot = CT{handles.id}.parsCrop.selectSpot;
     pars.eTable = CT{handles.id}.eTable;
     pars.saveData = false; 
     % check that rawDataNames are loaded
     if isempty(CT{handles.id}.rawDataNames(:,f))
        error(['Did not find raw dax file names for FOV ',num2str(f), '. Please re-run Load Exp. Table.']) 
     end
     %  was CT{handles.id}.regData{f}, now moved to use csv
     regData = LoadRegData(CT{handles.id}.saveFolder,'fov',f); 
     % main function
     [CT{handles.id}.fidSpots,CT{handles.id}.dataSpots] = ...
         ChrTracer3p3_CropSpots( CT{handles.id}.rawDataNames(:,f),regData,...
              CT{handles.id}.eTable,CT{handles.id}.lociXY,'parameters',pars);
          
     CT{handles.id}.currStep = 'Fit Spots'; % Auto advance to next step
     
%-------------------------------- FIT SPOTS ------------------------------%  
elseif strcmp(currStepName,'Fit Spots') 
    pars = CT{handles.id}.parsFit;
    s = CT{handles.id}.parsCrop.selectSpot;
    pars.eTable = CT{handles.id}.eTable;  
    pars.lociXY = CT{handles.id}.lociXY;
    spotTableTemp = ChrTracer3_FitSpots(CT{handles.id}.fidSpots,CT{handles.id}.dataSpots,s,'parameters',pars,'veryverbose',true);
    CT{handles.id}.spotTableTemp = spotTableTemp;


%------------------------ SELECT ALL SPOTS -------------------------------%
elseif strcmp(currStepName,'Select All Spots') 
% Fit all spots
    ChrTracer_SpotSelector('id',handles.id,'fov',CT{handles.id}.currFOV)
   

%-------------------------- FIT ALL ---------------------------------%
elseif strcmp(currStepName,'Fit All') 
    pars = CT{handles.id}.parsFitAll;
    
    % remove fovs which lack ROI spots or are not requested.
    spotFiles = ls([CT{handles.id}.saveFolder,'fov*_selectSpots.csv']);
    fovsWithSpots = str2num(spotFiles(:,4:6))'; %#ok<ST2NM>
    
    for f=fovsWithSpots
        disp(['Processing data for FOV ',num2str(f)]);       
        tableSaveName = [CT{handles.id}.saveFolder, 'fov',num2str(f,'%03d'),'_AllFits.csv'];
        if exist(tableSaveName,'file')
            if strcmp(pars.overwrite,'Ask')
                answer = questdlg([tableSaveName ' exists, overwrite?'], ... % question
                    'Prompt ', ...  % pop-up label
                    'Yes','No','No'); % op1 op2 default
            else
                answer = pars.overwrite;
            end
        else
            answer = 'Yes';
        end
        % Handle response
        switch answer
            case 'Yes'
                % crop spots
                spotXY = LoadSpotXY(CT{handles.id}.saveFolder,'fov',f); % CT{handles.id}.allLociXY{f}; % Load from CSV             
                regData = LoadRegData(CT{handles.id}.saveFolder,'fov',f); 
                [CT{handles.id}.fidiSpots,CT{handles.id}.dataSpots] = ...
                ChrTracer3p3_CropAllSpots(CT{handles.id}.rawDataNames(:,f),...
                                      regData,...
                                      CT{handles.id}.eTable,...
                                      spotXY,...
                                    'parameters',CT{handles.id}.parsCrop,...
                                    'showPlots',false,...
                                    'numParallel',pars.numParallel);

                % fit spots
                spotDataTable = ChrTracer3_FitAllSpots(...
                    CT{handles.id}.fidiSpots,CT{handles.id}.dataSpots,...
                    'parameters',CT{handles.id}.parsFit,...
                    'eTable',CT{handles.id}.eTable,...
                    'numParallel',1,... pars.numParallel,... 
                    'lociXY',spotXY,...
                    'showPlots',false,...
                    'fov',f); % 

                % save data table
                writetable(spotDataTable,tableSaveName);
                disp(['wrote ',tableSaveName]);

            case 'No'
                disp(['skipping fov ',num2str(f)]);
        end
    end
end
% set(handles.ButtonRunStep,'String',CT{handles.id}.currStep);
SetupSteps(hObject,eventdata,handles);
guidata(hObject, handles);
%-------------------------------------------------------------------------% 


% % ---- Called by Select Spots
% %  currently replaced by manual selection
% function SelectSpots(handles,pars)
%     global CT
%         im = LoadDaxFromEtable(pars.eTableXLS,...
%                     'dataType','fiducial',...
%                     'fov',pars.fov,...
%                     'hybNumber',pars.refHyb);
%     fileName = ['fov',num2str(f,'%03d'),'_selectSpots.csv'];
%     saveName = [CT{handles.id}.saveFolder,fileName];
%     if exist(saveName,'file')
%         answer = questdlg('found a previous spot map in save folder.', ... % question
%                         'Skip', ...  % pop-up label
%                         'Load','Recompute','Load'); % op1 op2 default
%         if strcmp(answer,'Load')
%             spots = readtable(saveName);
%             spots = spots{:,:};
%             CT{handles.id}.lociXY = spots;
%             figure(1); clf; imagesc(im); hold on;
%             plot(spots(:,1),spots(:,2),'yo');
%             text(spots(:,1)+2,spots(:,2),cellstr(num2str( (1:size(spots,1))')),'color','w');
%         else
%             CT{handles.id}.lociXY = AutoSelectSpots(im,'parameters',pars,'showPlots',true,'numberSpots',true);
%         end
%     else
%         CT{handles.id}.lociXY = AutoSelectSpots(im,'parameters',pars,'showPlots',true,'numberSpots',true);      
%     end
%     figure(1); title(['fov ',num2str(f)]);


% --------------------------------------------------------
function GetNewPars(hObject, eventdata, handles)
% hObject    handle to ButtonParameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CT
currStepName = CT{handles.id}.currStep;
if strcmp(currStepName,'Load Ex. Table')
    CT{handles.id}.parsLoadExpTable = SimpleParameterGUI(CT{handles.id}.parsLoadExpTable);
elseif strcmp(currStepName,'Fix Global Drift')
    CT{handles.id}.parsFixDrift = SimpleParameterGUI(CT{handles.id}.parsFixDrift); 
elseif strcmp(currStepName,'Validate Drift Fix')
    CT{handles.id}.parsFovOverlay = SimpleParameterGUI(CT{handles.id}.parsFovOverlay);
elseif strcmp(currStepName,'Save Aligned Data')
    CT{handles.id}.parsSaveRegData = SimpleParameterGUI(CT{handles.id}.parsSaveRegData);
elseif strcmp(currStepName,'MemMap Data')
    CT{handles.id}.parsMemMapData = SimpleParameterGUI(CT{handles.id}.parsMemMapData);
elseif strcmp(currStepName,'Select Spots')
    CT{handles.id}.parsSelectSpots = SimpleParameterGUI(CT{handles.id}.parsSelectSpots); 
elseif strcmp(currStepName,'Crop Spots')
    CT{handles.id}.parsCrop = SimpleParameterGUI(CT{handles.id}.parsCrop); 
elseif strcmp(currStepName,'Fit Spots')
    CT{handles.id}.parsFit = SimpleParameterGUI(CT{handles.id}.parsFit);
elseif strcmp(currStepName,'Select All Spots')
    CT{handles.id}.parsSelectAllSpots = SimpleParameterGUI(CT{handles.id}.parsSelectAllSpots,'maxFieldsPerPrompt',15);
    CT{handles.id}.parsSelectSpots =  JoinStructures(CT{handles.id}.parsSelectAllSpots,CT{handles.id}.parsSelectSpots);
elseif strcmp(currStepName,'Fit All')
    CT{handles.id}.parsFitAll = SimpleParameterGUI(CT{handles.id}.parsFitAll);

else
    disp(['parameter GUI not yet written for step ', num2str(CT{handles.id}.currStep)]);
end

% --- Executes on button press in ButtonParameters.
function ButtonParameters_Callback(hObject, eventdata, handles)
GetNewPars(hObject, eventdata, handles)

% --------------------------------------------------------------------
function MenuStepPars_Callback(hObject, eventdata, handles)
% hObject    handle to MenuStepPars (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
GetNewPars(hObject, eventdata, handles)


% % Obsolete, moved to SpotSelector
% function SaveSpotTable(lociXY,f,handles)
%     global CT
%     fileName = ['fov',num2str(f,'%03d'),'_selectSpots.csv'];
%     saveName = [CT{handles.id}.saveFolder,fileName];
%     saveTable = table(lociXY(:,1),lociXY(:,2));
%     saveTable.Properties.VariableNames = {'locusX','locusY'};
%     answer = 'Yes';
%     if exist(saveName,'file')
%          answer = questdlg('Overwrite existing spot table?', ... % question
%                         'Overwrite?', ...  % pop-up label
%                         'Yes','No','Yes'); % op1 op2 default   
%     end
%     if strcmp(answer,'Yes')
%         writetable(saveTable,saveName);
%         disp(['wrote ',saveName]);
%     elseif strcmp(answer,'No') % load existing table into memory
%         lociXY = readtable(saveName);
%         CT{handles.id}.lociXY=[lociXY{:,1},lociXY{:,2}];
%         disp('loaded previously saved spot table into memory');
%     end

function MenuSaveFOVdata_Callback(hObject, eventdata, handles)
    global CT
    % SaveSpotTable(CT{handles.id}.lociXY,CT{handles.id}.currFOV,handles);

function MenuLoadFOVdata_Callback(hObject, eventdata, handles)
    global CT
    [file,folder] = uigetfile({'*.csv'},'Select a spot file',CT{handles.id}.saveFolder);
    spots = readtable([folder,file]);
    CT{handles.id}.lociXY = spots{:,:};
    disp(['loaded ',file]);

% --- Internal Function, Clear all spot data from current FOV
function ClearFOVdata(handles)
% might be useful, needs updating

% --- Executes on button press in ButtonBackStep.
function ButtonBackStep_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonBackStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global CT
    % Use the step buttons to move through FOVs
    changeStep = true;
%     if strcmp(CT{handles.id}.currStep,'Select All Spots')
%         changeStep = false;
%         % selectFOVs = CT{handles.id}.parsSelectAllSpots.selectFOVs;
%         if CT{handles.id}.currFOV > 1
%            CT{handles.id}.currFOV = CT{handles.id}.currFOV - 1;
%            txtDir = ['Click "Select All Spots" to ID spots for FOV: ',num2str(CT{handles.id}.currFOV)];
%            set(handles.TextDir,'String',txtDir);
%         else
%            changeStep = true; 
%         end
%     end
    % change step
    if changeStep
        step = find(strcmp( CT{handles.id}.stepNames, CT{handles.id}.currStep )); 
        if step > 1
            step = step - 1;
        else
            warning('Already reset to step 1'); 
        end
        CT{handles.id}.currStep =  CT{handles.id}.stepNames{step};
        SetupSteps(hObject,eventdata,handles);
    end
    guidata(hObject, handles);


% --- Executes on button press in ButtonNextStep.
function ButtonNextStep_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonNextStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global CT
    changeStep = true;
%     if strcmp(CT{handles.id}.currStep,'Select All Spots')
%         changeStep = false;
%         selectFOVs = CT{handles.id}.parsSelectAllSpots.selectFOVs;
%         if any(isinf(selectFOVs))
%            selectFOVs = 1:size(CT{handles.id}.fiducialFrames,2); 
%         end
%         if CT{handles.id}.currFOV < max(selectFOVs)         
%            SaveSpotTable(CT{handles.id}.lociXY,CT{handles.id}.currFOV,handles);
%            CT{handles.id}.allLociXY{CT{handles.id}.currFOV} = CT{handles.id}.lociXY;
%            CT{handles.id}.currFOV = CT{handles.id}.currFOV + 1;
%            txtDir = ['Spots recorded for FOV ',num2str(CT{handles.id}.currFOV-1), '.  Click "Select All Spots" to ID spots for FOV: ',num2str(CT{handles.id}.currFOV)];
%            set(handles.TextDir,'String',txtDir);
%         elseif CT{handles.id}.currFOV == max(selectFOVs) 
%             SaveSpotTable(CT{handles.id}.lociXY,CT{handles.id}.currFOV,handles);
%             CT{handles.id}.allLociXY{CT{handles.id}.currFOV} = CT{handles.id}.lociXY;
%            
%            txtDir = ['Spots recorded for final FOV ',num2str(CT{handles.id}.currFOV), '.  Click "Next Step" to advance to "Fit All"'];
%            set(handles.TextDir,'String',txtDir);
%            changeStep = true; 
%         else
%            changeStep = true; 
%         end
%     end
    % change step
    if changeStep
        % auto-save step parameters
        currStepName = CT{handles.id}.currStep;
        saveName = [CT{handles.id}.saveFolder,'Pars ',currStepName,'.csv'];
        SaveParsAsCSV(hObject, eventdata, handles,saveName);
        
        % advance to next step
        step = find(strcmp( CT{handles.id}.stepNames, currStepName));
        numSteps = length(CT{handles.id}.stepNames);
        if step < numSteps
            step = step + 1;
        else
             warning('Already reached last step'); 
        end
        CT{handles.id}.currStep =  CT{handles.id}.stepNames{step};
        SetupSteps(hObject,eventdata,handles);     
    end
    guidata(hObject, handles);

function SetupSteps(hObject,eventdata,handles)
    global CT
    currStepName = CT{handles.id}.currStep;
    step = find(strcmp( CT{handles.id}.stepNames, CT{handles.id}.currStep ));
    set(handles.ButtonRunStep,'String',CT{handles.id}.currStep);
    set(handles.TextDir,'String',CT{handles.id}.stepDirs{step}); %#ok<FNDSB>
    if strcmp(currStepName,'Select Spots') || strcmp(currStepName,'Select All Spots') 
        set(handles.ButtonEditManually1,'Visible','Off');
        set(handles.ButtonEditManually1,'String','Remove Spots');
        set(handles.ButtonEditManually2,'Visible','Off');
        set(handles.ButtonEditManually2,'String','Add Spots');
    elseif strcmp(currStepName,'Fix Global Drift')
        set(handles.ButtonEditManually1,'Visible','Off'); % not active, switched off
        set(handles.ButtonEditManually1,'String','Color Correct');
    else
        set(handles.ButtonEditManually1,'Visible','Off');
        set(handles.ButtonEditManually2,'Visible','Off');
    end
    guidata(hObject, handles);




% --- Executes on button press in ButtonEditManually1.
function ButtonEditManually1_Callback(hObject, eventdata, handles)
% This button is only visible for steps which allow manual edits
% The function of the button can be different depending on which step of
% the analysis the GUI is currently on. 
    % 
    global CT 
    currButtonName = get(handles.ButtonEditManually1,'String');
    if strcmp(currButtonName,'Color Correct')
        disp('function not available yet');
    elseif strcmp(currButtonName,'Remove Spots')
        lociXY = CT{handles.id}.lociXY;
        fig = figure(1);
        selectDir = 'Use the mouse to select spots to add (left click). Right click when done. Tip: to zoom in, do so before selecting "Remove Spots"';
        set(handles.TextDir,'String',selectDir);
        spotsToRemove = SelectSpotsInFig(fig,'verbose',false);

        % remove existing data
        lineHandles = findobj(fig,'Type','Line');
        for j=1:length(lineHandles)
           delete(lineHandles(j)); 
        end
        textHandles = findobj(fig,'Type','Text');
        for j=1:length(textHandles)
           delete(textHandles(j)); 
        end

        % update plot and lociXY list
        lociXY(spotsToRemove,:) = [];
        figure(fig); 
        plot(lociXY(:,1),lociXY(:,2),'yo'); hold on;
        text(lociXY(:,1)+2,lociXY(:,2),cellstr(num2str( (1:size(lociXY,1))')),'color','w');
        CT{handles.id}.lociXY = lociXY;
    end



% --- Executes on button press in ButtonEditManually2.
function ButtonEditManually2_Callback(hObject, eventdata, handles)
% This button is only visible for steps which allow manual edits
% The function of the button can be different depending on which step of
% the analysis the GUI is currently on. 
    global CT 
    currButtonName = get(handles.ButtonEditManually2,'String');
    if  strcmp(currButtonName,'Add Spots')
        lociXY = CT{handles.id}.lociXY;
        fig = figure(1);
        selectDir = 'Use the mouse to select spots to add (left click). Right click when done. Tip: to zoom in, do so before selecting "Add Spots"';
        set(handles.TextDir,'String',selectDir);
        addedSpots = AddSpotsToFig(fig,'verbose',false);

        % remove existing data
        lineHandles = findobj(fig,'Type','Line');
        for j=1:length(lineHandles)
           delete(lineHandles(j)); 
        end
        textHandles = findobj(fig,'Type','Text');
        for j=1:length(textHandles)
           delete(textHandles(j)); 
        end
        % update plot and lociXY list
        lociXY = cat(1,lociXY,round(addedSpots));
        figure(fig); 
        plot(lociXY(:,1),lociXY(:,2),'yo'); hold on;
        text(lociXY(:,1)+2,lociXY(:,2),cellstr(num2str( (1:size(lociXY,1))')),'color','w');
        CT{handles.id}.lociXY = lociXY;
    end


% --------------------------------------------------------------------
function ParametersMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ParametersMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MenuSavePars_Callback(hObject, eventdata, handles)
global CT
[file,path] = uiputfile('*.csv','Save parameters as',CT{handles.id}.saveFolder);
saveName = [path,file];
SaveParsAsCSV(hObject, eventdata, handles,saveName)

function SaveParsAsCSV(hObject, eventdata, handles,saveName)
global CT
try
    currStepName = CT{handles.id}.currStep;
    if strcmp(currStepName,'Load Ex. Table')
        parsOut = struct2table(CT{handles.id}.parsLoadExpTable);
    elseif strcmp(currStepName,'Fix Global Drift')
        parsOut = struct2table(CT{handles.id}.parsFixDrift ); 
    elseif strcmp(currStepName,'Validate Drift Fix')
        parsOut = struct2table(CT{handles.id}.parsFovOverlay );
    elseif strcmp(currStepName,'Save Aligned Data')
        parsOut = struct2table(CT{handles.id}.parsSaveRegData );
    elseif strcmp(currStepName,'MemMap Data')
        parsOut = struct2table(CT{handles.id}.parsMemMapData);
    elseif strcmp(currStepName,'Select Spots')
        parsOut = struct2table(CT{handles.id}.parsSelectSpots ); 
    elseif strcmp(currStepName,'Crop Spots')
        parsOut = struct2table(CT{handles.id}.parsCrop); 
    elseif strcmp(currStepName,'Fit Spots')
        parsOut = struct2table(CT{handles.id}.parsFit);
    elseif strcmp(currStepName,'Select All Spots')
        parsOut = struct2table(CT{handles.id}.parsSelectAllSpots);
    elseif strcmp(currStepName,'Fit All')
        parsOut = struct2table(CT{handles.id}.parsFitAll);
    end
    writetable(parsOut,saveName);
catch 
    % some parameters have cell arrays which can't be converted to tables
    % later we should do something about that. 
end

% --------------------------------------------------------------------
function MenuLoadPars_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function EditSaveFolder_Callback(hObject, eventdata, handles)
% hObject    handle to EditSaveFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CT
CT{handles.id}.saveFolder = get(handles.EditSaveFolder,'String');
SetFigureSavePath(CT{handles.id}.saveFolder,'makeDir',true);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function EditSaveFolder_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --------------------------------------------------------------------
function AnalysisMenu_Callback(hObject, eventdata, handles)
% hObject    handle to AnalysisMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% ============== OBSOLETE Menu commands (to remove)
% need to remove from GUI using GUIDE:
function MenuAnalyzeAllFOV_Callback(hObject, eventdata, handles)
function MenuAnalyzeAll_Callback(hObject, eventdata, handles)
function MenuRegisterAll_Callback(hObject, eventdata, handles)
function MenuFitAllAnotated_Callback(hObject, eventdata, handles)
function MenuLoadExpTable_Callback(hObject, eventdata, handles)
%================================================================

% --------------------------------------------------------------------
function MenuFOVpars_Callback(hObject, eventdata, handles)



%================== OBSOLETE other internal functions (to remove)
% not embedded in GUI, should be able to just delete
function AutoRegister(hObject,eventdata,handles,f)
function AutoFitData(hObject,eventdata,handles,f)
function RemoveLocus(hObject,eventdata,handles)
%=============================================================

% ------------------------------------------------------------------
% NOT sure if this is still active:
function EditTextFolder_Callback(hObject, eventdata, handles) %#ok<*INUSD>


% --- Executes during object creation, after setting all properties.
function EditTextFolder_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% OBSOLETE:
% function EditTextFOV_Callback(hObject, eventdata, handles)
% global CT
% CT{handles.id}.parsFOV.fov = str2double(get(handles.EditTextFOV,'String'));
% disp(['current FOV = ',num2str(CT{handles.id}.parsFOV.fov)]);
% 
% 
% % --- Executes during object creation, after setting all properties.
% function EditTextFOV_CreateFcn(hObject, eventdata, handles)
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end


% Obsolete
% --------------------------------------------------------------------
% function eTable = ReadExperimentTable(hObject,eventdata,handles)
% global CT
% % load experiment table
% dataTypes = CT{handles.id}.parsLoadData.dataTypes;
% CT{handles.id}.expTableXLS = get(handles.EditTextFolder,'String'); 
% eTable = readtable(CT{handles.id}.expTableXLS);
% hybFolderID = false(height(eTable),1);
% if sum(  strcmp(dataTypes,'any') ) == 0
%     for i=1:length(dataTypes)
%         hybFolderID = hybFolderID | strcmp(eTable.DataType,dataTypes{i});
%     end
% else
%     hybFolderID = true(height(eTable),1);
% end
% eTable = eTable(hybFolderID,:);
% CT{handles.id}.eTable = eTable;    
    


% 
% % Obsolete, moved to Fid frames
% function imFrame = GetFidFrame(handles,f)
% global CT
% saveMaxProject = true;
% saveFolder = CT{handles.id}.saveFolder;
% refHyb = CT{handles.id}.parsFixDrift.refHybe;
% useDataChns = CT{handles.id}.parsSelectSpots.useDataChns;
% useFid = ~any(useDataChns);
% if useFid
%     imFrame = CT{handles.id}.fiducialFrames{refHyb,f};
% else
%     % --Use mean projection of data channel instead of fiducial data-----
%    if ~isempty(CT{handles.id}.dataAllFrames(:,:,f))
%         stackFrames = cat(3,CT{handles.id}.dataAllFrames{:,:,f});
%         imFrame = mean(stackFrames,3); % average seems to work better than max, we benefit from some overlap 
%    else
%        % if data isn't already loaded in RAM, let's try load it from the hard disk
%        maxFiles = ls([saveFolder,'max_fov',num2str(f,'%03d'),'_*data*.dax']);
%        if ~isempty(maxFiles)
%            maxFiles = cellstr(maxFiles);
%             datFlatFrames = cell(length(maxFiles),1);
%            for d= useDataChns % 1:length(maxFiles)
%             datFlatFrames{d} = ReadDax([saveFolder,maxFiles{d}],'verbose',false);
%            end
%            stackFrames = cat(3,datFlatFrames{:});
%            imFrame = mean(stackFrames,3); 
%        end
%        if isempty(maxFiles) 
%            % if there aren't any max projections of the data files, lets make them now 
%          daxFiles  = ls([saveFolder,'fov',num2str(f,'%03d'),'_*data*.dax']);
%          if ~isempty(daxFiles)
%              daxFiles = cellstr(daxFiles);
%          end
%          datFlatFrames = cell(length(daxFiles),1);
%          for d=useDataChns
%             [dax,infoFile] = ReadDax([saveFolder,daxFiles{d}],'verbose',false);
%             datFlatFrames{d} = max(dax,[],3);
%             % --- Save these max projections, we may want them later ----
%             if saveMaxProject
%                 infoOut = infoFile;
%                 infoOut.number_of_frames = 1;
%                 infoOut.localName = regexprep(['max_',daxFiles{d}],'.dax','.inf');
%                 WriteDAXFiles(datFlatFrames{d},infoOut,'verbose',true);
%             end
%          end
%         stackFrames = cat(3,datFlatFrames{:});
%         imFrame = mean(stackFrames,3); 
%        end
%    end
% end

