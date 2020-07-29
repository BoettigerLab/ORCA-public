function varargout = ChrTracer2(varargin)
% CHRTRACER2 MATLAB code for ChrTracer2.fig
%      CHRTRACER2, by itself, creates a new CHRTRACER2 or raises the existing
%      singleton*.
%
%      H = CHRTRACER2 returns the handle to a new CHRTRACER2 or the handle to
%      the existing singleton*.
%
%      CHRTRACER2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHRTRACER2.M with the given input arguments.
%
%      CHRTRACER2('Property','Value',...) creates a new CHRTRACER2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ChrTracer2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ChrTracer2_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help ChrTracer2

% Last Modified by GUIDE v2.5 09-Jul-2018 12:29:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ChrTracer2_OpeningFcn, ...
                   'gui_OutputFcn',  @ChrTracer2_OutputFcn, ...
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


% --- Executes just before ChrTracer2 is made visible.
function ChrTracer2_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ChrTracer2 (see VARARGIN)

% all data for the chromosome tracer is kept inside a structure array;
% so there can be multiple version of ChrTracer2 open, each version has its
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

% defaults
% (keep these)
CT{handles.id}.currentSpotID = 1; 
CT{handles.id}.spotTable = table();
CT{handles.id}.im = [];
CT{handles.id}.rejectPoints = [];

set(handles.ButtonRunStep,'String','Load Ex. Table');
% Choose default command line output for ChrTracer2
handles.output = hObject;

DefaultParmeters(hObject,eventdata,handles);

% set current FOV save folder
EditTextFOV_Callback(hObject, eventdata, handles)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ChrTracer2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ChrTracer2_OutputFcn(hObject, eventdata, handles) 
global CT;
varargout{1} = handles.output;



function DefaultParmeters(hObject,eventdata,handles)
% These default parameters are passed forth to the indicated function
% They can also be edited in the popup GUI
% 
% 
global CT scratchPath;
% FOV parameters ----------------------------------------------%
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true}; 
defaults(end+1,:) = {'veryverbose','boolean',false};
defaults(end+1,:) = {'showPlots', 'boolean', true}; 
defaults(end+1,:) = {'saveData', 'boolean', true}; 
defaults(end+1,:) = {'fov', 'integer', 1}; 
defaults(end+1,:) = {'lociXY','freeType',[]}; % for plotting previously selected spots 
defaults(end+1,:) = {'stopOnError','boolean',false};
defaults(end+1,:) = {'eTable','array',[]};
defaults(end+1,:) = {'numParallel', 'integer',1};
defaults(end+1,:) = {'saveFolder', 'string', scratchPath};
defaults(end+1,:) = {'dataFolder', 'string', ''};
defaults(end+1,:) = {'overwrite', 'boolean', false};
defaults(end+1,:) = {'separateFOVfolders', 'boolean', true};
CT{handles.id}.parsFOV  = ParseVariableArguments([],defaults,'ChrTracer');
% Load Data and Register ----------------------------------------------%
defaults = cell(0,3);
defaults(end+1,:) = {'filename', 'string', 'ConvZscan'}; 
defaults(end+1,:) = {'dataTypes','cell',{'any'}}; % hybe toe repeat expression
defaults(end+1,:) = {'alignmentBoxWidth', 'positive', 800}; % pixels.  size of region to use for frame alignment
defaults(end+1,:) = {'alignUpsample', 'positive', 1}; % upsample data for coarse alignment (slow, should be unnecessary);
defaults(end+1,:) = {'alignContrastLow', 'fraction', .8}; % low image threshold for contrast balance prior to coarse alignment
defaults(end+1,:) = {'alignContrastHigh', 'fraction', .999}; % high threshold  for contrast balance prior to coarse alignment
defaults(end+1,:) = {'refHybe','integer',1}; % hybe to use to start alignment
defaults(end+1,:) = {'previousAlignFrames','array',[]}; % previous align frames
defaults(end+1,:) = {'rotation','boolean',true}; % also compute and correct rotation angle
defaults(end+1,:) = {'maxD','positive',15}; % distance in pixels to match objects prior to computing rotation
defaults(end+1,:) = {'targetHybes','integer',[]};
defaults(end+1,:) = {'goodHybes','boolean',[]};
defaults(end+1,:) = {'corrAngles','freeType',[]}; % rotate angles to check for better alignment
defaults(end+1,:) = {'saveMaxProject','boolean',true}; % save maximum projection for quick loading in future
CT{handles.id}.parsLoadData = ParseVariableArguments([],defaults,'ChrTracer_LoadData');
% FOV summary ----------------------------------------------%
defaults = cell(0,3);
defaults(end+1,:) = {'contrastLow','fraction',0};
defaults(end+1,:) = {'contrastHigh','fraction',.999};
defaults(end+1,:) = {'gain','nonnegative',1};
defaults(end+1,:) = {'showFidXY','boolean',true};
defaults(end+1,:) = {'rescaleTile','positive',.1};
defaults(end+1,:) = {'colormap','colormap','hsv'};
defaults(end+1,:) = {'saveFigFile', 'boolean', true}; 
CT{handles.id}.parsFovOverlay  = ParseVariableArguments([],defaults,'ChrTracer_FOVsummaryPlots');
% Select ROI ----------------------------------------------%
defaults = cell(0,3);
defaults(end+1,:) = {'ROIselectMode',{'manual','auto'},'manual'};
defaults(end+1,:) = {'autoSelectThreshold','fraction',.992};
defaults(end+1,:) = {'autoSelectDownsample','positive',3};
defaults(end+1,:) = {'mergeFrames','string','median'};
defaults(end+1,:) = {'border','nonnegative',2};
defaults(end+1,:) = {'coarseBackgroundScale','nonnegative',0}; % 50
defaults(end+1,:) = {'fineBackgroundScale','nonnegative',0};   % 5
defaults(end+1,:) = {'gain','positive',1};
defaults(end+1,:) = {'useFid','boolean',true};
defaults(end+1,:) = {'displayContrastLow','fraction',0};
defaults(end+1,:) = {'displayContrastHigh','fraction',.999};
CT{handles.id}.parsROI = ParseVariableArguments([],defaults,'ChrTracer_ROIselect');
% Crop Spots ----------------------------------------------%
defaults = cell(0,3);
defaults(end+1,:) = {'selectSpot', 'positive', 1};
defaults(end+1,:) = {'boxWidth', 'positive', 16};
defaults(end+1,:) = {'goodHybes','array',[]};
defaults(end+1,:) = {'showFolderNames', 'boolean',false};
defaults(end+1,:) = {'fidGain','positive',1};
defaults(end+1,:) = {'datGain','positive',1};
CT{handles.id}.parsCrop  = ParseVariableArguments([],defaults,'ChrTracer_CropAndPlot');
% Fit Spots  ----------------------------------------------%
defaults = cell(0,3);
% Fiducial alignment defaults
defaults(end+1,:) = {'saveData','boolean',false}; 
defaults(end+1,:) = {'upsample','positive',4};  % 8 for accuracy 2 for speed
defaults(end+1,:) = {'maxShiftXY','positive',4};
defaults(end+1,:) = {'maxShiftZ','positive',6};
% Fiducial fitting defaults
defaults(end+1,:) = {'fidMinPeakHeight', 'positive', 200};
defaults(end+1,:) = {'fidMinSep', 'nonnegative', 8};  % Min separation between peaks in pixels.  Closer than this will be averaged
defaults(end+1,:) = {'fidKeepBrightest','integer',1};  % Max number of peaks to allow
% Data fitting defaults
defaults(end+1,:) = {'datBoxXY','positive',5}; % box radius in pixels
defaults(end+1,:) = {'datBoxZ','positive',7}; % box radius in pixels
defaults(end+1,:) = {'datMinPeakHeight', 'positive', 500};
defaults(end+1,:) = {'datMaxFitWidth', 'positive', 6};
defaults(end+1,:) = {'datMinSep', 'nonnegative', 3};
defaults(end+1,:) = {'datMinHBratio','nonnegative',1.2}; % peak value over background value
defaults(end+1,:) = {'datMinAHratio','nonnegative',.25}; % fitted height over background vs peak value
defaults(end+1,:) = {'datMaxUncert','nonnegative',2}; % pixels
defaults(end+1,:) = {'overwrite','boolean',true}; % override FOV pars overwrite
CT{handles.id}.parsFit = ParseVariableArguments([],defaults,'ChrTracer_FitSpots');
%-----------------Chromatic Aberation correction
defaults = cell(0,3);
defaults(end+1,:) = {'chromaticMapMat','string',''};
defaults(end+1,:) = {'applyToChannel','integer',0};
CT{handles.id}.parsColorCorrect = ParseVariableArguments([],defaults,'Color Correct');


% default Steps
CT{handles.id}.stepNames = {'Load Ex. Table';
                            'Validate Overlay';  % Save FOV Overlay
                            'Save Aligned Data';
                            'Select ROI';
                            'Crop Spots';
                            'Fit Spots'};
CT{handles.id}.currStep = 'Load Ex. Table';

% step directions
CT{handles.id}.stepDirs = {...
'Select an experiment table to load and specify a directory in which to save the data.';
'Select "Validate FOV Overlay" to see results of x,y drift correction.  If image looks white instead of color shifted, proceed to "Next Step".';
'Select "Save Aligned Data" to apply previously computed drift correction to all data files and save new aligned data files in the specified Save Folder.';
'Select "Step Pars" to select manual or auto ROI identification and to adjust automatic parameters. Then click "Select ROI".'
'Select "Step Pars" to adjust crop-box size or change spot number (see figure from ROI selection). Then click "Crop Spots".';
'Select "Step Pars" to adjust fit thresholds and fit parameters. Then click "Fit Spots".';
};

% --- Executes on button press in ButtonRunStep.
function ButtonRunStep_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to ButtonRunStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CT

currStepName = CT{handles.id}.currStep;

% % determine whether to use fov subfolders or not when saving data. 
% if CT{handles.id}.parsFOV.separateFOVfolders
%     f = str2double(get(handles.EditTextFOV,'String'));
%     fovSaveFolder = [CT{handles.id}.parsFOV.saveFolder,'fov',num2str(f,'%03d'),filesep];
%     if exist(fovSaveFolder,'dir')~=7
%         mkdir(fovSaveFolder);
%     end
%     CT{handles.id}.fovSaveFolder = fovSaveFolder;
% else
%     fovSaveFolder = CT{handles.id}.parsFOV.saveFolder;
% end

%---------------------------- LOAD DATA ----------------------------------%    
if strcmp(currStepName,'Load Ex. Table')
    if ~isempty(CT{handles.id}.parsFOV.lociXY)
        cprintf('red','Loading a new FOV will clear data from current FOV. Are you sure?')
        cprintf('red','Tip: you can save current data from the file menu');
        choice = input('Continue 1=Yes, 0=No. ?  ');
    else
        choice = true;
    end
    if choice 
        ClearFOVdata(handles); % cleanup
        CT{handles.id}.expTableXLS = get(handles.EditTextFolder,'String'); 
        CT{handles.id}.parsFOV.fov = str2double(get(handles.EditTextFOV,'String'));
        if isempty(CT{handles.id}.expTableXLS)
            error('Please enter an experiment table');
        end
        % pars = CT{handles.id};
        pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsLoadData);
        % now we actually load the table
        set(handles.TextDir,'String','Loading data and correcting x,y drift...'); guidata(hObject, handles);
        dataOut = ChrTracer2_LoadData(CT{handles.id}.expTableXLS,'parameters',pars);
        % store the data in the CT
        CT{handles.id}.fiducialFrames = dataOut.fiducialFrames;
        CT{handles.id}.regData = dataOut.regData;
        % update some of the parameters
        CT{handles.id}.parsFOV.goodHybes = dataOut.goodHybes;
        CT{handles.id}.parsFOV.saveFolder = dataOut.saveFolder;
        CT{handles.id}.parsFOV.dataFolder = dataOut.dataFolder;
        CT{handles.id}.parsFOV.eTable = dataOut.eTable;
        CT{handles.id}.numDataChns = dataOut.numDataChns;
        set(handles.EditSaveFolder,'String',CT{handles.id}.parsFOV.saveFolder);
        disp('Experiment Data Loaded'); 
         set(handles.TextDir,'String','Data loaded. Is x,y, drift correction okay?'); guidata(hObject, handles);
        CT{handles.id}.currStep = 'Validate Overlay';
        % Auto advance to next step
    else
        disp('Load table aborted'); 
    end
%---------------------------- SAVE FOV OVERLAYS --------------------------%    
elseif strcmp(currStepName,'Validate Overlay')
    pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsFovOverlay);
    ChrTracer2_FOVsummaryPlots(CT{handles.id}.fiducialFrames,'parameters',pars);
    CT{handles.id}.currStep = 'Save Aligned Data';
    disp('Select Save Aligned Data if the alignment looks good');
    disp('Otherwise select Back Step, change Parameters, and Reload Data');

%----------------------------- SAVE ALIGNED DATA --------------------------%   
elseif strcmp(currStepName,'Save Aligned Data') % 'Adjust Alignment'
    set(handles.TextDir,'String','Aligning data and building memory maps, please be patient...'); guidata(hObject, handles);
    pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsLoadData);
    [CT{handles.id}.fidMapData,CT{handles.id}.datMapData,...
     CT{handles.id}.fidFrames,CT{handles.id}.dataFrames] = ...
        ChrTracer2_SaveAlignedData(CT{handles.id}.parsFOV.eTable,CT{handles.id}.regData,'parameters',pars);
    CT{handles.id}.currStep = 'Select ROI';
%---------------------------- SELECT ROI ---------------------------------%
elseif strcmp(currStepName,'Select ROI') % SELECT ROI
    if CT{handles.id}.parsROI.useFid
        im = CT{handles.id}.fidFrames;
    else
        % a more robust selection using data if fiducial image is mostly
        % background. Note high background fiducial images can still make
        % good drift correction, but are poor for selecting spots. 
        nFrames = size(CT{handles.id}.dataFrames,3);
        im = cat(3,max(CT{handles.id}.dataFrames(:,:,1:floor(nFrames/3)),[],3),...
            max(CT{handles.id}.dataFrames(:,:,floor(nFrames/3)+1:2*floor(nFrames/3)),[],3),...
            max(CT{handles.id}.dataFrames(:,:,2*floor(nFrames/3)+1:end),[],3) );
    end  
    if isempty(im)
       set(handles.TextDir,'String','Error: fidMapData is empty. Run "Save Aligned Data" first and try again.'); guidata(hObject, handles);
       error('fidMapData is empty. Run "Save Aligned Data" first and try again');  
    end   
    pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsROI);
    [spot,im] = ChrTracer_ROIselect(im,'parameters',pars,'badSpots',CT{handles.id}.rejectPoints);
    CT{handles.id}.im = im;
    if strcmp(pars.ROIselectMode,'manual')
        CT{handles.id}.parsFOV.lociXY(CT{handles.id}.currentSpotID,:) = spot;
        CT{handles.id}.parsFOV.selectSpots = CT{handles.id}.currentSpotID;
    else
        numSpots =  size(spot,1);
        CT{handles.id}.parsFOV.lociXY = spot;
        CT{handles.id}.parsFOV.numSpots = numSpots;
    end
%----------------------------- CROP SPOTS --------------------------------%    
elseif strcmp(currStepName,'Crop Spots') % PLOT SPOTS
     pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsCrop);
     pars.currentSpot = CT{handles.id}.parsCrop.selectSpot;
     pars.eTable = CT{handles.id}.parsFOV.eTable;
     pars.saveData = false; 
     [CT{handles.id}.fidSpots,CT{handles.id}.dataSpots] = ...
         ChrTracer2_CropSpots(CT{handles.id}.fidMapData,CT{handles.id}.datMapData,'parameters',pars);
     CT{handles.id}.currStep = 'Fit Spots';
     SetupSteps(hObject,eventdata,handles);
     % Auto advance to next step
%-------------------------------- FIT SPOTS ------------------------------%  
elseif strcmp(currStepName,'Fit Spots') 
    pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsFit,'conflict','keepSecond');
    % pars.selectSpots = CT{handles.id}.currentSpotID;
    pars.selectSpots = CT{handles.id}.parsCrop.selectSpot;
    pars.showPlots = true;
    pars.saveData = false; 
    pars.numParallel = 1;
    spotTableTemp = ChrTracer2_FitAllSpots(CT{handles.id}.fidSpots,CT{handles.id}.dataSpots,'parameters',pars,'saveTable',false);
    CT{handles.id}.spotTableTemp = spotTableTemp;
    % Requires 'Next Step' to proceed
%--------------------------- COMPUTE MATRICES ----------------------------% 
elseif strcmp(currStepName,'Compute Matrices')
    disp('function not available yet'); 
%------------------------------- ACCEPT FITS -----------------------------% 
elseif strcmp(currStepName,'Accept Fits')
  disp('function not available yet');
end

% set(handles.ButtonRunStep,'String',CT{handles.id}.currStep);
SetupSteps(hObject,eventdata,handles);
guidata(hObject, handles);
%-------------------------------------------------------------------------% 




%======================= Batch Processing ================================%
function AnalyzeAllInFOV(hObject,eventdata,handles,f)
% Analyze all spots in FOV
global CT
CT{handles.id}.parsFOV.fov = f;
fovSaveFolder = [CT{handles.id}.parsFOV.saveFolder,'fov',num2str(f,'%03d'),filesep];
if exist(fovSaveFolder,'dir')~=7
    mkdir(fovSaveFolder);
end
% -------------------------- CLEAR UP -----------------------------------% 
ClearFOVdata(handles);

%---------------------------- LOAD DATA --------------------------%  
pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsLoadData);
saveFolder = CT{handles.id}.parsFOV.saveFolder;
savedExperimentTable = [saveFolder,filesep,'ExperimentLayout.xlsx'];
%if ~exist(savedExperimentTable,'file')~=0
%   system(['xcopy ',CT{handles.id}.expTableXLS,' ',savedExperimentTable]); 
%end
dataOut = ChrTracer2_LoadData(CT{handles.id}.expTableXLS,'parameters',pars);
% store the data in the CT
CT{handles.id}.fiducialFrames = dataOut.fiducialFrames;
CT{handles.id}.regData = dataOut.regData;
% update some of the parameters
CT{handles.id}.parsFOV.goodHybes = dataOut.goodHybes;
CT{handles.id}.parsFOV.eTable = dataOut.eTable;
CT{handles.id}.parsFOV.dataFolder = dataOut.dataFolder;
set(handles.EditSaveFolder,'String',saveFolder);
disp('Experiment Data Loaded'); 

%---------------------------- SAVE FOV OVERLAYS --------------------------%    
pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsFovOverlay);
CT{handles.id}.im = ChrTracer2_FOVsummaryPlots(CT{handles.id}.fiducialFrames,'parameters',pars);

%---------------------------- SAVE ALIGNED DATA---------------------------%
 pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsLoadData);
[CT{handles.id}.fidMapData,CT{handles.id}.datMapData,...
 CT{handles.id}.fidFrames,CT{handles.id}.dataFrames] = ...
 ChrTracer2_SaveAlignedData(pars.eTable,CT{handles.id}.regData,...
    'parameters',pars,...
    'loadFlatFid',CT{handles.id}.parsROI.useFid,...
    'loadFlatData',~CT{handles.id}.parsROI.useFid);
%---------------------------- SELECT ROI ---------------------------------%
if CT{handles.id}.parsROI.useFid
    im = CT{handles.id}.fidFrames;
else
    % a more robust selection using data if fiducial image is mostly
    % background. Note high background fiducial images can still make
    % good drift correction, but are poor for selecting spots. 
    nFrames = size(CT{handles.id}.dataFrames,3);
    im = cat(3,max(CT{handles.id}.dataFrames(:,:,1:floor(nFrames/3)),[],3),...
        max(CT{handles.id}.dataFrames(:,:,floor(nFrames/3)+1:2*floor(nFrames/3)),[],3),...
        max(CT{handles.id}.dataFrames(:,:,2*floor(nFrames/3)+1:end),[],3) );
end
pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsROI);
spot = ChrTracer_ROIselect(im,'parameters',pars,'badSpots',CT{handles.id}.rejectPoints);
numSpots =  size(spot,1);
CT{handles.id}.parsFOV.lociXY = spot;
CT{handles.id}.parsFOV.numSpots = numSpots;

%----------------------------- CROP SPOTS --------------------------------%   
 pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsCrop);
 pars.saveFolder = fovSaveFolder;
 pars.selectSpots = 1:size(CT{handles.id}.parsFOV.lociXY,1);
 pars.showPlots = false; % don't show plots when running in batch
[CT{handles.id}.fidSpots,CT{handles.id}.dataSpots] = ChrTracer2_CropAllSpots(CT{handles.id}.fidMapData,CT{handles.id}.datMapData,'parameters',pars);

%-------------------------------- FIT SPOTS ------------------------------%  
pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsFit,'conflict','keepFirst');
pars.saveFolder = fovSaveFolder;
pars.saveData = CT{handles.id}.parsFit.saveData;
pars.showPlots = false; % don't show plots when running in batch
pars.selectSpots = 1:numSpots;
spotDataTable = ChrTracer2_FitAllSpots(CT{handles.id}.fidSpots,CT{handles.id}.dataSpots,'parameters',pars,'saveTable',false);
tableSaveName = [CT{handles.id}.parsFOV.saveFolder, 'fov',num2str(pars.fov,'%03d'),'_AllFits.csv'];
writetable(spotDataTable,tableSaveName);
disp(['wrote ',tableSaveName]);

%-------------------------------- SAVE DATA ------------------------------% 
MenuSaveFOVdata_Callback(hObject, eventdata, handles);

%-------------------------------- CLEAN UP ------------------------------% 
% close down parpool to free up memory;
poolobj = gcp('nocreate');
delete(poolobj);






% ============== Menu commands
% --------------------------------------------------------------------
function MenuAnalyzeAllFOV_Callback(hObject, eventdata, handles)
%% Analyze all spots in all FOV
% hObject    handle to MenuAnalyzeAllFOV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CT
% set up some parallel processing
CT{handles.id}.parsFOV.numParallel = 20;
CT{handles.id}.parsROI.ROIselectMode = 'auto';
saveFolder = CT{handles.id}.parsFOV.saveFolder;
 
% determine number of FOV in data folder  
daxFiles1 = cellstr(ls([saveFolder,'*h001_fid.dax']));   
sourceFolder = fileparts(get(handles.EditTextFolder,'String'));
daxFiles2 = cellstr(ls([sourceFolder, filesep,CT{handles.id}.parsFOV.eTable.FolderName{1},filesep,'C*.dax']));
if length(daxFiles1) >= length(daxFiles2) %  this is to preseve backwards compatability with a now outdated form. Will be removed in a later release.
    daxFiles = daxFiles1;
else
    daxFiles = daxFiles2;
end
daxNums = 1:length(daxFiles);

% check for existing files
foundSpotTables = cellstr(ls([saveFolder,'fov*AllFits.csv']));
if ~isempty(foundSpotTables{1})
    fovComplete = unique(cellfun(@(x) str2double(x(4:6)), foundSpotTables)); 
else
    fovComplete = [];
end

% skip FOVs which already have ctFiles
% (later we can modify this with parameters for overwrite)
fin = intersect(daxNums,fovComplete);
if ~isempty(fin)
    disp(['found data for FOVs ',num2str(fin')]);
    disp('skipping these FOVs');
end
fovs = setdiff(daxNums,fovComplete);

for f=fovs
    AnalyzeAllInFOV(hObject,eventdata,handles,f);
end

% --------------------------------------------------------------
function MenuAnalyzeAll_Callback(hObject, eventdata, handles)
%% Analyze all spots in current FOV
% hObject    handle to MenuAnalyzeAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CT
CT{handles.id}.expTableXLS = get(handles.EditTextFolder,'String'); 
AnalyzeAllInFOV(hObject,eventdata,handles,CT{handles.id}.parsFOV.fov);

    




% --------------------------------------------------------------------
function MenuRegisterAll_Callback(hObject, eventdata, handles)
%% Register all data 
% this allows the reader to click through final max-projected images and
% fit spots one field of view at a time, before launching the batch fitting
% analysis. 
% hObject    handle to MenuRegisterAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CT
% determine number of FOV in data folder
sourceFolder = fileparts(get(handles.EditTextFolder,'String'));
rawDax = cellstr(ls([sourceFolder, filesep,CT{handles.id}.parsFOV.eTable.FolderName{1},filesep,'C*.dax']));
daxNums = 1:length(rawDax);
saveFolder = CT{handles.id}.parsFOV.saveFolder;
alignedData = cellstr(ls([saveFolder,'fov*h001_fid.dax']));
fovComplete = cellfun(@(x) str2double(x(4:6)),alignedData);
fovs = setdiff(daxNums,fovComplete);
for f=fovs
    AutoRegister(hObject,eventdata,handles,f);
end



% --------------------------------------------------------------------
function MenuFitAllAnotated_Callback(hObject, eventdata, handles)
%% Fit all data
% hObject    handle to MenuFitAllAnotated (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% check for existing files
global CT
saveFolder = CT{handles.id}.parsFOV.saveFolder;
CT{handles.id}.parsFOV.numParallel = 20;
alignedData = cellstr(ls([saveFolder,'fov*h001_fid.dax']));
if ~isempty(alignedData)
    fovAligned = cellfun(@(x) str2double(x(4:6)),alignedData);
else
    fovAligned = [];
end
foundSpotTables = cellstr(ls([saveFolder,'fov*AllFits.csv']));
if ~isempty(foundSpotTables{1})
    fovComplete = unique(cellfun(@(x) str2double(x(4:6)), foundSpotTables)); 
else
    fovComplete = [];
end
% skip FOVs which already have ctFiles
% (later we can modify this with parameters for overwrite)
fin = intersect(fovAligned,fovComplete);
if ~isempty(fin)
    disp(['found data for FOVs ',num2str(fin')]);
    disp('skipping these FOVs');
end
fovs = setdiff(daxNums,fovComplete);
for f=fovs
    CT{handles.id}.parsFOV.lociXY = CT{handles.id}.fovDatas{fov}.lociXY;% load files
    AutoFitData(hObject,eventdata,handles,f);
end






%================== Auto Register =====================================
function AutoRegister(hObject,eventdata,handles,f)
%% Register all data
% Run commands to load and register data
%           Register all
global CT
CT{handles.id}.parsFOV.fov = f;
fovSaveFolder = [CT{handles.id}.parsFOV.saveFolder,'fov',num2str(f,'%03d'),filesep];
if exist(fovSaveFolder,'dir')~=7
    mkdir(fovSaveFolder);
end
% -------------------------- CLEAR UP -----------------------------------% 
ClearFOVdata(handles);

%---------------------------- LOAD DATA --------------------------%  
pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsLoadData);
saveFolder = CT{handles.id}.parsFOV.saveFolder;
savedExperimentTable = [saveFolder,filesep,'ExperimentLayout.xlsx'];
if ~exist(savedExperimentTable,'file')~=0
   system(['xcopy ',CT{handles.id}.expTableXLS,' ',savedExperimentTable]); 
end
dataOut = ChrTracer2_LoadData(CT{handles.id}.expTableXLS,'parameters',pars);
% store the data in the CT
CT{handles.id}.fiducialFrames = dataOut.fiducialFrames;
CT{handles.id}.regData = dataOut.regData;
% update some of the parameters
CT{handles.id}.parsFOV.goodHybes = dataOut.goodHybes;
CT{handles.id}.parsFOV.eTable = dataOut.eTable;
CT{handles.id}.parsFOV.dataFolder = dataOut.dataFolder;
set(handles.EditSaveFolder,'String',saveFolder);
disp('Experiment Data Loaded'); 

%---------------------------- SAVE FOV OVERLAYS --------------------------%    
pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsFovOverlay);
CT{handles.id}.im = ChrTracer2_FOVsummaryPlots(CT{handles.id}.fiducialFrames,'parameters',pars);

%---------------------------- SAVE ALIGNED DATA---------------------------%
 pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsLoadData);
[CT{handles.id}.fidMapData,CT{handles.id}.datMapData,...
 CT{handles.id}.fidFrames,CT{handles.id}.dataFrames] = ...
 ChrTracer2_SaveAlignedData(pars.eTable,CT{handles.id}.regData,...
    'parameters',pars,...
    'loadFlatFid',CT{handles.id}.parsROI.useFid,...
    'loadFlatData',~CT{handles.id}.parsROI.useFid,...
    'memMapData',false,'verbose',false);
%=========================================================================%

%==================== Auto Fit ========================================
function AutoFitData(hObject,eventdata,handles,f)
%% Auto fit all data
% Run commands to auto fit all data
%           Register all
global CT
fovSaveFolder = [CT{handles.id}.parsFOV.saveFolder,'fov',num2str(f,'%03d'),filesep];
if exist(fovSaveFolder,'dir')~=7
    mkdir(fovSaveFolder);
end
%---------------------------- READ ALIGNED DATA---------------------------%
 pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsLoadData);
[CT{handles.id}.fidMapData,CT{handles.id}.datMapData,...
 CT{handles.id}.fidFrames,CT{handles.id}.dataFrames] = ...
 ChrTracer2_SaveAlignedData(pars.eTable,CT{handles.id}.regData,...
    'parameters',pars,...
    'loadFlatFid',CT{handles.id}.parsROI.useFid,...
    'loadFlatData',~CT{handles.id}.parsROI.useFid,...
    'memMapData',true,'verbose',false);
%----------------------------- CROP SPOTS --------------------------------%   
 pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsCrop);
 pars.saveFolder = fovSaveFolder;
 pars.selectSpots = 1:size(CT{handles.id}.parsFOV.lociXY,1);
 pars.showPlots = false; % don't show plots when running in batch
[CT{handles.id}.fidSpots,CT{handles.id}.dataSpots] = ChrTracer2_CropAllSpots(CT{handles.id}.fidMapData,CT{handles.id}.datMapData,'parameters',pars);
%-------------------------------- FIT SPOTS ------------------------------%  
pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsFit,'conflict','keepFirst');
pars.saveFolder = fovSaveFolder;
pars.saveData = CT{handles.id}.parsFit.saveData;
pars.showPlots = false; % don't show plots when running in batch
pars.selectSpots = 1:numSpots;
spotDataTable = ChrTracer2_FitAllSpots(CT{handles.id}.fidSpots,CT{handles.id}.dataSpots,'parameters',pars,'saveTable',false);
tableSaveName = [CT{handles.id}.parsFOV.saveFolder, 'fov',num2str(pars.fov,'%03d'),'_AllFits.csv'];
writetable(spotDataTable,tableSaveName);
disp(['wrote ',tableSaveName]);
%-------------------------------- SAVE DATA ------------------------------% 
MenuSaveFOVdata_Callback(hObject, eventdata, handles);
%-------------------------------- CLEAN UP ------------------------------% 
% close down parpool to free up memory;
poolobj = gcp('nocreate');
delete(poolobj);











% --- Executes on button press in ButtonParameters.
function ButtonParameters_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonParameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CT
currStepName = CT{handles.id}.currStep;
if strcmp(currStepName,'Load Ex. Table')
    CT{handles.id}.parsLoadData = ChrTracer_ParameterGUI(CT{handles.id}.parsLoadData);
elseif strcmp(currStepName,'Validate Overlay')
    CT{handles.id}.parsFovOverlay = ChrTracer_ParameterGUI(CT{handles.id}.parsFovOverlay); 
elseif strcmp(currStepName,'Save Aligned Data')
    CT{handles.id}.parsLoadData = ChrTracer_ParameterGUI(CT{handles.id}.parsLoadData);
elseif strcmp(currStepName,'Select ROI')
    CT{handles.id}.parsROI = ChrTracer_ParameterGUI(CT{handles.id}.parsROI); 
elseif strcmp(currStepName,'Crop Spots')
    CT{handles.id}.parsCrop = ChrTracer_ParameterGUI(CT{handles.id}.parsCrop); 
elseif strcmp(currStepName,'Fit Spots')
    CT{handles.id}.parsFit = ChrTracer_ParameterGUI(CT{handles.id}.parsFit);
elseif strcmp(currStepName,'Register and Plot')
    ChrTracer_LinkSpotsParsGUI(handles.id);
else
    disp(['parameter GUI not yet written for step ', num2str(CT{handles.id}.currStep)]);
end




% --- Internal Function, Clear all spot data from current FOV
function ClearFOVdata(handles)
global CT;
CT{handles.id}.parsFOV.lociXY = []; 
CT{handles.id}.currentSpotID = 1; % external, like step
CT{handles.id}.spotTable = table();
CT{handles.id}.im = [];
CT{handles.id}.rejectPoints = [];
CT{handles.id}.fiducialFrames = [];
CT{handles.id}.fidMapData = [];
CT{handles.id}.datMapData = [];
CT{handles.id}.fidSpots = [];
CT{handles.id}.dataSpots = [];
CT{handles.id}.regData = [];

% ---- Internal Function, remove last point
function RemoveLocus(hObject,eventdata,handles)
global CT
CT{handles.id}.parsFOV.lociXY(end,:) = [];
CT{handles.id}.stepName = 'Select ROI';
set(handles.ButtonRunStep,'String','Crop Spots');
guidata(hObject, handles);


% --- Executes on button press in ButtonBackStep.
function ButtonBackStep_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonBackStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CT
step = find(strcmp( CT{handles.id}.stepNames, CT{handles.id}.currStep )); 
if step > 1
    step = step - 1;
else
    warning('Already reset to step 1'); 
end
CT{handles.id}.currStep =  CT{handles.id}.stepNames{step};
SetupSteps(hObject,eventdata,handles);
guidata(hObject, handles);


% --- Executes on button press in ButtonNextStep.
function ButtonNextStep_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonNextStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CT
step = find(strcmp( CT{handles.id}.stepNames, CT{handles.id}.currStep ));
numSteps = length(CT{handles.id}.stepNames);
if step < numSteps
    step = step + 1;
else
     warning('Already reached last step'); 
end
CT{handles.id}.currStep =  CT{handles.id}.stepNames{step};
SetupSteps(hObject,eventdata,handles);
guidata(hObject, handles);

function SetupSteps(hObject,eventdata,handles)
global CT
currStepName = CT{handles.id}.currStep;
step = find(strcmp( CT{handles.id}.stepNames, CT{handles.id}.currStep ));
set(handles.ButtonRunStep,'String',CT{handles.id}.currStep);
set(handles.TextDir,'String',CT{handles.id}.stepDirs{step}); %#ok<FNDSB>
if strcmp(currStepName,'Crop Spots')
    set(handles.ButtonEditManually1,'Visible','On');
    set(handles.ButtonEditManually1,'String','Mark Bad Spot');
elseif strcmp(currStepName,'Select ROI')
    set(handles.ButtonEditManually1,'Visible','On');
    set(handles.ButtonEditManually1,'String','Next FOV');
    set(handles.ButtonEditManually2,'Visible','On');
    set(handles.ButtonEditManually2,'String','Record spots');
elseif strcmp(currStepName,'Validate Overlay');
    set(handles.ButtonEditManually1,'Visible','On');
    set(handles.ButtonEditManually1,'String','Color Correct');
else
    set(handles.ButtonEditManually1,'Visible','Off');
    set(handles.ButtonEditManually2,'Visible','Off');
end
guidata(hObject, handles);


% --------------------------------------------------------------------
function ParametersMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ParametersMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MenuFOVpars_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFOVpars (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CT
CT{handles.id}.parsFOV = ChrTracer_ParameterGUI(CT{handles.id}.parsFOV);

% --------------------------------------------------------------------
function MenuStepPars_Callback(hObject, eventdata, handles)
% hObject    handle to MenuStepPars (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CT
currStepName = CT{handles.id}.currStep;
if strcmp(currStepName,'Load Ex. Table')
    CT{handles.id}.parsLoadData = ChrTracer_ParameterGUI(CT{handles.id}.parsLoadData);
elseif strcmp(currStepName,'Validate Overlay')
    CT{handles.id}.parsFovOverlay = ChrTracer_ParameterGUI(CT{handles.id}.parsFovOverlay); 
elseif strcmp(currStepName,'Save Aligned Data')
    CT{handles.id}.parsLoadData = ChrTracer_ParameterGUI(CT{handles.id}.parsLoadData);
elseif strcmp(currStepName,'Select ROI')
    CT{handles.id}.parsROI = ChrTracer_ParameterGUI(CT{handles.id}.parsROI); 
elseif strcmp(currStepName,'Crop Spots')
    CT{handles.id}.parsCrop = ChrTracer_ParameterGUI(CT{handles.id}.parsCrop); 
elseif strcmp(currStepName,'Fit Spots')
    ChrTracer_PSFparsGUI(handles.id);
elseif strcmp(currStepName,'Register and Plot')
    CT{handles.id}.parsLink = ChrTracer_ParameterGUI(CT{handles.id}.parsLink); 
    % ChrTracer_LinkSpotsParsGUI(handles.id);
else
    disp(['parameter GUI not yet written for step ', num2str(CT{handles.id}.step)]);
end



% --------------------------------------------------------------------
function MenuLoadFOVdata_Callback(hObject, eventdata, handles)
% hObject    handle to MenuLoadFOVdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CT
if isempty(CT{handles.id}.parsFOV.saveFolder)
    startfolder = pwd;
else
    startfolder = CT{handles.id}.parsFOV.saveFolder;
end
disp('Current images and spot-tables will merge with loaded ones');
[fileName,pathName,filterIndex] = uigetfile({'*.mat','Matlab data (*.mat)';...
    '*.*','All Files (*.*)'},'Select FOV data',startfolder);
if filterIndex ~= 0 % loading operation was not canceled
    load([pathName,fileName]); %#ok<*LOAD>
    try
        disp('loaded CT data');
        disp(ctData.parsFOV.eTable);
        % merge spotTable, linkTable, loci, and cell array of zoom-in-spots.
        CT{handles.id}.spotTable = cat(1,CT{handles.id}.spotTable,ctData.spotTable);
        CT{handles.id}.cy5LinkTable = cat(1,CT{handles.id}.cy5LinkTable,ctData.cy5LinkTable);
        CT{handles.id}.parsFOV.lociXY = cat(1,CT{handles.id}.parsFOV.lociXY, ctData.parsFOV.lociXY);
        CT{handles.id}.fidSpots = cat(1,CT{handles.id}.fidSpots,ctData.fidSpots);
        CT{handles.id}.dataSpots = cat(1,CT{handles.id}.dataSpots,ctData.dataSpots);
        CT{handles.id}.currentSpotID = length(CT{handles.id}.currentSpotID) + 1;
    catch
        error(['error reading from ' pathName,fileName]); 
    end
end
    
% --------------------------------------------------------------------
function MenuSaveFOVdata_Callback(hObject, eventdata, handles)
% hObject    handle to MenuSaveFOVdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CT
ctData = CT{handles.id};
ctData.fiducialFrames = [];
ctData.fidSpots = [];
ctData.dataSpots = [];
saveFolder = ctData.parsFOV.saveFolder;
fov = ctData.parsFOV.fov;
disp(['saving data as ',saveFolder,'fov',num2str(fov,'%03d'),'_ctData.mat']);
save([saveFolder,'fov',num2str(fov,'%03d'),'_ctData.mat'],'ctData'); 


% --------------------------------------------------------------------
function MenuLoadExpTable_Callback(hObject, eventdata, handles)
% hObject    handle to MenuLoadExpTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CT
if isempty(CT{handles.id}.parsFOV.saveFolder)
    startfolder = pwd;
else
    startfolder = CT{handles.id}.parsFOV.saveFolder;
end
[fileName,pathName,filterIndex] = uigetfile({'*.xslx','Excel Table (*.xslx)';...
    '*.*','All Files (*.*)'},'Select Experiment Table',startfolder);
if filterIndex ~= 0 % loading operation was not canceled
    expTableXLS = [pathName,filesep,fileName]; %  'E:\Kay\2017-02-06_en\ExperimentLayout.xlsx';
    if CT{handles.id}.verbose
        disp(['loading ',expTableXLS]);
    end
    set(handles.EditTextFolder,'String',expTableXLS);
    CT{handles.id}.expTableXLS = expTableXLS;
    guidata(hObject,handles); 
    ReadExperimentTable(hObject,eventdata,handles)
end


% --------------------------------------------------------------------
function MenuSavePars_Callback(hObject, eventdata, handles)
% hObject    handle to MenuSavePars (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global CT
parsData.parsFOV = CT{handles.id}.parsFOV;
parsData.parsLoadData = CT{handles.id}.parsLoadData;
parsData.parsFovOverlay = CT{handles.id}.parsFovOverlay;
parsData.parsROI = CT{handles.id}.parsROI;
parsData.parsCrop = CT{handles.id}.parsCrop;
parsData.parsFit = CT{handles.id}.parsFit;
parsData.parsLink = CT{handles.id}.parsLink;
saveFolder = parsData.parsFOV.saveFolder;
fov = parsData.parsFOV.fov;
disp(['saving data as ',saveFolder,'fov',num2str(fov,'%03d'),'_parsData.mat']);
save([saveFolder,'fov',num2str(fov,'%03d'),'_parsData.mat'],'parsData'); 


% --------------------------------------------------------------------
function MenuLoadPars_Callback(hObject, eventdata, handles)
% hObject    handle to MenuLoadPars (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CT
if isempty(CT{handles.id}.parsFOV.saveFolder)
    startfolder = pwd;
else
    startfolder = CT{handles.id}.parsFOV.saveFolder;
end

[fileName,pathName,filterIndex] = uigetfile({'*parsData.mat','Matlab Parameters (*.mat)';...
    '*.*','All Files (*.*)'},'Select Saved Parameters',startfolder);
if filterIndex ~= 0 % loading operation was not canceled 
    load([pathName,fileName]);
    try
        disp(['loaded ', pathName,fileName]);
        CT{handles.id}.parsFOV = parsData.parsFOV;
        CT{handles.id}.parsLoadData = parsData.parsLoadData;
        CT{handles.id}.parsFovOverlay = parsData.parsFovOverlay;
        CT{handles.id}.parsROI = parsData.parsROI;
        CT{handles.id}.parsCrop = parsData.parsCrop;
        CT{handles.id}.parsFit = parsData.parsFit;
        CT{handles.id}.parsLink = parsData.parsLink; 
    catch
        error(['error reading from ' pathName,fileName]); 
    end
end


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function EditTextFolder_Callback(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to EditTextFolder (see GCBO)
% eventdata  reserved - to be les.defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditTextFolder as text
%        str2double(get(hObject,'String')) returns contents of EditTextFolder as a double


% --- Executes during object creation, after setting all properties.
function EditTextFolder_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
function AnalysisMenu_Callback(hObject, eventdata, handles)
% hObject    handle to AnalysisMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function EditTextFOV_Callback(hObject, eventdata, handles)
% hObject    handle to EditTextFOV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CT
CT{handles.id}.parsFOV.fov = str2double(get(handles.EditTextFOV,'String'));
disp(['current FOV = ',num2str(CT{handles.id}.parsFOV.fov)]);
% CT{handles.id}.currStep = 'Load Ex. Table';
% set(handles.ButtonRunStep,'String',CT{handles.id}.currStep);
% CT{handles.id}.parsLoadData.targetHybes = []; % clear target hybes
% guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function EditTextFOV_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in ButtonEditManually1.
function ButtonEditManually1_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonEditManually1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 
% This button is only visible for steps which allow manual edits
% The function of the button can be different depending on which step of
% the analysis the GUI is currently on. 
% 
global CT 
currButtonName = get(handles.ButtonEditManually1,'String');
if strcmp(currButtonName,'Color Correct')
    CT{handles.id}.parsColorCorrect = SimpleParameterGUI(CT{handles.id}.parsColorCorrect);
    if fid~=0
        disp(['using ',CT{handles.id}.parsColorCorrect.chromaticMapMat,...
            ' for chromatic correction to channel ',...
            num2str(CT{handles.id}.parsColorCorrect.applyToChannel)]);
    end
elseif strcmp(currButtonName,'Next FOV')
   CT{handles.id}.parsFOV.lociXY = [];
   if CT{handles.id}.parsROI.useFid
       fov = str2double(get(handles.EditTextFOV,'String'))+1;
       set(handles.EditTextFOV,'String',num2str(fov)); guidata(hObject,handles);
       saveFolder = CT{handles.id}.parsFOV.saveFolder;
       maxFiles = cellstr(ls([saveFolder,'max_fov',num2str(fov,'%03d'),'*','fid.dax']));
       if ~isempty(maxFiles)
           daxFiles = strcat(saveFolder,maxFiles);
           daxFrames = cell(length(daxFiles),1);
           for h=1 % :length(daxFiles)
               daxFrames{h} = ReadDax(daxFiles{h},'verbose',false);
           end
       else
           daxFiles = cellstr(ls([saveFolder,'fov',num2str(fov,'%03d'),'*','fid.dax']));
           if ~isempty(daxFiles)
               daxFrames = cell(length(daxFiles),1);
               for h=1 % :length(daxFiles)
                   daxFrames{h} = max(ReadDax(daxFiles{h}),[],3,'verbose',false);
               end
           else
               disp(['fov',num2str(fov,'%03d'),'*','fid.dax',' not found.']);
           end
       end
       daxFrames = cat(3,daxFrames{:});
       CT{handles.id}.fiducialFrames = daxFrames;
   else
       fov = str2double(get(handles.EditTextFOV,'String'))+1;
       set(handles.EditTextFOV,'String',num2str(fov)); guidata(hObject,handles);
       saveFolder = CT{handles.id}.parsFOV.saveFolder;
       maxFiles = cellstr(ls([saveFolder,'max_fov',num2str(fov,'%03d'),'*','data*.dax']));
       if ~isempty(maxFiles)
           daxFiles = strcat(saveFolder,maxFiles);
           daxFrames = cell(length(daxFiles),1);
           for h=1:length(daxFiles)
               daxFrames{h} = ReadDax(daxFiles{h},'verbose',false);
           end
       else
           daxFiles = cellstr(ls([saveFolder,'fov',num2str(fov,'%03d'),'*','fid.dax']));
           if ~isempty(daxFiles)
               daxFrames = cell(length(daxFiles),1);
               for h=1:length(daxFiles)
                   daxFrames{h} = max(ReadDax(daxFiles{h}),[],3,'verbose',false);
               end
           else
               disp(['fov',num2str(fov,'%03d'),'*','fid.dax',' not found.']);
           end
       end
       daxFrames = cat(3,daxFrames{:});
       CT{handles.id}.dataFrames = daxFrames;
       daxFrames = max(daxFrames,[],3);
   end
   high = CT{handles.id}.parsROI.displayContrastHigh;
   figure(1); clf; imagesc(IncreaseContrast(daxFrames,'high',high)); colormap(gray);
elseif strcmp(currButtonName,'Mark Bad Spot')
   CT{handles.id}.rejectPoints = ...
       cat(1,CT{handles.id}.rejectPoints,...
       CT{handles.id}.parsFOV.lociXY(CT{handles.id}.currentSpotID,:) );
end



% --- Executes on button press in ButtonEditManually2.
function ButtonEditManually2_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonEditManually2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CT 
currButtonName = get(handles.ButtonEditManually2,'String');
if  strcmp(currButtonName,'Record spots')
    fov = str2double(get(handles.EditTextFOV,'String'));
    CT{handles.id}.fovDatas(fov).lociXY = CT{handles.id}.parsFOV.lociXY;
    nSpots = size(CT{handles.id}.parsFOV.lociXY,1);
    disp(['recorded ',num2str(nSpots),' for FOV ',num2str(fov)]);
end


function EditSaveFolder_Callback(hObject, eventdata, handles)
% hObject    handle to EditSaveFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CT
CT{handles.id}.parsFOV.saveFolder = get(handles.EditSaveFolder,'String');
SetFigureSavePath(CT{handles.id}.parsFOV.saveFolder,'makeDir',true);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function EditSaveFolder_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function eTable = ReadExperimentTable(hObject,eventdata,handles)
global CT
% load experiment table
dataTypes = CT{handles.id}.parsLoadData.dataTypes;
CT{handles.id}.expTableXLS = get(handles.EditTextFolder,'String'); 
eTable = readtable(CT{handles.id}.expTableXLS);
hybFolderID = false(height(eTable),1);
if sum(  strcmp(dataTypes,'any') ) == 0
    for i=1:length(dataTypes)
        hybFolderID = hybFolderID | strcmp(eTable.DataType,dataTypes{i});
    end
else
    hybFolderID = true(height(eTable),1);
end
eTable = eTable(hybFolderID,:);
CT{handles.id}.parsFOV.eTable = eTable;    
    

