function varargout = ChrTracer(varargin)
% CHRTRACER MATLAB code for ChrTracer.fig
%      CHRTRACER, by itself, creates a new CHRTRACER or raises the existing
%      singleton*.
%
%      H = CHRTRACER returns the handle to a new CHRTRACER or the handle to
%      the existing singleton*.
%
%      CHRTRACER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHRTRACER.M with the given input arguments.
%
%      CHRTRACER('Property','Value',...) creates a new CHRTRACER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ChrTracer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ChrTracer_OpeningFcn via varargin.
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

% Edit the above text to modify the response to help ChrTracer

% Last Modified by GUIDE v2.5 26-May-2017 14:25:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ChrTracer_OpeningFcn, ...
                   'gui_OutputFcn',  @ChrTracer_OutputFcn, ...
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


% --- Executes just before ChrTracer is made visible.
function ChrTracer_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ChrTracer (see VARARGIN)

% all data for the chromosome tracer is kept inside a structure array;
% so there can be multiple version of ChrTracer open, each version has its
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
CT{handles.id}.datTable = table();
CT{handles.id}.datTableTemp = table(); 
CT{handles.id}.spotTable = table();
CT{handles.id}.im = [];
CT{handles.id}.rejectPoints = [];

set(handles.ButtonRunStep,'String','Load Ex. Table');
% Choose default command line output for ChrTracer
handles.output = hObject;

DefaultParmeters(hObject,eventdata,handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ChrTracer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ChrTracer_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;



function DefaultParmeters(hObject,eventdata,handles)
% These default parameters are passed forth to the indicated function
% They can also be edited in the popup GUI
% 
% 
global CT scratchPath;
% FOV parameters
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true}; 
defaults(end+1,:) = {'veryverbose','boolean',false};
defaults(end+1,:) = {'showPlots', 'boolean', true}; 
defaults(end+1,:) = {'saveData', 'boolean', true}; 
defaults(end+1,:) = {'fov', 'integer', 1}; 
defaults(end+1,:) = {'firstSpot', 'integer', 1};
defaults(end+1,:) = {'lastSpot', 'integer', []};
defaults(end+1,:) = {'lociXY','freeType',[]}; % for plotting previously selected spots 
defaults(end+1,:) = {'stopOnError','boolean',false};
defaults(end+1,:) = {'eTable','array',[]};
defaults(end+1,:) = {'numParallel', 'integer',1};
defaults(end+1,:) = {'saveFolder', 'string', scratchPath};
defaults(end+1,:) = {'dataFolder', 'string', ''};
parsFOV = ParseVariableArguments([],defaults,'ChrTracer');
% Load Data and Register
defaults = cell(0,3);
defaults(end+1,:) = {'filename', 'string', 'ConvZscan'}; 
defaults(end+1,:) = {'dataTypes','cell',{'H','T','R'}};
defaults(end+1,:) = {'alignmentBoxWidth', 'positive', inf}; % pixels.  size of region to use for frame alignment
defaults(end+1,:) = {'alignUpsample', 'positive', 1}; % upsample data for coarse alignment (slow, should be unnecessary);
defaults(end+1,:) = {'alignContrastLow', 'fraction', .8}; % low image threshold for contrast balance prior to coarse alignment
defaults(end+1,:) = {'alignContrastHigh', 'fraction', .999}; % high threshold  for contrast balance prior to coarse alignment
defaults(end+1,:) = {'alignToFirst', 'boolean', true}; % align to first image? if false will align to previous non-empty image data
defaults(end+1,:) = {'minFracObj','fraction',.65}; % min fraction of objects from previous hybe found in target hybe to be acceptable 
defaults(end+1,:) = {'refHybe','integer',1}; % hybe to use to start alignment
defaults(end+1,:) = {'previousAlignFrames','array',[]}; % previous align frames
defaults(end+1,:) = {'rotation','boolean',true}; % also compute and correct rotation angle
defaults(end+1,:) = {'maxD','positive',15}; % distance in pixels to match objects prior to computing rotation
defaults(end+1,:) = {'targetHybes','integer',[]};
defaults(end+1,:) = {'goodHybes','boolean',[]};
defaults(end+1,:) = {'flattenBackground','nonnegative',0};
defaults(end+1,:) = {'corrAngles','freeType',[]}; % rotate angles to check for better alignment
parsLoadData = ParseVariableArguments([],defaults,'ChrTracer_LoadData');
% FOV summary
defaults = cell(0,3);
defaults(end+1,:) = {'contrastLow','fraction',0};
defaults(end+1,:) = {'contrastHigh','fraction',1};
defaults(end+1,:) = {'fidGain','nonnegative',1};
defaults(end+1,:) = {'datGain','nonnegative',1}; 
defaults(end+1,:) = {'saveFig','boolean',true};
defaults(end+1,:) = {'showFidXY','boolean',true};
defaults(end+1,:) = {'showFidXZ','boolean',false};
defaults(end+1,:) = {'showDatXY','boolean',true};
defaults(end+1,:) = {'showDatXZ','boolean',false};
defaults(end+1,:) = {'rescaleTile','positive',.1};
parsFovOverlay = ParseVariableArguments([],defaults,'CHRTracer_FOVsummaryPolots');
% Select ROI
defaults = cell(0,3);
defaults(end+1,:) = {'ROIselectMode',{'manual','auto'},'manual'};
defaults(end+1,:) = {'autoSelectThreshold','fraction',.992};
defaults(end+1,:) = {'autoSelectDownsample','positive',3};
defaults(end+1,:) = {'mergeFrames','string','median'};
defaults(end+1,:) = {'border','nonnegative',10};
defaults(end+1,:) = {'coarseBackgroundScale','nonnegative',0}; % 50
defaults(end+1,:) = {'fineBackgroundScale','nonnegative',0};   % 5
defaults(end+1,:) = {'gain','positive',1};
parsROI = ParseVariableArguments([],defaults,'ChrTracer_ROIselect');
% Plot Spots
defaults = cell(0,3);
defaults(end+1,:) = {'boxWidth', 'positive', 16};
defaults(end+1,:) = {'showFolderNames', 'boolean',false};
defaults(end+1,:) = {'currentSpot','integer',CT{handles.id}.currentSpotID};
defaults(end+1,:) = {'fidGain','nonnegative',.5};
defaults(end+1,:) = {'datGain','nonnegative',.5}; 
parsCrop = ParseVariableArguments([],defaults,'ChrTracer_CropAndPlot');
% Fit Spots
defaults = cell(0,3);
% defaults(end+1,:) = {'batch','boolean',false};  % 
defaults(end+1,:) = {'fidMinPeakHeight', 'positive', 400};
defaults(end+1,:) = {'fidCameraBackground', 'nonnegative', 0};
defaults(end+1,:) = {'fidPeakBlur', 'nonnegative', .5};
defaults(end+1,:) = {'fidTroubleshoot', 'boolean', false};
defaults(end+1,:) = {'fidMaxFitWidth', 'positive', 8};
defaults(end+1,:) = {'fidMinSep', 'nonnegative', 3};
defaults(end+1,:) = {'datMinPeakHeight', 'positive', 800};
defaults(end+1,:) = {'datCameraBackground', 'nonnegative', 0};
defaults(end+1,:) = {'datPeakBlur', 'nonnegative', .5};
defaults(end+1,:) = {'datTroubleshoot', 'boolean', false};
defaults(end+1,:) = {'datMaxFitWidth', 'positive', 6};
defaults(end+1,:) = {'datMinSep', 'nonnegative', 3};
defaults(end+1,:) = {'reloadFigs','boolean',false};
defaults(end+1,:) = {'nmXYpix','positive',154};
defaults(end+1,:) = {'nmZpix','positive',100};
defaults(end+1,:) = {'fidKeepBrightest','integer',1};
defaults(end+1,:) = {'fidRelativeHeight','fraction',0};
defaults(end+1,:) = {'fidMinHBratio','nonnegative',1.2}; % peak value over background value
defaults(end+1,:) = {'fidMinAHratio','nonnegative',.25}; % fitted height over background vs peak value
defaults(end+1,:) = {'fidMaxUncert','nonnegative',2}; % pixels
defaults(end+1,:) = {'datKeepBrightest','integer',1};
defaults(end+1,:) = {'datRelativeHeight','fraction',0};
defaults(end+1,:) = {'datMinHBratio','nonnegative',1.2}; % peak value over background value
defaults(end+1,:) = {'datMinAHratio','nonnegative',.25}; % fitted height over background vs peak value
defaults(end+1,:) = {'datMaxUncert','nonnegative',2}; % pixels
parsFit = ParseVariableArguments([],defaults,'ChrTracer_FitPSFs');
% Link and Render
defaults = cell(0,3);
defaults(end+1,:) = {'renderInterpMethod',{'pchip','spline','linear','nearest','next','previous','cubic'},'pchip'};
defaults(end+1,:) = {'renderTubeRadius','positive',25};
defaults(end+1,:) = {'renderBallRadius','positive',30};
defaults(end+1,:) = {'minFracFidPoints','fraction',.5};
defaults(end+1,:) = {'alignToPrevious','boolean',false};
defaults(end+1,:) = {'refHybe','integer',1};
defaults(end+1,:) = {'fixUnique','boolean',false};
parsLink = ParseVariableArguments([],defaults,'ChrTracer_SpotTableToMap');

% Combine in CT
CT{handles.id}.parsFOV = parsFOV;
CT{handles.id}.parsLoadData = parsLoadData;
CT{handles.id}.parsFovOverlay = parsFovOverlay;
CT{handles.id}.parsROI = parsROI; 
CT{handles.id}.parsCrop = parsCrop; 
CT{handles.id}.parsFit = parsFit;
CT{handles.id}.parsLink = parsLink;

% default Steps
CT{handles.id}.stepNames = {'Load Ex. Table';
                            'Save FOV Overlay';
                            'Adjust Alignment';
                            'Select ROI';
                            'Plot Spots';
                            'Fit PSFs';
                            'Register and Plot';
                            'Accept Fits'};
CT{handles.id}.currStep = 'Load Ex. Table';


% --- Executes on button press in ButtonRunStep.
function ButtonRunStep_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to ButtonRunStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CT

currStepName = CT{handles.id}.currStep;

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
        expTableXLS = get(handles.EditTextFolder,'String'); 
        CT{handles.id}.parsFOV.fov = str2double(get(handles.EditTextFOV,'String'));
        if isempty(expTableXLS)
            error('Please enter an experiment table');
        end
        % pars = CT{handles.id};
        pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsLoadData);
        % now we actually load the table
        [dataOut,parsOut] = ChrTracer_LoadData(expTableXLS,'parameters',pars);
        CT{handles.id}.rawFiducialFrames = dataOut.rawFiducialFrames;
        CT{handles.id}.rawDataFrames = dataOut.rawDataFrames;
        CT{handles.id}.fiducialFrames = dataOut.fiducialFrames;
        CT{handles.id}.dataFrames = dataOut.dataFrames;
        CT{handles.id}.parsFOV.goodHybes = dataOut.goodHybes;
        CT{handles.id}.parsFOV.eTable = dataOut.eTable;
        CT{handles.id}.parsFOV.saveFolder = dataOut.saveFolder;
        CT{handles.id}.parsFOV.dataFolder = dataOut.dataFolder;
        set(handles.EditSaveFolder,'String',CT{handles.id}.parsFOV.saveFolder);
        disp('Experiment Data Loaded'); 
        CT{handles.id}.currStep = 'Save FOV Overlay';
        % Auto advance to next step
    else
        disp('Load table aborted'); 
    end
    
elseif strcmp(currStepName,'Save FOV Overlay')
    pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsFovOverlay);
    CT{handles.id}.im = ChrTracer_FOVsummaryPlots(CT{handles.id}.fiducialFrames,...
        CT{handles.id}.dataFrames,'parameters',pars);
    CT{handles.id}.currStep = 'Adjust Alignment';
    disp('Open parameters to select "target frames" to realign');
    disp('Select Next Step if Alignment already looks good');
    
elseif strcmp(currStepName,'Adjust Alignment')
    pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsLoadData);
    targetHybes = CT{handles.id}.parsLoadData.targetHybes;
    if ~isempty(pars.goodHybes)
        % if goodHybes was left empty, we don't want to pass it to the function
        % this would be interpreted as everything is good.
        CT{handles.id}.parsFOV.goodHybes = pars.goodHybes;  
    end
    if ~isempty(targetHybes)
        [fiducialAlignFrames,dataAlignFrames,goodHybes] = RegisterImages(...
            CT{handles.id}.rawFiducialFrames,CT{handles.id}.rawDataFrames,...
        'parameters',pars,'previousAlignFrames',CT{handles.id}.fiducialFrames,...
        'goodHybes',CT{handles.id}.parsFOV.goodHybes);
        % update 
        CT{handles.id}.fiducialFrames = fiducialAlignFrames;
        CT{handles.id}.dataFrames(targetHybes,:) = dataAlignFrames(targetHybes,:);
        CT{handles.id}.parsFOV.goodHybes = goodHybes;
    end
    figure();
    ref = max(fiducialAlignFrames{CT{handles.id}.parsLoadData.refHybe},[],3);
    new = cellfun(@(x) max(x,[],3),fiducialAlignFrames(targetHybes),'UniformOutput',false);
    Ncolor(  20*cat(3,ref, new{:}) );
    CT{handles.id}.currStep = 'Save FOV Overlay';

elseif strcmp(currStepName,'Select ROI')
    pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsROI);
    spot = ChrTracer_ROIselect(CT{handles.id}.fiducialFrames,...
        'parameters',pars,...
        'im',CT{handles.id}.im,...
        'badSpots',CT{handles.id}.rejectPoints);
    if strcmp(pars.ROIselectMode,'manual')
        CT{handles.id}.parsFOV.lociXY(CT{handles.id}.currentSpotID,:) = spot;
        CT{handles.id}.parsFOV.firstSpot = CT{handles.id}.currentSpotID;
        CT{handles.id}.parsFOV.lastSpot = CT{handles.id}.currentSpotID;
        CT{handles.id}.parsCrop.currentSpot= CT{handles.id}.currentSpotID;
    else
        numSpots =  size(spot,1);
        CT{handles.id}.parsFOV.lociXY = spot;
        CT{handles.id}.parsFOV.numSpots = numSpots;
        CT{handles.id}.parsFOV.firstSpot = 1;
        CT{handles.id}.parsFOV.lastSpot = numSpots;
    end
    CT{handles.id}.currStep = 'Plot Spots'; 
    SetupSteps(hObject,eventdata,handles);
    % Auto advance to next step 

elseif strcmp(currStepName,'Plot Spots')
     pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsCrop);
     if pars.numParallel < 2
         CT{handles.id}.parsFOV.numParallel = 1;
         pars.numParallel = 1;
         pars.firstSpot = CT{handles.id}.parsCrop.currentSpot;
         pars.lastSpot = CT{handles.id}.parsCrop.currentSpot;
     end
     [CT{handles.id}.fidSpots,CT{handles.id}.dataSpots,parsOut] = ...
         ChrTracer_CropSpots(CT{handles.id}.fiducialFrames,CT{handles.id}.dataFrames,...   ChrTracer_CropAndPlot
             'parameters',pars);
     CT{handles.id}.currStep = 'Fit PSFs';
     SetupSteps(hObject,eventdata,handles);
     % Auto advance to next step

elseif strcmp(currStepName,'Fit PSFs')
    pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsFit,'conflict','keepFirst');
    if pars.numParallel < 2 % stepwise mode is for 1 spot at a time. 
        pars.firstSpot = CT{handles.id}.parsCrop.currentSpot;
        pars.lastSpot = CT{handles.id}.parsCrop.currentSpot;
    end
    [spotTableTemp,parsOut] = ChrTracer_FitSpots(...
        CT{handles.id}.fidSpots,CT{handles.id}.dataSpots,...
        'parameters',pars);
    CT{handles.id}.spotTableTemp = spotTableTemp;
    % Requires 'Next Step' to proceed
    
elseif strcmp(currStepName,'Register and Plot')
    pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsLink);
    if pars.numParallel < 2  % stepwise mode is for 1 spot at a time. 
        pars.firstSpot = CT{handles.id}.parsCrop.currentSpot;
        pars.lastSpot = CT{handles.id}.parsCrop.currentSpot;
    end
    dat_table = ChrTracer_UniqueSpotTableToMap(CT{handles.id}.spotTableTemp,'parameters',pars); 
    CT{handles.id}.datTableTemp = dat_table;
    % Requires 'Next Step' to proceed
    
elseif strcmp(currStepName,'Accept Fits')
    disp('saving fits...');
    CT{handles.id}.spotTable = cat(1,CT{handles.id}.spotTable,...
        CT{handles.id}.spotTableTemp);
     CT{handles.id}.datTable = cat(1,CT{handles.id}.datTableTemp,...
         CT{handles.id}.datTableTemp);
    CT{handles.id}.currentSpotID = CT{handles.id}.currentSpotID+1;
    CT{handles.id}.currStep = 'Select ROI';
end

set(handles.ButtonRunStep,'String',CT{handles.id}.currStep);
guidata(hObject, handles);

% --- Executes on button press in ButtonParameters.
function ButtonParameters_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonParameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CT
currStepName = CT{handles.id}.currStep;
if strcmp(currStepName,'Load Ex. Table')
    CT{handles.id}.parsLoadData = ChrTracer_ParameterGUI(CT{handles.id}.parsLoadData);
elseif strcmp(currStepName,'Save FOV Overlay')
    CT{handles.id}.parsFovOverlay = ChrTracer_ParameterGUI(CT{handles.id}.parsFovOverlay); 
elseif strcmp(currStepName,'Adjust Alignment')
    CT{handles.id}.parsLoadData = ChrTracer_ParameterGUI(CT{handles.id}.parsLoadData);
elseif strcmp(currStepName,'Select ROI')
    CT{handles.id}.parsROI = ChrTracer_ParameterGUI(CT{handles.id}.parsROI); 
elseif strcmp(currStepName,'Plot Spots')
    CT{handles.id}.parsCrop = ChrTracer_ParameterGUI(CT{handles.id}.parsCrop); 
elseif strcmp(currStepName,'Fit PSFs')
    ChrTracer_PSFparsGUI(handles.id);
elseif strcmp(currStepName,'Register and Plot')
    ChrTracer_LinkSpotsParsGUI(handles.id);
else
    disp(['parameter GUI not yet written for step ', num2str(CT{handles.id}.step)]);
end




% --- Internal Function, Clear all spot data from current FOV
function ClearFOVdata(handles)
global CT;
CT{handles.id}.parsFOV.lociXY = []; 
CT{handles.id}.currentSpotID = 1; % external, like step
CT{handles.id}.spotTable = table();
CT{handles.id}.im = [];
CT{handles.id}.rejectPoints = [];

% ---- Internal Function, remove last point
function RemoveLocus(hObject,eventdata,handles)
global CT
CT{handles.id}.parsFOV.lociXY(end,:) = [];
CT{handles.id}.stepName = 'Select ROI';
set(handles.ButtonRunStep,'String','Plot Spots');
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
set(handles.ButtonRunStep,'String',CT{handles.id}.currStep);
if strcmp(currStepName,'Fit PSFs')
    set(handles.ButtonEditManually1,'String','Remove Points');
    set(handles.ButtonEditManually2,'String','Add Points');
    set(handles.ButtonEditManually1,'Visible','On');
    set(handles.ButtonEditManually2,'Visible','On');
elseif strcmp(currStepName,'Plot Spots')
    set(handles.ButtonEditManually1,'Visible','On');
    set(handles.ButtonEditManually1,'String','Mark Bad Spot');
    set(handles.ButtonEditManually2,'Visible','Off');
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
elseif strcmp(currStepName,'Save FOV Overlay')
    CT{handles.id}.parsFovOverlay = ChrTracer_ParameterGUI(CT{handles.id}.parsFovOverlay); 
elseif strcmp(currStepName,'Adjust Alignment')
    CT{handles.id}.parsLoadData = ChrTracer_ParameterGUI(CT{handles.id}.parsLoadData);
elseif strcmp(currStepName,'Select ROI')
    CT{handles.id}.parsROI = ChrTracer_ParameterGUI(CT{handles.id}.parsROI); 
elseif strcmp(currStepName,'Plot Spots')
    CT{handles.id}.parsCrop = ChrTracer_ParameterGUI(CT{handles.id}.parsCrop); 
elseif strcmp(currStepName,'Fit PSFs')
    ChrTracer_PSFparsGUI(handles.id);
elseif strcmp(currStepName,'Register and Plot')
    CT{handles.id}.parsLink = ChrTracer_ParameterGUI(CT{handles.id}.parsLink); 
    % ChrTracer_LinkSpotsParsGUI(handles.id);
else
    disp(['parameter GUI not yet written for step ', num2str(CT{handles.id}.step)]);
end

% --------------------------------------------------------------------
function MenuAnalyzeAll_Callback(hObject, eventdata, handles)
% Analyze all spots in current FOV
% hObject    handle to MenuAnalyzeAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global CT

% Typically an experiment table will already be loaded.  If it is not, we
% load one. 
if isempty(CT{handles.id}.parsFOV.eTable)
    expTableXLS = get(handles.EditTextFolder,'String'); 
    CT{handles.id}.parsFOV.fov = str2double(get(handles.EditTextFOV,'String'));
    if isempty(expTableXLS)
        error('Please enter an experiment table');
    end
    pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsLoadData);
    % now we actually load the table
    [dataOut,parsOut] = ChrTracer_LoadData(expTableXLS,'parameters',pars);
    CT{handles.id}.fiducialFrames = dataOut.fiducialFrames;
    CT{handles.id}.dataFrames = dataOut.dataFrames;
    CT{handles.id}.parsFOV.goodHybes = dataOut.goodHybes;
    disp('Experiment Data Loaded'); 
end

% Confirm that we want to clear the existing spot tables. 
if ~isempty(CT{handles.id}.parsFOV.lociXY) 
   cprintf('red','This will overwrite existing spot data');
   choice = input('Do you wish to continue? 1=yes, 0=no?  ');
else 
    choice = true;
end
    
if choice
    % Add some parallel processors
    CT{handles.id}.parsFOV.numParallel = 16;
    CT{handles.id}.parsROI.ROIselectMode = 'auto';
    
    % Save FOV summary plots
    fovFile = [CT{handles.id}.parsFOV.saveFolder,'fov',...
        num2str(CT{handles.id}.parsFOV.fov,'%03d'),'_fid_overlayFig.png'];
    if exist(fovFile,'file')~=2
        disp('saving FOV summary plots');
    pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsFovOverlay);
    CT{handles.id}.im = ChrTracer_FOVsummaryPlots(CT{handles.id}.fiducialFrames,...
        CT{handles.id}.dataFrames,'parameters',pars);
    end
    % autoROI select
    disp('finding spots in current fov');
    pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsROI);
    spot = ChrTracer_ROIselect(CT{handles.id}.fiducialFrames,'parameters',pars,'im',CT{handles.id}.im);
    numSpots =  size(spot,1);
    CT{handles.id}.parsFOV.lociXY = spot;
    CT{handles.id}.parsFOV.numSpots = numSpots;
    CT{handles.id}.parsFOV.firstSpot = 1;
    CT{handles.id}.parsFOV.lastSpot = numSpots;
    
    % Crop
    disp('cropping spots');
     pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsCrop);
     pars.firstSpot = 1; 
     pars.lastSpot = numSpots;
     [CT{handles.id}.fidSpots,CT{handles.id}.dataSpots] = ChrTracer_CropSpots(...
             CT{handles.id}.fiducialFrames,CT{handles.id}.dataFrames,...
             'parameters',pars);

    % Fit 
    pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsFit,'conflict','keepFirst');
    % ----------------------------------------------------% 
    pars.showPlots = false;  % prevent memory overload and increase speed
    % -----------------------------------------------------% 
    spotTableTemp = ChrTracer_FitSpots(...
        CT{handles.id}.fidSpots,CT{handles.id}.dataSpots,...
        'parameters',pars);
    % Store in global
    CT{handles.id}.spotTable = spotTableTemp;
    
     % register
     pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsLink,'conflict','keepFirst');
     CT{handles.id}.datTable = ChrTracer_UniqueSpotTableToMap(CT{handles.id}.spotTable,'parameters',pars); 
    
    % save data
    MenuSaveFOVdata_Callback(hObject, eventdata, handles)
end




% --------------------------------------------------------------------
function MenuAnalyzeAllFOV_Callback(hObject, eventdata, handles)
% hObject    handle to MenuAnalyzeAllFOV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CT

% set up some parallel processing
CT{handles.id}.parsFOV.numParallel = 20;
CT{handles.id}.parsROI.ROIselectMode = 'auto';

% load experiment table
dataTypes = CT{handles.id}.parsLoadData.dataTypes;
expTableXLS = get(handles.EditTextFolder,'String'); 
eTable = readtable(expTableXLS);
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

% determine number of FOV in data folder
dataFolder = fileparts(expTableXLS);
currFolder = [dataFolder,filesep,eTable.FolderName{1},filesep];
daxFiles =  cellstr(ls([currFolder,'ConvZScan','*.dax']));
daxNums = 1:length(daxFiles);

% check for existing files
saveFolder = CT{handles.id}.parsFOV.saveFolder;
foundSpotTables = ls([saveFolder,'fov*datTable.csv']);
if ~isempty(foundSpotTables)
    fovComplete = str2num(foundSpotTables(:,4:6))'; %#ok<ST2NM>
else
    fovComplete = [];
end

% skip FOVs which already have ctFiles
% (later we can modify this with parameters for overwrite)
fin = intersect(daxNums,fovComplete);
if ~isempty(fin)
    disp(['found data for FOVs ',num2str(fin)]);
    disp('skipping these FOVs');
end
fovs = setdiff(daxNums,fovComplete);

for f=fovs
    
    CT{handles.id}.parsFOV.fov = f;
    pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsLoadData);
    % now we actually load the table
    dataOut = ChrTracer_LoadData(expTableXLS,'parameters',pars);
    CT{handles.id}.fiducialFrames = dataOut.fiducialFrames;
    CT{handles.id}.dataFrames = dataOut.dataFrames;
    CT{handles.id}.parsFOV.goodHybes = dataOut.goodHybes;
    disp('Experiment Data Loaded'); 
    
    % Save FOV summary plots if not already saved
    fovFile = [CT{handles.id}.parsFOV.saveFolder,'fov',...
        num2str(CT{handles.id}.parsFOV.fov,'%03d'),'_fid_overlayFig.png'];
    if exist(fovFile,'file')~=2
        pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsFovOverlay);
        CT{handles.id}.im = ChrTracer_FOVsummaryPlots(CT{handles.id}.fiducialFrames,...
            CT{handles.id}.dataFrames,'parameters',pars);
    end
    % autoROI select
    pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsROI);
    [spot] = ChrTracer_ROIselect(CT{handles.id}.fiducialFrames,'parameters',pars,'im',CT{handles.id}.im);
    numSpots =  size(spot,1);
    CT{handles.id}.parsFOV.lociXY = spot;
    CT{handles.id}.parsFOV.numSpots = numSpots;
    CT{handles.id}.parsFOV.firstSpot = 1;
    CT{handles.id}.parsFOV.lastSpot = numSpots;
    
    % Crop
    % skip if all spots already exist
    cropTiles = cellstr( ls([CT{handles.id}.parsFOV.saveFolder,num2str(CT{handles.id}.parsFOV.fov),'*_dat_tileXY.fig']) );
    if length(cropTiles) < CT{handles.id}.parsFOV.numSpots
         pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsCrop);
         pars.firstSpot = 1; 
         pars.lastSpot = numSpots;
         [CT{handles.id}.fidSpots,CT{handles.id}.dataSpots] = ChrTracer_CropSpots(...
                 CT{handles.id}.fiducialFrames,CT{handles.id}.dataFrames,...
                 'parameters',pars);
    end

    % clear large data
    CT{handles.id}.fiducialFrames = [];
    CT{handles.id}.dataFrames = [];
    CT{handles.id}.rawFiducialFrames = [];
    CT{handles.id}.rawDataFrames = [];
    
    % Fit
    % skip if all spots already exits
    spottables = cellstr( ls([CT{handles.id}.parsFOV.saveFolder,num2str(CT{handles.id}.parsFOV.fov),'*_spottable.csv']) );
    if length(spottables) < CT{handles.id}.parsFOV.numSpots
        pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsFit,'conflict','keepFirst');
        %------------------------------------------------------% 
        pars.showPlots = false;  % prevent memory overload and increase speed
        %------------------------------------------------------% 
        [spotTableTemp] = ChrTracer_FitSpots(...
            CT{handles.id}.fidSpots,CT{handles.id}.dataSpots,...
            'parameters',pars);
    else
        allTables = cell(numSpots,1);
        for s=1:numSpots
            allTables{s} = readtable([CT{1}.parsFOV.saveFolder, spottables{s}]);
        end
        spotTableTemp = cat(1,allTables{:});
    end
    
    % Store in global
    CT{handles.id}.spotTable = spotTableTemp;
     
     % register
     pars = JoinStructures(CT{handles.id}.parsFOV,CT{handles.id}.parsLink,'conflict','keepFirst');
     CT{handles.id}.datTable = ChrTracer_UniqueSpotTableToMap(CT{handles.id}.spotTable,'parameters',pars); 
   
    % save data
    MenuSaveFOVdata_Callback(hObject, eventdata, handles);
    
    % close down parpool to free up memory;
    poolobj = gcp('nocreate');
    delete(poolobj);
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
    load([pathName,fileName]);
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
ctData.dataFrames = [];
ctData.rawFiducialFrames = [];
ctData.rawDataFrames = [];
saveFolder = ctData.parsFOV.saveFolder;
fov = ctData.parsFOV.fov;
disp(['saving data as ',saveFolder,'fov',num2str(fov,'%03d'),'_ctData.mat']);
save([saveFolder,'fov',num2str(fov,'%03d'),'_ctData.mat'],'ctData'); 

tableName = [saveFolder,'fov',num2str(fov,'%03d'),'_AllSpotTable.csv'];
disp(['Writing table: ',tableName]); 
writetable(ctData.spotTable,tableName); 

tableName = [saveFolder,'fov',num2str(fov,'%03d'),'_datTable.csv'];
disp(['Writing table: ',tableName]); 
writetable(ctData.datTable,tableName); 
disp('Data saved'); 



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
    CT{handles.id} = ChrTracer_LoadData(expTableXLS,'parameters',CT{handles.id});
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
% hObject    handle to EditTextFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
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
CT{handles.id}.currStep = 'Load Ex. Table';
set(handles.ButtonRunStep,'String',CT{handles.id}.currStep);
CT{handles.id}.parsLoadData.targetHybes = []; % clear target hybes
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function EditTextFOV_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditTextFOV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [spots,hybs,axesHandles,figHandle,compFigHandle,compAxeHandles,isfid,isXY] = GetSpotData()
% Return user input selection of spot

    % This part should become an external function
    [spots,hybs,axesHandles] = GetPixelInfo();
    figHandle = gcf;
    
    if ~isempty(strfind(figHandle.Name,'dat_tile'))
        isfid = false;
        dataType = 'dat_tile';
    elseif ~isempty(strfind(figHandle.Name,'fid_tile'))
        isfid = true;
        dataType = 'fid_tile';
    else 
        warning('could not determine data type'); 
    end
    
    if ~isempty(strfind(figHandle.Name,'_tileXZ'))
        isXY =false;
    elseif ~isempty(strfind(figHandle.Name,'_tileXY'))
        isXY = true;
    else 
        warning('could not determine projection type'); 
    end
    
    allFigs = findobj('Type','Figure');
    if isXY
        figId = StringFind({allFigs.Name},dataType,'boolean',true) & ...
            StringFind({allFigs.Name},'_tileXZ','boolean',true); 
    else
        figId = StringFind({allFigs.Name},dataType,'boolean',true) & ...
            StringFind({allFigs.Name},'_tileXY','boolean',true); 
    end
    compFigHandle = allFigs(figId);
    if length(compFigHandle) > 1
        cprintf([1 0 0],'Error: found multiple figures with the same name');
        disp(compFigHandle);
        error('Close extra figure and try again. If fails, close all figures and rerun Fit PSFS');        
    end
        
    numHybes = length(figHandle.Children);
    compAxeHandles = cell(length(hybs),1);
    for h=1:length(hybs)
        compAxeHandles{h} = compFigHandle.Children(numHybes+1-hybs(h));
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
currStepName = CT{handles.id}.currStep;
if strcmp(currStepName,'Plot Spots')
   CT{handles.id}.rejectPoints = ...
       cat(1,CT{handles.id}.rejectPoints,...
       CT{handles.id}.parsFOV.lociXY(CT{handles.id}.currentSpotID,:) );
end


if strcmp(currStepName,'Fit PSFs')
    
    [spots,hybs,axesHandles,figHandle,compFigHandle,compAxeHandles,isfid,isXY] = GetSpotData();
    % spots - x,y, coordinate in selected subplot
    % hybs - number of subplot which was selected
    % axesHandles - handle of selected subplot
    % figHandle - handle of selected figure
    % compFigHandle - handle of companion projection figure (e.g. xz projection if xy figure was selected)
    % compAxesHandles - handle to companion subplot in companion figure
    % isfid - true if selected figure was the fiducial plot. otherwise false 
    % isXY - true if selected figure was the XY projection. otherwise false 
    
    % targetSpots = false(length(hybs),1);
    spotTableTemp = CT{handles.id}.spotTableTemp;
    targetSpots = false(height(CT{handles.id}.spotTableTemp),1); 
    
    % remove all points in selected hybe (panel)
    figure(figHandle);
    for n=1:length(axesHandles)
        allLines = findobj(axesHandles{n},'Type','Line');
        for i=1:length(allLines)
            allLines(i).Color = 'r';
        end
        figure(compFigHandle);
        allLines = findobj(compAxeHandles{n},'Type','Line');
        for i=1:length(allLines)
            allLines(i).Color = 'r'; % recolor points red if removed.
        end
        % delete(allLines); % remove points
        delPoint = find(spotTableTemp.panel == hybs(n) & spotTableTemp.isfid == isfid);
        if isempty(delPoint)
            disp(['no spots to remove from hyb ',num2str(hybs(n))]);
        elseif length(delPoint) > 1
            disp(['multiple points returned in hyb ',num2str(hybs(n)), ' is an error.']);
        else
            targetSpots(delPoint) = true;
        end
    end
    disp(spotTableTemp(targetSpots,:));
    disp(['deleting entry for hybs ',num2str(hybs') ' from spotTableTemp']);
    CT{handles.id}.spotTableTemp(targetSpots,:) = [];
    
        
    if CT{handles.id}.parsFOV.saveData
       SaveFigure(figHandle,'formats',{'png','fig'},'overwrite',true); 
    end
end

% --- Executes on button press in ButtonEditManually2.
function ButtonEditManually2_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonEditManually2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CT 
currStepName = CT{handles.id}.currStep;

  %================ ADD SPOTS =====================
if strcmp(currStepName,'Fit PSFs')
    % parameter parsing
    s = CT{handles.id}.currentSpotID;
    pars.lociXY = CT{handles.id}.parsFOV.lociXY;
    pars.fov = CT{handles.id}.parsFOV.fov;
    pars.goodHybes = CT{handles.id}.parsFOV.goodHybes;
    pars.eTable = CT{handles.id}.parsFOV.eTable;
    parsFit = CT{handles.id}.parsFit;
    pars.nmXYpix = parsFit.nmXYpix;
    pars.nmZpix = parsFit.nmZpix;
    
    if isempty(pars.goodHybes)
        [~,~,~,numHybes] = size(CT{handles.id}.dataFrames);
        pars.goodHybes = true(numHybes,1);
    end

    numGoodHybes = sum(pars.goodHybes);
    hybNum = find(pars.goodHybes);
    if ~isempty(pars.eTable)
        dataType = pars.eTable.DataType;
        bitNum = pars.eTable.Bit;
    end

    [spots,hybs,axesHandles,figHandle,compFigHandle,compAxeHandles,isfid,isXY] = GetSpotData();
   
    if isfid
        spotStack = cat(4,CT{handles.id}.fidSpots{s}{pars.goodHybes});
        % load fit parameters from Step Pars
        pars.minPeakHeight = parsFit.fidMinPeakHeight;
        pars.peakBlur = parsFit.fidPeakBlur;
        pars.cameraBackground = parsFit.fidCameraBackground;
        pars.maxFitWidth = parsFit.fidMaxFitWidth;
        pars.minSep = parsFit.fidMinSep;
        pars.keepBrightest = parsFit.fidKeepBrightest;
        pars.relativeHeight = parsFit.fidRelativeHeight;
        pars.minHBratio = parsFit.fidMinHBratio;
        pars.minAHratio = parsFit.fidMinAHratio;
        pars.maxUncert = parsFit.fidMaxUncert;
        pars.xyUnitConvert = pars.nmXYpix;
        pars.zUnitConvert = pars.nmZpix;
        pars.troubleshoot = parsFit.fidTroubleshoot;
    else
        spotStack = cat(4,CT{handles.id}.dataSpots{s}{pars.goodHybes});
        % load fit parameters from Step Pars
        pars.minPeakHeight = parsFit.datMinPeakHeight;
        pars.peakBlur = parsFit.datPeakBlur;
        pars.cameraBackground = parsFit.datCameraBackground;
        pars.maxFitWidth = parsFit.datMaxFitWidth;
        pars.minSep = parsFit.datMinSep;
        pars.keepBrightest = parsFit.datKeepBrightest;
        pars.relativeHeight = parsFit.datRelativeHeight;
        pars.minHBratio = parsFit.datMinHBratio;
        pars.minAHratio = parsFit.datMinAHratio;
        pars.maxUncert = parsFit.datMaxUncert;
        pars.xyUnitConvert = pars.nmXYpix;
        pars.zUnitConvert = pars.nmZpix;
        pars.troubleshoot = parsFit.datTroubleshoot;
    end
    
    for n=1:length(hybs)
        spot = round(spots(n,:));
        hyb = hybNum(hybs(n));
        if isXY
            prof1 = squeeze(spotStack(spot(2),spot(1),:,hyb));  % figure(11); clf; plot(prof1);
            [~,z_i] = max(prof1); 
            x_i = spot(1); 
            y_i = spot(2);
        else
            prof1 = squeeze(spotStack(:,spot(1),spot(2),hyb));   % figure(11); clf; plot(prof1);
            [~,y_i] = max(prof1);
            x_i = spot(1);
            z_i = spot(2); 
        end
       
        im3D = squeeze(spotStack(:,:,:,hyb));
        dTable = FitPsf3D(im3D,'seedPoint',[x_i,y_i,z_i],'parameters',pars);

        if isXY
            figure(figHandle);
            subplot(axesHandles{n});
            hold on; plot(dTable.x/pars.nmXYpix,dTable.y/pars.nmXYpix,'wo','MarkerSize',15);
            figure(compFigHandle);
            subplot(compAxeHandles{n});
            hold on; plot(dTable.x/pars.nmXYpix,dTable.z/pars.nmZpix,'wo','MarkerSize',15);
        else
            figure(compFigHandle);
            subplot(compAxeHandles{n});
            hold on; plot(dTable.x/pars.nmXYpix,dTable.y/pars.nmXYpix,'wo','MarkerSize',15);
            figure(figHandle);
            subplot(axesHandles{n});
            hold on; plot(dTable.x/pars.nmXYpix,dTable.z/pars.nmZpix,'wo','MarkerSize',15);
        end


        % keep track of hybe, spot number, and is fiduial. 
        dTable.panel = hybs(n)*ones(length(dTable.x),1);
        dTable.hybe = hyb*ones(length(dTable.x),1);
        dTable.dataType = repmat(dataType{hyb},length(dTable.x),1);
        dTable.bit = repmat(bitNum(hyb),length(dTable.x),1);
        dTable.s = s*ones(length(dTable.x),1);
        dTable.isfid = isfid & true(length(dTable.x),1);
        dTable.fov = pars.fov*ones(length(dTable.x),1);
        dTable.locusX = pars.lociXY(s,1)*ones(length(dTable.x),1);
        dTable.locusY = pars.lociXY(s,2)*ones(length(dTable.x),1);
        disp(dTable)

        % Add to spotTableTemp
        CT{handles.id}.spotTableTemp = cat(1,CT{handles.id}.spotTableTemp,dTable);
    end
    
    if CT{handles.id}.parsFOV.saveData
       SaveFigure(figHandle,'formats',{'png','fig'},'overwrite',true); 
    end
    
end




function EditSaveFolder_Callback(hObject, eventdata, handles)
% hObject    handle to EditSaveFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CT
CT{handles.id}.parsFOV.saveFolder = get(handles.EditSaveFolder,'String');
SetFigureSavePath(CT{handles.id}.parsFOV.saveFolder,'makeDir',true);
guidata(hObject,handles);

% Hints: get(hObject,'String') returns contents of EditSaveFolder as text
%        str2double(get(hObject,'String')) returns contents of EditSaveFolder as a double


% --- Executes during object creation, after setting all properties.
function EditSaveFolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditSaveFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
