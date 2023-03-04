function varargout = ChrTracer_SpotSelector(varargin)
% CHRTRACER_SPOTSELECTOR MATLAB code for ChrTracer_SpotSelector.fig
%      CHRTRACER_SPOTSELECTOR, by itself, creates a new CHRTRACER_SPOTSELECTOR or raises the existing
%      singleton*.
%
%      H = CHRTRACER_SPOTSELECTOR returns the handle to a new CHRTRACER_SPOTSELECTOR or the handle to
%      the existing singleton*.
%
%      CHRTRACER_SPOTSELECTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHRTRACER_SPOTSELECTOR.M with the given input arguments.
%
%      CHRTRACER_SPOTSELECTOR('Property','Value',...) creates a new CHRTRACER_SPOTSELECTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ChrTracer_SpotSelector_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ChrTracer_SpotSelector_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ChrTracer_SpotSelector

% Last Modified by GUIDE v2.5 13-Jun-2019 14:17:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ChrTracer_SpotSelector_OpeningFcn, ...
                   'gui_OutputFcn',  @ChrTracer_SpotSelector_OutputFcn, ...
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


% --- Executes just before ChrTracer_SpotSelector is made visible.
function ChrTracer_SpotSelector_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ChrTracer_SpotSelector (see VARARGIN)
global CT
% Choose default command line output for ChrTracer_SpotSelector
handles.output = hObject;
% parse default parameters
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true}; 
defaults(end+1,:) = {'veryverbose', 'boolean', false}; 
defaults(end+1,:) = {'id', 'integer', 1}; 
defaults(end+1,:) = {'fov', 'integer', 1}; 
defaults(end+1,:) = {'hybsToLoad', 'integer', 1}; 
defaults(end+1,:) = {'dataType',{'fiducial','data'},'fiducial'};
defaults(end+1,:) = {'hybsToLoad', 'integer', 1}; 
pars = ParseVariableArguments(varargin,defaults,mfilename); 
% update handles structure with id to match CT
handles.id = pars.id;
set(handles.TextID,'String',['id ',num2str(handles.id)]);
% set defaults in edit boxes
set(handles.EditFOV,'String',num2str(pars.fov));
set(handles.EditDispHybs,'String',num2str(pars.hybsToLoad));
CT{handles.id}.parsSpotSelector = pars;
CT{handles.id}.parsSpotSelector.bkd = [];
% load selected spot pars from CT
importPars = CT{handles.id}.parsSelectAllSpots;
if isfield(importPars,'fov')
    importPars = rmfield(importPars,'fov');
end
CT{handles.id}.parsSpotSelector.autofitPars  = importPars; 
LoadImage(handles); % load image into SpotSelector
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes ChrTracer_SpotSelector wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ChrTracer_SpotSelector_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;


  

function imHybs=LoadImage(handles)
%  need speed test. Hopefully this load from disk approach is fast enough
%  to be convienent. Otherwise we can revert to using the information
%  loaded in RAM in CT (though I'd prefer to purge that info). 
    global CT
    set(handles.TextDir,'String','loading data, please wait...');  pause(.1); % guidata(handles.output, handles);
    eTableXLS = CT{handles.id}.expTableXLS;
    % update 'Source Data' Menu
    eTable = readtable(eTableXLS);
    [fidChn,dataChns] = GetChnNamesFromTable(eTable); % building names
    datLabeled = strcat('Chn ', dataChns',' Hybs:'); 
    sourceDat = cat(1,{['Fid. ',fidChn{1},' Hybs:']},{'All Data Hybs:'}, datLabeled{:});
    set(handles.popupmenu1,'String',sourceDat); % updates
    % find hybs to load
    hybsToLoad = CT{handles.id}.parsSpotSelector.hybsToLoad;
    if isinf(hybsToLoad)
        hybsToLoad = 1:height(eTable); 
    end
    fov = CT{handles.id}.parsSpotSelector.fov;
    imHybs = cell(length(hybsToLoad),1);
    try
        % load all requested hybs of indicated dataType for this FOV
        % This dax loader will auto align the data
        temp = LoadDaxFromEtable(eTableXLS,...
                    'dataType',CT{handles.id}.parsSpotSelector.dataType,...
                    'daxRootDefault',CT{handles.id}.parsLoadExpTable.daxRootDefault,...
                    'fov',fov,...
                    'hybNumber',hybsToLoad,...
                    'fixDrift',true,...
                    'driftFolder',CT{handles.id}.saveFolder,...
                    'simplifyOutput',false);
        selDataChns = get(handles.popupmenu1,'Value')-2;
        if selDataChns==-1 % fiducial
            selDataChns = 1; % it will be the only channel
        elseif selDataChns == 0 % all
            selDataChns = 1:size(temp,3);
        end
        imHybs = squeeze(temp(:,fov,selDataChns));

        % max project these hybs. 
        if length(imHybs) > 1
            im = max(cat(3,imHybs{:}),[],3);
        else
            im = imHybs{1}; % this is faster than cat on nothing. 
        end
        CT{handles.id}.currFOVim = im;
        set(handles.TextDir,'String','image loaded'); % guidata(hObject, handles);        
        high = CT{handles.id}.parsSpotSelector.autofitPars.displayContrastHigh; % CT{handles.id}.parsSpotSelector.autofitPars.imContrastHigh;
        low = CT{handles.id}.parsSpotSelector.autofitPars.displayContrastLow; % CT{handles.id}.parsSpotSelector.autofitPars.imContrastLow;
        dispIm = IncreaseContrast(im,'high',high,'low',low);
        figure(1); clf; imagesc(dispIm); colormap(gray); colorbar;  
    catch er
        textDir = {['failed to load data from fov ',num2str(fov),'.  '],
                    'Maybe you have finished processing all FOV.',
                    'Exit ChrTracer_SpotSelector to return to ChrTracer'};
        set(handles.TextDir,'String',textDir); % guidata(hObject, handles);
        disp(er.message);
    end

% --- Executes on button press in ButtonPars.
function ButtonPars_Callback(hObject, eventdata, handles)
global CT
CT{handles.id}.parsSpotSelector.autofitPars = SimpleParameterGUI(CT{handles.id}.parsSpotSelector.autofitPars);


% --- Executes on button press in ButtonFindSpots.
function ButtonFindSpots_Callback(hObject, eventdata, handles)
    global CT
    % import parameters
    pars = CT{handles.id}.parsSpotSelector.autofitPars;
    if ~isfield(CT{handles.id},'currFOVim')
        LoadImage(handles);
    end
    im = CT{handles.id}.currFOVim;
    f = CT{handles.id}.parsSpotSelector.fov;
    % check for background detect
    backgroundCorrect=CT{handles.id}.parsSpotSelector.autofitPars.backgroundCorrect;
    switch backgroundCorrect
        case 'none'
            CT{handles.id}.parsSpotSelector.bkd = [];
        otherwise % case 'median'
            if isempty(CT{handles.id}.parsSpotSelector.bkd)
                bkd = CT_FlattenBackground(handles);
            else
               bkd = CT{handles.id}.parsSpotSelector.bkd;
            end
            im = makeuint(double(im)./bkd,16);
    end
    % apply contrast
    high = CT{handles.id}.parsSpotSelector.autofitPars.displayContrastHigh; % CT{handles.id}.parsSpotSelector.autofitPars.imContrastHigh;
    low = CT{handles.id}.parsSpotSelector.autofitPars.displayContrastLow; % CT{handles.id}.parsSpotSelector.autofitPars.imContrastLow;
    dispIm = IncreaseContrast(im,'high',high,'low',low);
    figure(1); clf; imagesc(dispIm); colormap(gray); colorbar; 
    % check if spots already exist   
    fileName = ['fov',num2str(f,'%03d'),'_selectSpots.csv'];
    saveName = [CT{handles.id}.saveFolder,fileName];
    if exist(saveName,'file')
        answer = questdlg('found a previous spot map in save folder.', ... % question
                        'Skip', ...  % pop-up label
                        'Load','Recompute','Load'); % op1 op2 default
        if strcmp(answer,'Load')
            spots = readtable(saveName);
            spots = spots{:,:};
            CT{handles.id}.lociXY = spots;
%             figure(1); hold on; % loading is handled below
%             plot(spots(:,1),spots(:,2),'yo');
%             text(spots(:,1)+2,spots(:,2),cellstr(num2str( (1:size(spots,1))')),'color','w');
        else
            CT{handles.id}.lociXY = AutoSelectSpots(im,'parameters',pars,'showPlots',false,'numberSpots',true);
        end
    else
        CT{handles.id}.lociXY = AutoSelectSpots(im,'parameters',pars,'showPlots',false,'numberSpots',true);      
    end
    figure(1); title(['fov ',num2str(f)]);
    spots = CT{handles.id}.lociXY;
    hold on; plot(spots(:,1),spots(:,2),'yo','MarkerSize',10);
    text(spots(:,1)+12,spots(:,2),cellstr(num2str( (1:size(spots,1))')),'color','w');
    if pars.densityFilterNeibs > 0
        [~,dist] = knnsearch(spots,spots,'K',10);
        distScore = false(size(dist));
        distScore(dist<pars.densityFilterMinDist) = true;
        badSpts = sum(distScore,2)>pars.densityFilterNeibs;
        plot(spots(badSpts,1),spots(badSpts,2),'ro','MarkerSize',10);
        CT{handles.id}.lociXY(badSpts,:) = [];
    end
    
    
    

    
function bkd = CT_FlattenBackground(handles)
    global CT
    % get all files
    backgroundCorrect=CT{handles.id}.parsSpotSelector.autofitPars.backgroundCorrect;
    orig = LoadDaxFromEtable(CT{handles.id}.expTableXLS,...
                        'dataType',CT{handles.id}.parsSpotSelector.dataType,...
                        'fov',inf,...
                        'hybNumber',1);
    [flat,bkd] = FlattenBackground(orig,'backgroundCorrect',backgroundCorrect,'showPlots',true);
    CT{handles.id}.parsSpotSelector.bkd = bkd;
    
% -- save 
function SaveSpotTable(lociXY,f,handles)
    global CT
    fileName = ['fov',num2str(f,'%03d'),'_selectSpots.csv'];
    saveName = [CT{handles.id}.saveFolder,fileName];
    saveTable = table(lociXY(:,1),lociXY(:,2));
    saveTable.Properties.VariableNames = {'locusX','locusY'};
    answer = 'Yes';
    if exist(saveName,'file')
         answer = questdlg('Overwrite existing spot table?', ... % question
                        'Overwrite?', ...  % pop-up label
                        'Yes','No','Yes'); % op1 op2 default   
    end
    if strcmp(answer,'Yes')
        writetable(saveTable,saveName);
        disp(['wrote ',saveName]);
    elseif strcmp(answer,'No') % load existing table into memory
        lociXY = readtable(saveName);
        CT{handles.id}.lociXY=[lociXY{:,1},lociXY{:,2}];
        disp('loaded previously saved spot table into memory');
    end
    
    
% --- Executes on button press in ButtonSaveSpots.
function ButtonSaveSpots_Callback(hObject, eventdata, handles)
global CT
f = CT{handles.id}.parsSpotSelector.fov;
SaveSpotTable(CT{handles.id}.lociXY,f,handles);
disp('data saved, advancing to next FOV...');
ButtonNext_Callback(hObject, eventdata, handles);

% --- Executes on button press in ButtonAdd.
function ButtonAdd_Callback(hObject, eventdata, handles)
    global CT
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

% --- Executes on button press in ButtonRemove.
function ButtonRemove_Callback(hObject, eventdata, handles)   
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
    RemoveSpots(handles,spotsToRemove);

% --- Executes on button press in ButtonRemoveInArea.
function ButtonRemoveInArea_Callback(hObject, eventdata, handles)
    % update plot and lociXY list
    global CT;
    % spotsToRemove = SpotsInPolygon();
%     high = CT{handles.id}.parsSpotSelector.autofitPars.imContrastHigh;
%     low = CT{handles.id}.parsSpotSelector.autofitPars.imContrastLow;
%     dispIm = IncreaseContrast(CT{handles.id}.currFOVim,'high',high,'low',low);
    fig = figure(1); 
    regIDs = SegmentationGUI(fig,...
        'MarkerTags',{'remove'},...
        'Points',CT{handles.id}.lociXY,...
        'MarkerShape','o','MarkerSize',5,'MarkerDefaultColor',[1 1 0]); 
    spotsToRemove = ~cellfun(@isempty,regIDs);
    RemoveSpots(handles,spotsToRemove);
    
% --- Executes on button press in ButtonRemoveIDs.
function ButtonRemoveIDs_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonRemoveIDs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global CT
    spotsToRemove = str2num(get(handles.EditRemoveIDs,'String')); %#ok<ST2NM>
    CT{handles.id}.parsSpotSelector.removeIDs  =spotsToRemove;
    fig = figure(1);
    RemoveSpots(handles,spotsToRemove);

%  -- internal function, 
%  remove spots from fig and CT list in memory
function RemoveSpots(handles,spotsToRemove)
    global CT
    lociXY = CT{handles.id}.lociXY;
    lociXY(spotsToRemove,:) = [];
    LoadImage(handles); hold on;
    plot(lociXY(:,1),lociXY(:,2),'yo'); 
    text(lociXY(:,1)+2,lociXY(:,2),cellstr(num2str( (1:size(lociXY,1))')),'color','w');
    CT{handles.id}.lociXY = lociXY;
    

% --- Executes on button press in ButtonNext.
function ButtonNext_Callback(hObject, eventdata, handles)
global CT
f = CT{handles.id}.parsSpotSelector.fov;
f = f +1;
set(handles.EditFOV,'String',num2str(f));
CT{handles.id}.parsSpotSelector.fov = f;
guidata(hObject, handles);
LoadImage(handles);


% --- Executes on button press in ButtonPrevious.
function ButtonPrevious_Callback(hObject, eventdata, handles)
global CT
f = CT{handles.id}.parsSpotSelector.fov;
if f>1
    f = f -1;
else
    disp('already at fov 1');
end
set(handles.EditFOV,'String',num2str(f));
CT{handles.id}.parsSpotSelector.fov = f;
guidata(hObject, handles);
LoadImage(handles);




% --- Triggered by modification of Edit box EditFOV
function EditFOV_Callback(hObject, eventdata, handles)
% updates FOV upon edit of field FOV
global CT
CT{handles.id}.parsSpotSelector.fov = str2double(get(handles.EditFOV,'String'));
LoadImage(handles);

% --- Triggered by modification of Edit box EditRemoveIDs
function EditRemoveIDs_Callback(hObject, eventdata, handles)
% updates on modification of RemoveIDs
global CT
CT{handles.id}.parsSpotSelector.removeIDs = str2num(get(handles.EditRemoveIDs,'String')); %#ok<ST2NM>



function EditDispHybs_Callback(hObject, eventdata, handles)
global CT
CT{handles.id}.parsSpotSelector.hybsToLoad = str2num(get(handles.EditDispHybs,'String')); %#ok<ST2NM>
LoadImage(handles);

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
global CT
% contents = cellstr(get(hObject,'String'));
% currentSelect = contents{get(hObject,'Value')};
if get(hObject,'Value') == 1
        CT{handles.id}.parsSpotSelector.dataType = 'fiducial';
else
        CT{handles.id}.parsSpotSelector.dataType = 'data';
end
LoadImage(handles);
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function EditFOV_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function EditRemoveIDs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function EditDispHybs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
