function varargout = EmbryoGUI(varargin)
% EMBRYOGUI MATLAB code for EmbryoGUI.fig
%      EMBRYOGUI, by itself, creates a new EMBRYOGUI or raises the existing
%      singleton*.
%
%      H = EMBRYOGUI returns the handle to a new EMBRYOGUI or the handle to
%      the existing singleton*.
%
%      EMBRYOGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EMBRYOGUI.M with the given input arguments.
%
%      EMBRYOGUI('Property','Value',...) creates a new EMBRYOGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EmbryoGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EmbryoGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EmbryoGUI

% Last Modified by GUIDE v2.5 22-Jan-2018 08:40:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EmbryoGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @EmbryoGUI_OutputFcn, ...
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


% --- Executes just before EmbryoGUI is made visible.
function EmbryoGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EmbryoGUI (see VARARGIN)


global EG
if isempty(EG) || isempty(EG{1})
    EG = cell(1,1);
else
    EG = [EG;cell(1,1)];
end
id = length(EG);
set(handles.instance,'String',['inst id',num2str(id)]);
handles.id = id;

% default parameters
EG{handles.id}.markDNA = false;
EG{handles.id}.markRegions = false;
EG{handles.id}.selectedChannels = [1,2]; % default
EG{handles.id}.mapRange = [100,450];

% Choose default command line output for EmbryoGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EmbryoGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = EmbryoGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function EditFile_Callback(hObject, eventdata, handles)
% hObject    handle to EditFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditFile as text
%        str2double(get(hObject,'String')) returns contents of EditFile as a double


% --- Executes during object creation, after setting all properties.
function EditFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ButtonLoad.
function ButtonLoad_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
global EG
fname = get(handles.EditFile,'String');
load(fname,'embDat'); 
% load channel names (if specified)
if isfield(embDat,'chnNames')
    set(handles.ListboxChannel,'String',char(embDat.chnNames));
else
    set(handles.ListboxChannel,'String',num2str( (1:size(embDat.im,3))' ) );
end
% load region segmentation data (if any)
blank = false;
if isfield(embDat,'embPolygonData')
    if ~isempty(embDat.embPolygonData)  
        regionAnnotations = embDat.embPolygonData;
        set(handles.ListboxRegions,'String',char(embDat.embPolygonData.Tag));
        % Convert embSegmentIDs to a matrix if necessary
        %    this allows each spot to be associated with multiple regions
        if ~isfield(embDat,'regionMatrix')
            nRegs = height(embDat.embPolygonData);
            nSpots = size(embDat.maps,3);
            regNames = embDat.embPolygonData.Tag;
            regionMatrix = false(nSpots,nRegs); 
            for r=1:nRegs
               regionMatrix(:,r) = strcmp(embDat.embSegmentIDs,regNames{r});
            end  
        else
            regionMatrix = embDat.regionMatrix;
        end
    else
        blank = true;
    end
else
    blank=true;
end
if blank
    regionMatrix = [];
    regionAnnotations = table();
end
% Export variables to GUI-global EG
EG{handles.id}.embDat = embDat;
EG{handles.id}.regionMatrix = regionMatrix;
EG{handles.id}.regionAnnotations = regionAnnotations;
% Update Image
ButtonUpdate_Callback(hObject, eventdata, handles);


function DisplayImage(hObject,eventdata,handles) %#ok<*INUSL>
% display the current image in imC
global EG
selected = EG{handles.id}.selectedChannels;
nChns = length(selected);
chnNames = cellstr(get(handles.ListboxChannel,'String')); 
cmap = [1 1 1; hsv(nChns-1)];
% Display image
figure(100+handles.id); clf; imagesc(EG{handles.id}.imC);
% Add a colorbar key
colormap(cmap);
h = colorbar;
set(h,'YTick',(1/nChns:1/(nChns):1)-.5/nChns,'YTickLabel', chnNames(selected));
% Optional: if requested, mark all DNA
if EG{handles.id}.markDNA
    hold on;
    xy=EG{handles.id}.embDat.spotXY;
    plot(xy(:,1),xy(:,2),'w+','MarkerSize',2);  
end
% Optional: if requested, draw boundaries around selected regions
if EG{handles.id}.markRegions
    hold on;
    AllRegions = EG{handles.id}.regionAnnotations;
    MarkerTags = cellstr(get(handles.ListboxRegions,'String')); 
    MarkerTags = MarkerTags(EG{handles.id}.selectedRegions);  
    polygonData = AllRegions(EG{handles.id}.selectedRegions,:);  
    nRegs = height(polygonData);
    cmap2 = hsv(nRegs); 
    for p=1:nRegs
        plot([polygonData.Positions{p}(:,1); polygonData.Positions{p}(1,1)],...
             [polygonData.Positions{p}(:,2); polygonData.Positions{p}(1,2)],...
            'color',cmap2(strcmp(polygonData.Tag{p},MarkerTags),:));
        xc = nanmean(polygonData.Positions{p}(:,1));
        yc = nanmean(polygonData.Positions{p}(:,2));
        text(xc,yc,polygonData.Tag{p},'color',cmap2(strcmp(polygonData.Tag{p},MarkerTags),:));
    end
end

% --- Executes on button press in ButtonUpdate.
function ButtonUpdate_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global EG
EG{handles.id}.selectedChannels = get(handles.ListboxChannel,'Value');
EG{handles.id}.selectedRegions = get(handles.ListboxRegions,'Value');
embDat = EG{handles.id}.embDat;
selected = EG{handles.id}.selectedChannels;
nChns = length(selected);
cmap = [1 1 1; hsv(nChns-1)];
imO = IncreaseContrast(embDat.im(:,:,selected),'low',.55,'high',.999995);
% imO = IncreaseContrast(embDat.im(:,:,selected),'low',.5,'high',.99999);
imC = Ncolor(imO,'colormap',cmap);
EG{handles.id}.imC = imC;
DisplayImage(hObject,eventdata,handles)





% --- Executes on button press in ButtonSelectRegions.
function ButtonSelectRegions_Callback(hObject, eventdata, handles)
% Add regions using the ZC segmentation GUI
global EG
regions = inputdlg('Enter region names as comma separated list',...
    'select names',1,{'region1, region2'});
if ~isempty(regions)
    regNames = strsplit( regions{1},', ');
    % call the ZC Segmentation GUI to draw regions manually
    [segIDs,regionAnnotations] = SegmentationGUI(EG{handles.id}.imC,...
        'Points', EG{handles.id}.embDat.spotXY,'MarkerTags',regNames,...
        'figHandle',100+handles.id);
    % update RegionMatrix truth table: regionMatrix
    nSpots = size(EG{handles.id}.embDat.spotXY,1);
    nRegs = height(regionAnnotations);
    newNames = regionAnnotations.Tag;
    addregionMatrix = false(nSpots,nRegs);
    for r = 1:nRegs
       addregionMatrix(:,r) = strcmp(segIDs,newNames{r});
    end
    % export new variables to the GUI global
    EG{handles.id}.regionMatrix = cat(2,EG{handles.id}.regionMatrix ,addregionMatrix);
    EG{handles.id}.regionAnnotations = cat(1,EG{handles.id}.regionAnnotations, regionAnnotations);
    % update Regions list
    currentRegions = cellstr(get(handles.ListboxRegions,'String')); 
    updatedRegions = cat(1,currentRegions,newNames);
    set(handles.ListboxRegions,'String',char(updatedRegions));
    % update GUI
    guidata(hObject, handles);
end

% --- Executes on button press in ButtonRemoveRegions.
function ButtonRemoveRegions_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonRemoveRegions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global EG;
currentRegions = cellstr(get(handles.ListboxRegions,'String')); 
selectedRegions = get(handles.ListboxRegions,'Value');
updatedRegions = currentRegions;
updatedRegions(selectedRegions) = []; % remove from cell array
EG{handles.id}.regionAnnotations(selectedRegions,:) =[]; % remove from table
EG{handles.id}.regionMatrix(:,selectedRegions,:) =[]; % remove from matrix
set(handles.ListboxRegions,'String',char(updatedRegions)); % update GUI
guidata(hObject, handles);

% --- Executes on button press in ButtonCompareMaps.
function ButtonCompareMaps_Callback(hObject, eventdata, handles)
% Plots distance maps for the selected regions
global EG
sel = get(handles.ListboxRegions,'Value');
currentRegions = cellstr(get(handles.ListboxRegions,'String'));
selLabels = currentRegions(sel);
figure(200+handles.id); clf; 
for r=1:length(sel)
    subplot(1,length(sel),r); 
    reg = EG{handles.id}.regionMatrix(:,sel(r)); 
    im = nanmedian(EG{handles.id}.embDat.maps(:,:,reg),3);
    imagesc(im); caxis(EG{handles.id}.mapRange); colorbar;
    title(selLabels{r});
end
colormap(flipud(parula));

% --- Executes on button press in ToggleMarkRegions.
function ToggleMarkRegions_Callback(hObject, eventdata, handles)
global EG
toggled = get(handles.ToggleMarkRegions,'Value');
if toggled
    EG{handles.id}.markRegions = true;
else
    EG{handles.id}.markRegions = false;
end

% --- Executes on button press in ToggleMarkDNA.
function ToggleMarkDNA_Callback(hObject, eventdata, handles)
global EG
set(handles.ListboxRegions,'Max',50);
toggled = get(handles.ToggleMarkDNA,'Value');
if toggled
    EG{handles.id}.markDNA = true;
else
    EG{handles.id}.markDNA = false;
end

% --------------------------------------------------------------------
function MenuSaveData_Callback(hObject, eventdata, handles)
% hObject    handle to MenuSaveData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
saveName = get(handles.EditFile,'String');
SaveData(hObject,eventdata,handles,saveName)

% --------------------------------------------------------------------
function MenuSaveAs_Callback(hObject, eventdata, handles)
% hObject    handle to MenuSaveAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
saveName = get(handles.EditFile,'String');
% [folder,file,type] = fileparts(saveName);
[filename,pathname,zeroIfQuit] = uiputfile('*.mat','Save Data As',saveName);
if zeroIfQuit~=0
    SaveData(hObject,eventdata,handles,[pathname, filename])
end
% --------------------------------------------------------------------
function SaveData(hObject,eventdata,handles,saveName)
global EG
currentRegions = cellstr(get(handles.ListboxRegions,'String')); % these should already be in the updated regionAnnotation table 
embDat = EG{handles.id}.embDat;
embDat.embPolygonData = EG{handles.id}.regionAnnotations; % update table of recording names and boundaries of selected regions 
embDat.regionMatrix = EG{handles.id}.regionMatrix; % Matrix of which regions each spot is a member of 
disp(['writing ',saveName]);
save(saveName,'embDat');


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global EG;

choice = questdlg('Close embryo GUI?','Confirm','Yes','No','Yes');
if strcmp(choice,'Yes')
    % Hint: delete(hObject) closes the figure
    try
        figure(100+handles.id); close;
        if length(EG) == handles.id
            EG(handles.id) = [];
        else
            EG{handles.id} = 'empty';
        end
    catch
    end
    delete(hObject);
end

%============================
%% Auto Generated
%============================

% --------------------------------------------------------------------
function MenuFile_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in ListboxChannel.
function ListboxChannel_Callback(hObject, eventdata, handles)
% hObject    handle to ListboxChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ListboxChannel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ListboxChannel


% --- Executes during object creation, after setting all properties.
function ListboxChannel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ListboxChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in ListboxRegions.
function ListboxRegions_Callback(hObject, eventdata, handles)
% hObject    handle to ListboxRegions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ListboxRegions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ListboxRegions


% --- Executes during object creation, after setting all properties.
function ListboxRegions_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% To delete

% --- Executes on button press in CheckDNA.
function CheckDNA_Callback(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to CheckDNA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckDNA
