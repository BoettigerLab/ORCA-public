function varargout = ChrDataViewer(varargin)
% CHRDATAVIEWER MATLAB code for ChrDataViewer.fig
%      CHRDATAVIEWER, by itself, creates a new CHRDATAVIEWER or raises the existing
%      singleton*.
%
%      H = CHRDATAVIEWER returns the handle to a new CHRDATAVIEWER or the handle to
%      the existing singleton*.
%
%      CHRDATAVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHRDATAVIEWER.M with the given input arguments.
%
%      CHRDATAVIEWER('Property','Value',...) creates a new CHRDATAVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ChrDataViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ChrDataViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ChrDataViewer

% Last Modified by GUIDE v2.5 09-Nov-2017 15:15:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ChrDataViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @ChrDataViewer_OutputFcn, ...
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


% --- Executes just before ChrDataViewer is made visible.
function ChrDataViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ChrDataViewer (see VARARGIN)

global DV
if isempty(DV)
    DV = cell(1,1);
else
    DV = [DV;cell(1,1)];
end
id = length(DV);
set(handles.DVinstance,'String',['inst id',num2str(id)]);
handles.id = id;

DV{handles.id}.dataFolder = '';
DV{handles.id}.currentCell = 1;

% Choose default command line output for ChrDataViewer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ChrDataViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ChrDataViewer_OutputFcn(hObject, eventdata, handles)  %#ok<*INUSL>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



%=========================================================================%
% --- Executes on button press in ButtonLoadData.
function ButtonLoadData_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonLoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DV

showXY = get(handles.CheckDisplayXY,'Value');
showXZ = get(handles.CheckDisplayXZ,'Value');
loadDat = get(handles.CheckLoadData,'Value');
loadFid = get(handles.CheckLoadFid,'Value');
loadFits = get(handles.CheckLoadFits,'Value');
showLevels = get(handles.CheckDisplayLevels,'Value');
showMap = get(handles.CheckDisplayMap,'Value');
showPolymer = get(handles.CheckDisplayPolymer,'Value');
showOverlay = get(handles.CheckDisplayOverlays,'Value');

if isempty(DV{handles.id}.dataFolder)
  DV{handles.id}.dataFolder = get(handles.EditDataFolder,'String');
  if isempty(DV{handles.id}.dataFolder)
      error('Please specify a data folder containing i4d files');
  end
end
% dataFolder = 'J:\2017-09-19_IMR90_chr21_CT2out\'; 
% fovFolder = [dataFolder,'fov001\'];
% fovFolder = 'J:\2017-09-19_IMR90_chr21_CT2out\fov001\'

fovFolder = DV{handles.id}.dataFolder; 
c = DV{handles.id}.currentCell;
% c = 1; % 8 which cell to look at; 

i4ds = strcat(fovFolder,cellstr(ls([fovFolder,'*AlignedData.i4d'])));
DV{handles.id}.numCells = length(i4ds); 

if loadDat
    datMat = ReadImage4D(i4ds{c},'showPlots',false);
    nIms = size(datMat,4);
end
if loadFits
    fits = readtable( regexprep(i4ds{c},'AlignedData.i4d','fits.csv') );
else
    fits = [];
end
if loadFid
    fidMat = ReadImage4D(regexprep(i4ds{c},'AlignedData','AlignedFid'),'showPlots',false);
    nIms = size(datMat,4);
end
tileLabels = cellstr( num2str((1:nIms)') );

if showOverlay
    figure(103); clf;
    if loadFid && loadDat
        subplot(1,2,1); Ncolor((2/nIms)*IncreaseContrast(squeeze(max(fidMat(:,:,:,:),[],3))));
        subplot(1,2,2); Ncolor((2/nIms)*IncreaseContrast(squeeze(max(datMat,[],3))));
    elseif loadFid
        Ncolor((2/nIms)*IncreaseContrast(squeeze(max(fidMat(:,:,:,:),[],3))));
    elseif loadDat
        Ncolor((2/nIms)*IncreaseContrast(squeeze(max(datMat,[],3))));
    end
end

figure(101); clf;
if showXY
    if loadFid && loadDat
        subplot(2,1,1); PlotProjection4D(fidMat,'fits',[],'projection','xy','tileLabels',tileLabels); axis image;
        subplot(2,1,2); PlotProjection4D(datMat,'fits',fits,'projection','xy','tileLabels',tileLabels); axis image;
    elseif loadFid
        PlotProjection4D(fidMat,'fits',[],'projection','xy','tileLabels',tileLabels); axis image;
    elseif loadDat
        PlotProjection4D(datMat,'fits',fits,'projection','xy','tileLabels',tileLabels); axis image;
    end
end

figure(102); clf;
if showXZ
    if loadFid && loadDat
        subplot(2,1,1); PlotProjection4D(fidMat,'fits',[],'projection','xz','tileLabels',tileLabels); axis image;
        subplot(2,1,2); PlotProjection4D(datMat,'fits',fits,'projection','xz','tileLabels',tileLabels); axis image;
    elseif loadFid
        PlotProjection4D(fidMat,'fits',[],'projection','xz','tileLabels',tileLabels); axis image;
    elseif loadDat
        PlotProjection4D(datMat,'fits',fits,'projection','xz','tileLabels',tileLabels); axis image;
    end
end

if showLevels
    pixMax = zeros(nIms,1);
    for i=1:nIms
        datIm = datMat(:,:,:,i);
        pixMax(i) = max(datIm(:));
    end
    figure(105); clf; bar(pixMax);
end

if showMap
    
end


if showPolymer
    
end

%========================================================================%

function EditDataFolder_Callback(hObject, eventdata, handles) %#ok<*DEFNU,*INUSD>
% hObject    handle to EditDataFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditDataFolder as text
%        str2double(get(hObject,'String')) returns contents of EditDataFolder as a double


% --- Executes during object creation, after setting all properties.
function EditDataFolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditDataFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditCellNum_Callback(hObject, eventdata, handles)
% hObject    handle to EditCellNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global DV
DV{handles.id}.currentCell = str2double(get(handles.EditCellNum,'String'));
% Hints: get(hObject,'String') returns contents of EditCellNum as text
%        str2double(get(hObject,'String')) returns contents of EditCellNum as a double


% --- Executes during object creation, after setting all properties.
function EditCellNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditCellNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function SliderCellNum_Callback(hObject, eventdata, handles)
% hObject    handle to SliderCellNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function SliderCellNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SliderCellNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in CheckDisplayXZ.
function CheckDisplayXZ_Callback(hObject, eventdata, handles)
% hObject    handle to CheckDisplayXZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckDisplayXZ


% --- Executes on button press in CheckDisplayXY.
function CheckDisplayXY_Callback(hObject, eventdata, handles)
% hObject    handle to CheckDisplayXY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckDisplayXY


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in CheckLoadFid.
function CheckLoadFid_Callback(hObject, eventdata, handles)
% hObject    handle to CheckLoadFid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckLoadFid


% --- Executes on button press in CheckLoadData.
function CheckLoadData_Callback(hObject, eventdata, handles)
% hObject    handle to CheckLoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckLoadData


% --- Executes on button press in CheckLoadFits.
function CheckLoadFits_Callback(hObject, eventdata, handles)
% hObject    handle to CheckLoadFits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckLoadFits


% --- Executes on button press in ButtonKeep.
function ButtonKeep_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonKeep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ButtonReject.
function ButtonReject_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonReject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ButtonRefitAll.
function ButtonRefitAll_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonRefitAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ButtonRefitData.
function ButtonRefitData_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonRefitData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in CheckDisplayLevels.
function CheckDisplayLevels_Callback(hObject, eventdata, handles)
% hObject    handle to CheckDisplayLevels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckDisplayLevels


% --- Executes on button press in CheckDisplayMap.
function CheckDisplayMap_Callback(hObject, eventdata, handles)
% hObject    handle to CheckDisplayMap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckDisplayMap


% --- Executes on button press in CheckDisplayPolymer.
function CheckDisplayPolymer_Callback(hObject, eventdata, handles)
% hObject    handle to CheckDisplayPolymer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckDisplayPolymer

