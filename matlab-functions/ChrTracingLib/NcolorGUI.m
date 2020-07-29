function varargout = NcolorGUI(varargin)
% NCOLORGUI MATLAB code for NcolorGUI.fig
%      NCOLORGUI, by itself, creates a new NCOLORGUI or raises the existing
%      singleton*.
%
%      H = NCOLORGUI returns the handle to a new NCOLORGUI or the handle to
%      the existing singleton*.
%
%      NCOLORGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NCOLORGUI.M with the given input arguments.
%
%      NCOLORGUI('Property','Value',...) creates a new NCOLORGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NcolorGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NcolorGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NcolorGUI

% Last Modified by GUIDE v2.5 20-Jun-2018 19:22:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NcolorGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @NcolorGUI_OutputFcn, ...
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


% --- Executes just before NcolorGUI is made visible.
function NcolorGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NcolorGUI (see VARARGIN)
if length(varargin) < 1
    error('No input.  Please pass an HxWxN matrix to NcolorGUI');
end

% setup global GUI data struct
global NG
if isempty(NG) 
     NG = cell(1,1);
     NG{1} = struct();
end
 id = 1; % fixed instance
% % multiple instances
% if isempty(NG) || isempty(NG{1})
%     NG = cell(1,1);
% else
%     NG = [NG;cell(1,1)];
% end
% id = length(NG);

% load defaults into current values
im = varargin{1};
defaults = cell(0,3);
defaults(end+1,:) = {'names','cell',{}};
defaults(end+1,:) = {'colormap','colormap','hsv'};  % No option for type 'table'
defaults(end+1,:) = {'figHandle','freeType',100}; %  
defaults(end+1,:) = {'instanceID','freeType',id};
defaults(end+1,:) = {'contrastMax','freeType',.999};
defaults(end+1,:) = {'contrastMin','freeType',.5};
defaults(end+1,:) = {'restart','boolean',true};
pars = ParseVariableArguments(varargin(2:end),defaults,mfilename);

id = pars.instanceID;
if ~ishandle(pars.figHandle)
    pars.figHandle = pars.figHandle + pars.instanceID;
end
set(handles.instance,'String',['inst id',num2str(id)]);
handles.id = id;

[ymax,xmax,nChns] = size(im);
NG{handles.id}.im = im;
NG{handles.id}.pars  = pars;
% set defaults for parameters
NG{handles.id}.fov = [1,xmax,1,ymax];
% if ~isfield(NG{handles.id},'fov') || pars.restart
%     NG{handles.id}.fov = [1,xmax,1,ymax];
% end
if ~isfield(NG{handles.id},'selectedChannels')  || pars.restart
    NG{handles.id}.selectedChannels = 1;
end
if ~isfield(NG{handles.id},'tuneChannels')  || pars.restart
    NG{handles.id}.tuneChannels = 1;
end
% Setup default contrast options
if isempty(pars.contrastMin)
    pars.contrastMin = 0.5;
end
if isempty(pars.contrastMax)
    pars.contrastMax = .9999;
end
if length(pars.contrastMin) < nChns
    NG{handles.id}.contrastMin = pars.contrastMin.*ones(nChns,1);
else
    NG{handles.id}.contrastMin = pars.contrastMin;
end
if length(pars.contrastMax) < nChns
    NG{handles.id}.contrastMax = pars.contrastMax.*ones(nChns,1);
else
    NG{handles.id}.contrastMax = pars.contrastMax;
end
% setup colormap
NG{handles.id}.cmap = GetColorMap(pars.colormap,nChns);


% setup default output
NG{handles.id}.imO = 0;
NG{handles.id}.imC = 0;
if isempty(pars.figHandle)
    pars.figHandle = 100 + handles.id;
end
handles.figureH = figure(pars.figHandle);
handles.figureH.Name = 'NcolorGUI';

% set up sliders
n=1;
set(handles.MinSlider,'Value',NG{handles.id}.contrastMin(n));
set(handles.MaxSlider,'Value',NG{handles.id}.contrastMax(n));
set(handles.EditMin,'String',num2str(NG{handles.id}.contrastMin(n)));
set(handles.EditMax,'String',num2str(NG{handles.id}.contrastMax(n)));
set(handles.ListboxDisplay,'Max',50);
set(handles.ListboxTune,'Max',1);
set(handles.ListboxDisplay,'Value',NG{handles.id}.selectedChannels);

% setup names
if ~isempty(pars.names)
  set(handles.ListboxDisplay,'String',char(pars.names));
  set(handles.ListboxTune,'String',char(pars.names));
else
   set(handles.ListboxDisplay,'String',num2str( (1:nChns)' ) ); 
   set(handles.ListboxTune,'String',num2str( (1:nChns)' ) ); 
   NG{handles.id}.pars.names = cellstr(num2str( (1:nChns)'));
end

% Choose default command line output for NcolorGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
NG{handles.id}.NcolorGUI_handle = handles.figure1;


ButtonUpdate_Callback(hObject, eventdata, handles)

% UIWAIT makes NcolorGUI wait for user response (see UIRESUME)
if nargout > 0 % brilliant, only freezes if requested  output
    uiwait(handles.figure1);
end


% --- Outputs from this function are returned to the command line.
function varargout = NcolorGUI_OutputFcn(hObject, eventdata, handles)  %#ok<*INUSL>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global NG
% Get default command line output from handles structure
varargout{1} = NG{handles.id}.imO;
varargout{2} = NG{handles.id}.imC;
varargout{3} = [NG{handles.id}.contrastMin,NG{handles.id}.contrastMax];
% delete(handles.figure1);  % Close window




% --- Executes on button press in ButtonUpdate.
function ButtonUpdate_Callback(hObject, eventdata, handles)
% read in the selected options and update
global NG
im = NG{handles.id}.im;
selected = get(handles.ListboxDisplay,'Value');
NG{handles.id}.selectedChannels = selected;

% get FOV (allows hold zoom to work)
xmin = NG{handles.id}.fov(1);  
xmax = NG{handles.id}.fov(2); 
ymin = NG{handles.id}.fov(3); 
ymax = NG{handles.id}.fov(4); 
im = im(ymin:ymax,xmin:xmax,:);
% apply selected contrast
[h,w,totChns] = size(im);
contrastMin = NG{handles.id}.contrastMin;% (selected);
contrastMax = NG{handles.id}.contrastMax;% (selected);
imO = zeros(h,w,totChns,'uint16');
for n=1:totChns
    if any(n==selected)
        imO(:,:,n) = IncreaseContrast(im(:,:,n),...
                        'low',contrastMin(n),...
                        'high',contrastMax(n));
    end
end
% apply selected colormap
cmap = NG{handles.id}.cmap;
imC = Ncolor(imO,'colormap',cmap);
NG{handles.id}.imO = imO;
NG{handles.id}.imC = imC;
DisplayImage(hObject,eventdata,handles);

% --- Main display function
% display the current image in imC
function DisplayImage(hObject,eventdata,handles)
global NG
% Display image
%figure(100+handles.id); clf;
figure(handles.figureH); clf;
imagesc(NG{handles.id}.imC);
% Add color bar key
%  Remove inactive names from colorbar for clarity
selected = NG{handles.id}.selectedChannels;
nChns = size( NG{handles.id}.im,3);
isSelected = false(nChns,1);
isSelected(selected) = true;
chnNames = cellstr(get(handles.ListboxDisplay,'String')); 
chnNames(~isSelected) = {' '};
colormap(NG{handles.id}.cmap);
h = colorbar;
set(h,'YTick',(1/nChns:1/(nChns):1)-.5/nChns,'YTickLabel', chnNames);



% --- Executes on slider movement.
function MaxSlider_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
global NG
n = get(handles.ListboxTune,'Value');
NG{handles.id}.contrastMax(n) = get(handles.MaxSlider,'Value');
set(handles.EditMax,'String',num2str(NG{handles.id}.contrastMax(n)));
guidata(hObject, handles);
ButtonUpdate_Callback(hObject, eventdata, handles);

% --- Executes on slider movement.
function MinSlider_Callback(hObject, eventdata, handles)
global NG
n = get(handles.ListboxTune,'Value');
NG{handles.id}.contrastMin(n) = get(handles.MinSlider,'Value');
set(handles.EditMin,'String',num2str(NG{handles.id}.contrastMin(n)));
guidata(hObject, handles);
ButtonUpdate_Callback(hObject, eventdata, handles);

% -- Executes on edit of text
function EditMin_Callback(hObject, eventdata, handles)
global NG
n = get(handles.ListboxTune,'Value');
NG{handles.id}.contrastMin(n) = str2double(get(handles.EditMin,'String'));
set(handles.MinSlider,'Value',NG{handles.id}.contrastMin(n));
guidata(hObject, handles);
ButtonUpdate_Callback(hObject, eventdata, handles);

% -- Executes on edit of text
function EditMax_Callback(hObject, eventdata, handles)
global NG
n = get(handles.ListboxTune,'Value');
NG{handles.id}.contrastMax(n) = str2double(get(handles.EditMax,'String'));
set(handles.MaxSlider,'Value',NG{handles.id}.contrastMax(n));
guidata(hObject, handles);
ButtonUpdate_Callback(hObject, eventdata, handles);



% --- Executes on button press in ButtonExportAndClose.
function ButtonExportAndClose_Callback(hObject, eventdata, handles)
if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    uiresume(handles.figure1);  % Triggers SegmentationGUI_OutputFcn
else
    delete(handles.figure1);
end

% --- Executes on selection change in ListboxDisplay.
function ListboxDisplay_Callback(hObject, eventdata, handles)
ButtonUpdate_Callback(hObject, eventdata, handles);

% --- Executes on selection change in ListboxTune.
function ListboxTune_Callback(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to ListboxTune (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global NG
n = get(handles.ListboxTune,'Value');
contrastMax = NG{handles.id}.contrastMax(n);
contrastMin = NG{handles.id}.contrastMin(n);
set(handles.MinSlider,'Value',contrastMin);
set(handles.MaxSlider,'Value',contrastMax);
set(handles.EditMin,'String',num2str(contrastMin));
set(handles.EditMax,'String',num2str(contrastMax));



% Hints: contents = cellstr(get(hObject,'String')) returns ListboxTune contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ListboxTune




%% Menu options
% --------------------------------------------------------------------
function MenuFile_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MenuEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MenuEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MenuColormap_Callback(hObject, eventdata, handles)
% hObject    handle to MenuColormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global NG
cmap = NG{handles.id}.cmap;
nChns = size(cmap,1);
choices = inputdlg({'values'},'select values',nChns,{num2str(cmap)});
if ~isempty(choices)
    cmapOut = str2num(choices{1}); %#ok<ST2NM>
    NG{handles.id}.cmap = cmapOut;
end

% --------------------------------------------------------------------
function MenuSavePNG_Callback(hObject, eventdata, handles)
% hObject    handle to MenuSavePNG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('sorry not available yet!');

% --------------------------------------------------------------------
function MenuSaveMat_Callback(hObject, eventdata, handles)
% hObject    handle to MenuSaveMat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('sorry not available yet!');










%% GUI display setup
%%=========================================================================
% --- Executes during object creation, after setting all properties.
function EditMax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function EditMin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function MinSlider_CreateFcn(hObject, eventdata, handles)
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes during object creation, after setting all properties.
function MaxSlider_CreateFcn(hObject, eventdata, handles)
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function ListboxDisplay_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function ListboxTune_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: delete(hObject) closes the figure
global NG;
try
    figure(100+handles.id); close;
    if length(NG) == handles.id
        NG(handles.id) = [];
    else
        NG{handles.id} = 'empty';
    end
catch
end
delete(hObject);


% --- Executes on button press in ToggleZoom.
function ToggleZoom_Callback(hObject, eventdata, handles)
% hObject    handle to ToggleZoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global NG
im = NG{handles.id}.im;
isOn = get(handles.ToggleZoom,'Value');
if isOn
    figure(100+handles.id);
    xlims = round(get(gca,'Xlim'));
    ylims = round(get(gca,'Ylim'));
    xmin = xlims(1);
    xmax = xlims(2);
    ymin = ylims(1);
    ymax = ylims(2);
else
   xmin = 1;
   ymin = 1;
   [ymax,xmax,~] = size(im); 
   disp('reset');
end
NG{handles.id}.fov = [xmin,xmax,ymin,ymax];


% --------------------------------------------------------------------
function ToolZoomIn_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ToolZoomIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ToolZoomOut_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ToolZoomOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
