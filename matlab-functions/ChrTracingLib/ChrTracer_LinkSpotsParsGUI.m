function varargout = ChrTracer_LinkSpotsParsGUI(varargin)
% CHRTRACER_LINKSPOTSPARSGUI MATLAB code for ChrTracer_LinkSpotsParsGUI.fig
%      CHRTRACER_LINKSPOTSPARSGUI, by itself, creates a new CHRTRACER_LINKSPOTSPARSGUI or raises the existing
%      singleton*.
%
%      H = CHRTRACER_LINKSPOTSPARSGUI returns the handle to a new CHRTRACER_LINKSPOTSPARSGUI or the handle to
%      the existing singleton*.
%
%      CHRTRACER_LINKSPOTSPARSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHRTRACER_LINKSPOTSPARSGUI.M with the given input arguments.
%
%      CHRTRACER_LINKSPOTSPARSGUI('Property','Value',...) creates a new CHRTRACER_LINKSPOTSPARSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ChrTracer_LinkSpotsParsGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ChrTracer_LinkSpotsParsGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ChrTracer_LinkSpotsParsGUI

% Last Modified by GUIDE v2.5 14-Mar-2017 15:28:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ChrTracer_LinkSpotsParsGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ChrTracer_LinkSpotsParsGUI_OutputFcn, ...
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


% --- Executes just before ChrTracer_LinkSpotsParsGUI is made visible.
function ChrTracer_LinkSpotsParsGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ChrTracer_LinkSpotsParsGUI (see VARARGIN)

% Choose default command line output for ChrTracer_LinkSpotsParsGUI
handles.output = hObject;

if ~isempty(varargin)
    id = varargin{1}; 
    handles.id = id;
else
    handles.id = 1;
end

% Choose default command line output for ChrTracer_PSFparsGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
SetupDefaults(hObject,eventdata,handles)

% UIWAIT makes ChrTracer_LinkSpotsParsGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ChrTracer_LinkSpotsParsGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Internal function
function SetupDefaults(hObject,eventdata,handles)
global CT;
passedStruct = CT{handles.id}.parsLink;
passedTags = fieldnames(passedStruct);
localTags = fieldnames(handles); 
for i=1:length(localTags)
    k = StringFind(passedTags,localTags{i},'exactly',true);
    if ~isempty(k)
        value = passedStruct.(passedTags{k});
        if islogical(value)
            set(handles.(localTags{i}),'Value',value);
        else
            set(handles.(localTags{i}),'String',num2str(value));
        end
    end
end
CT{handles.id}.parsLink = passedStruct;
guidata(hObject, handles);


% --- Executes on button press in ButtonUpdate.
function ButtonUpdate_Callback(hObject, eventdata, handles)  %#ok<*INUSL,*DEFNU>
% hObject    handle to ButtonUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CT;
passedStruct = CT{handles.id}.parsLink;
passedTags = fieldnames(passedStruct);
localTags = fieldnames(handles); 
for i=1:length(localTags)
    k = StringFind(passedTags,localTags{i},'exactly',true);
    if ~isempty(k)
        if islogical(passedStruct.(passedTags{k}))
            value = logical(get(handles.(localTags{i}),'Value'));
        elseif ischar(passedStruct.(passedTags{k}))
            value = get(handles.(localTags{i}),'String');
        elseif isfloat(passedStruct.(passedTags{k}))
            value = get(handles.(localTags{i}),'String');
        end
        if islogical(value)
            passedStruct.(passedTags{k}) = value;
        elseif ischar(passedStruct.(passedTags{k}))
            passedStruct.(passedTags{k}) = value;
        elseif isfloat(passedStruct.(passedTags{k}))
            passedStruct.(passedTags{k}) = str2double(value);
        end
    end
end
CT{handles.id}.parsLink = passedStruct;
disp('updated parameters');

% --- Executes on button press in ButtonCancel.
function ButtonCancel_Callback(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to ButtonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(ChrTracer_LinkSpotsParsGUI);




function maxDistFromStart_Callback(hObject, eventdata, handles)
% hObject    handle to maxDistFromStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxDistFromStart as text
%        str2double(get(hObject,'String')) returns contents of maxDistFromStart as a double


% --- Executes during object creation, after setting all properties.
function maxDistFromStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxDistFromStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxJump_Callback(hObject, eventdata, handles)
% hObject    handle to maxJump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxJump as text
%        str2double(get(hObject,'String')) returns contents of maxJump as a double


% --- Executes during object creation, after setting all properties.
function maxJump_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxJump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function renderTubeRadius_Callback(hObject, eventdata, handles)
% hObject    handle to renderTubeRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of renderTubeRadius as text
%        str2double(get(hObject,'String')) returns contents of renderTubeRadius as a double


% --- Executes during object creation, after setting all properties.
function renderTubeRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to renderTubeRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function renderBallRadius_Callback(hObject, eventdata, handles)
% hObject    handle to renderBallRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of renderBallRadius as text
%        str2double(get(hObject,'String')) returns contents of renderBallRadius as a double


% --- Executes during object creation, after setting all properties.
function renderBallRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to renderBallRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function renderInterpMethod_Callback(hObject, eventdata, handles)
% hObject    handle to renderInterpMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of renderInterpMethod as text
%        str2double(get(hObject,'String')) returns contents of renderInterpMethod as a double


% --- Executes during object creation, after setting all properties.
function renderInterpMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to renderInterpMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
