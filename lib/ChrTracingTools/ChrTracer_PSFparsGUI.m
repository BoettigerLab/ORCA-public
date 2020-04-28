function varargout = ChrTracer_PSFparsGUI(varargin)
% CHRTRACER_PSFPARSGUI MATLAB code for ChrTracer_PSFparsGUI.fig
%      CHRTRACER_PSFPARSGUI, by itself, creates a new CHRTRACER_PSFPARSGUI or raises the existing
%      singleton*.
%
%      H = CHRTRACER_PSFPARSGUI returns the handle to a new CHRTRACER_PSFPARSGUI or the handle to
%      the existing singleton*.
%
%      CHRTRACER_PSFPARSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHRTRACER_PSFPARSGUI.M with the given input arguments.
%
%      CHRTRACER_PSFPARSGUI('Property','Value',...) creates a new CHRTRACER_PSFPARSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ChrTracer_PSFparsGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ChrTracer_PSFparsGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% This is a simple parameters GUI
% There are only two functions of note:
% SetupDefaults
%    This function scans fields from the input global (in this case CT) and
%    finds matching fields in the GUI handles.  
%    The GUI fields are updated to match the current values in the input
%    global.  Booleans are checkboxes. All other variables are numeric.
%    Strings could be allowed by sorting out all the fields with the word
%    "name" or "string" in the fieldname for example. 
% ButtonOK_Callback
%   This function updates values in the input global variable with the
%   values entered in the GUI. 
% 

% Edit the above text to modify the response to help ChrTracer_PSFparsGUI

% Last Modified by GUIDE v2.5 15-Mar-2017 18:34:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ChrTracer_PSFparsGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ChrTracer_PSFparsGUI_OutputFcn, ...
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


% --- Executes just before ChrTracer_PSFparsGUI is made visible.
function ChrTracer_PSFparsGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ChrTracer_PSFparsGUI (see VARARGIN)

% "Opening Func" actually called at opening AND closing!
% when called at closing there is no varargin input...  
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
if ~isfield(handles,'closing')
    SetupDefaults(hObject,eventdata,handles)
end

% UIWAIT makes ChrTracer_PSFparsGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ChrTracer_PSFparsGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Internal function
function SetupDefaults(hObject,eventdata,handles)
global CT; 
passedStruct = CT{handles.id}.parsFit;
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
CT{handles.id}.parsFit = passedStruct;
guidata(hObject, handles);

% --- Executes on button press in ButtonOK.
function ButtonOK_Callback(hObject, eventdata, handles) %#ok<*INUSL,*DEFNU>
% hObject    handle to ButtonOK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CT;
passedStruct = CT{handles.id}.parsFit;
passedTags = fieldnames(passedStruct);
localTags = fieldnames(handles); 
for i=1:length(localTags)
    k = StringFind(passedTags,localTags{i},'exactly',true);
    if ~isempty(k)
        if islogical(passedStruct.(passedTags{k}))
            value = logical(get(handles.(localTags{i}),'Value'));
        else
            value = get(handles.(localTags{i}),'String');
        end
        if islogical(value)
            passedStruct.(passedTags{k}) = value;
        else
            passedStruct.(passedTags{k}) = str2double(value);
        end
    end
end
CT{handles.id}.parsFit = passedStruct;


% --- Executes on button press in ButtonCancel.
function ButtonCancel_Callback(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to ButtonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.closing = true;
guidata(hObject, handles);
close(ChrTracer_PSFparsGUI);


function datMinPeakHeight_Callback(hObject, eventdata, handles)
% hObject    handle to datMinPeakHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of datMinPeakHeight as text
%        str2double(get(hObject,'String')) returns contents of datMinPeakHeight as a double


% --- Executes during object creation, after setting all properties.
function datMinPeakHeight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datMinPeakHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function datPeakBlur_Callback(hObject, eventdata, handles)
% hObject    handle to datPeakBlur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of datPeakBlur as text
%        str2double(get(hObject,'String')) returns contents of datPeakBlur as a double


% --- Executes during object creation, after setting all properties.
function datPeakBlur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datPeakBlur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function datCameraBackground_Callback(hObject, eventdata, handles)
% hObject    handle to datCameraBackground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of datCameraBackground as text
%        str2double(get(hObject,'String')) returns contents of datCameraBackground as a double


% --- Executes during object creation, after setting all properties.
function datCameraBackground_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datCameraBackground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in datTroubleshoot.
function datTroubleshoot_Callback(hObject, eventdata, handles)
% hObject    handle to datTroubleshoot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of datTroubleshoot



function datMaxFitWidth_Callback(hObject, eventdata, handles)
% hObject    handle to datMaxFitWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of datMaxFitWidth as text
%        str2double(get(hObject,'String')) returns contents of datMaxFitWidth as a double


% --- Executes during object creation, after setting all properties.
function datMaxFitWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datMaxFitWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fidMinPeakHeight_Callback(hObject, eventdata, handles)
% hObject    handle to fidMinPeakHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fidMinPeakHeight as text
%        str2double(get(hObject,'String')) returns contents of fidMinPeakHeight as a double


% --- Executes during object creation, after setting all properties.
function fidMinPeakHeight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fidMinPeakHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fidPeakBlur_Callback(hObject, eventdata, handles)
% hObject    handle to fidPeakBlur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fidPeakBlur as text
%        str2double(get(hObject,'String')) returns contents of fidPeakBlur as a double


% --- Executes during object creation, after setting all properties.
function fidPeakBlur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fidPeakBlur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fidCameraBackground_Callback(hObject, eventdata, handles)
% hObject    handle to fidCameraBackground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fidCameraBackground as text
%        str2double(get(hObject,'String')) returns contents of fidCameraBackground as a double


% --- Executes during object creation, after setting all properties.
function fidCameraBackground_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fidCameraBackground (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fidTroubleshoot.
function fidTroubleshoot_Callback(hObject, eventdata, handles)
% hObject    handle to fidTroubleshoot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fidTroubleshoot



function fidMaxFitWidth_Callback(hObject, eventdata, handles)
% hObject    handle to fidMaxFitWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fidMaxFitWidth as text
%        str2double(get(hObject,'String')) returns contents of fidMaxFitWidth as a double


% --- Executes during object creation, after setting all properties.
function fidMaxFitWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fidMaxFitWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function datMinSep_Callback(hObject, eventdata, handles)
% hObject    handle to datMinSep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of datMinSep as text
%        str2double(get(hObject,'String')) returns contents of datMinSep as a double


% --- Executes during object creation, after setting all properties.
function datMinSep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datMinSep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fidMinSep_Callback(hObject, eventdata, handles)
% hObject    handle to fidMinSep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fidMinSep as text
%        str2double(get(hObject,'String')) returns contents of fidMinSep as a double


% --- Executes during object creation, after setting all properties.
function fidMinSep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fidMinSep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function datKeepBrightest_Callback(hObject, eventdata, handles)
% hObject    handle to datKeepBrightest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of datKeepBrightest as text
%        str2double(get(hObject,'String')) returns contents of datKeepBrightest as a double


% --- Executes during object creation, after setting all properties.
function datKeepBrightest_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datKeepBrightest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function datRelativeHeight_Callback(hObject, eventdata, handles)
% hObject    handle to datRelativeHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of datRelativeHeight as text
%        str2double(get(hObject,'String')) returns contents of datRelativeHeight as a double


% --- Executes during object creation, after setting all properties.
function datRelativeHeight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datRelativeHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function datMinHBratio_Callback(hObject, eventdata, handles)
% hObject    handle to datMinHBratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of datMinHBratio as text
%        str2double(get(hObject,'String')) returns contents of datMinHBratio as a double


% --- Executes during object creation, after setting all properties.
function datMinHBratio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datMinHBratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function datMinAHratio_Callback(hObject, eventdata, handles)
% hObject    handle to datMinAHratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of datMinAHratio as text
%        str2double(get(hObject,'String')) returns contents of datMinAHratio as a double


% --- Executes during object creation, after setting all properties.
function datMinAHratio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datMinAHratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function datMaxUncert_Callback(hObject, eventdata, handles)
% hObject    handle to datMaxUncert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of datMaxUncert as text
%        str2double(get(hObject,'String')) returns contents of datMaxUncert as a double


% --- Executes during object creation, after setting all properties.
function datMaxUncert_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datMaxUncert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fidKeepBrightest_Callback(hObject, eventdata, handles)
% hObject    handle to fidKeepBrightest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fidKeepBrightest as text
%        str2double(get(hObject,'String')) returns contents of fidKeepBrightest as a double


% --- Executes during object creation, after setting all properties.
function fidKeepBrightest_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fidKeepBrightest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fidRelativeHeight_Callback(hObject, eventdata, handles)
% hObject    handle to fidRelativeHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fidRelativeHeight as text
%        str2double(get(hObject,'String')) returns contents of fidRelativeHeight as a double


% --- Executes during object creation, after setting all properties.
function fidRelativeHeight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fidRelativeHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fidMinHBratio_Callback(hObject, eventdata, handles)
% hObject    handle to fidMinHBratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fidMinHBratio as text
%        str2double(get(hObject,'String')) returns contents of fidMinHBratio as a double


% --- Executes during object creation, after setting all properties.
function fidMinHBratio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fidMinHBratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fidMinAHratio_Callback(hObject, eventdata, handles)
% hObject    handle to fidMinAHratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fidMinAHratio as text
%        str2double(get(hObject,'String')) returns contents of fidMinAHratio as a double


% --- Executes during object creation, after setting all properties.
function fidMinAHratio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fidMinAHratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fidMaxUncert_Callback(hObject, eventdata, handles)
% hObject    handle to fidMaxUncert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fidMaxUncert as text
%        str2double(get(hObject,'String')) returns contents of fidMaxUncert as a double


% --- Executes during object creation, after setting all properties.
function fidMaxUncert_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fidMaxUncert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
