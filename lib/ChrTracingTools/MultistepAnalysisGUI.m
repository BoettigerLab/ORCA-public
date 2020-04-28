function varargout = MultistepAnalysisGUI(varargin)
% MULTISTEPANALYSISGUI MATLAB code for MultistepAnalysisGUI.fig
%      MULTISTEPANALYSISGUI, by itself, creates a new MULTISTEPANALYSISGUI or raises the existing
%      singleton*.
%
%      H = MULTISTEPANALYSISGUI returns the handle to a new MULTISTEPANALYSISGUI or the handle to
%      the existing singleton*.
%
%      MULTISTEPANALYSISGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MULTISTEPANALYSISGUI.M with the given input arguments.
%
%      MULTISTEPANALYSISGUI('Property','Value',...) creates a new MULTISTEPANALYSISGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MultistepAnalysisGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MultistepAnalysisGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MultistepAnalysisGUI

% Last Modified by GUIDE v2.5 25-Jun-2019 17:38:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MultistepAnalysisGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @MultistepAnalysisGUI_OutputFcn, ...
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


% --- Executes just before MultistepAnalysisGUI is made visible.
function MultistepAnalysisGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MultistepAnalysisGUI (see VARARGIN)

% Choose default command line output for MultistepAnalysisGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes MultistepAnalysisGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MultistepAnalysisGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%===============
% Edit the figure file to create the buttons etc
%   just ignore this m-file for running the figure, all of that is handled
%   by the class which calls this figure layout.
% That class and its functions can be nested in other classes.  
%   Any functions specified here, in contrast, would only be available
%   here.

% --- Executes on button press in ButtonPars.
function ButtonPars_Callback(hObject, eventdata, handles)

% --- Executes on button press in ButtonRun.
function ButtonRun_Callback(hObject, eventdata, handles)

% --- Executes on button press in ButtonNext.
function ButtonNext_Callback(hObject, eventdata, handles)

% --- Executes on button press in ButtonBack.
function ButtonBack_Callback(hObject, eventdata, handles) %#ok<*INUSL>

