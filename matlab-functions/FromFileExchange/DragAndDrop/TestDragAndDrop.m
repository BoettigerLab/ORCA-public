function varargout = TestDragAndDrop(varargin)
% TESTDRAGANDDROP MATLAB code for TestDragAndDrop.fig
%      TESTDRAGANDDROP, by itself, creates a new TESTDRAGANDDROP or raises the existing
%      singleton*.
%
%      H = TESTDRAGANDDROP returns the handle to a new TESTDRAGANDDROP or the handle to
%      the existing singleton*.
%
%      TESTDRAGANDDROP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TESTDRAGANDDROP.M with the given input arguments.
%
%      TESTDRAGANDDROP('Property','Value',...) creates a new TESTDRAGANDDROP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TestDragAndDrop_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TestDragAndDrop_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TestDragAndDrop

% Last Modified by GUIDE v2.5 15-Apr-2014 11:04:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TestDragAndDrop_OpeningFcn, ...
                   'gui_OutputFcn',  @TestDragAndDrop_OutputFcn, ...
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


% --- Executes just before TestDragAndDrop is made visible.
function TestDragAndDrop_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TestDragAndDrop (see VARARGIN)

% Choose default command line output for TestDragAndDrop
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes TestDragAndDrop wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TestDragAndDrop_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
