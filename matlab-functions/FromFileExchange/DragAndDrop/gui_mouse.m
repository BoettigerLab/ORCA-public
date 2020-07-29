function varargout = gui_mouse(varargin)
% GUI_MOUSE M-file for gui_mouse.fig
%      GUI_MOUSE, by itself, creates a new GUI_MOUSE or raises the existing
%      singleton*.
%
%      H = GUI_MOUSE returns the handle to a new GUI_MOUSE or the handle to
%      the existing singleton*.
%
%      GUI_MOUSE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_MOUSE.M with the given input arguments.
%
%      GUI_MOUSE('Property','Value',...) creates a new GUI_MOUSE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_mouse_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_mouse_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_mouse

% Last Modified by GUIDE v2.5 24-Jan-2012 21:44:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_mouse_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_mouse_OutputFcn, ...
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


% --- Executes just before gui_mouse is made visible.
function gui_mouse_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_mouse (see VARARGIN)

% Choose default command line output for gui_mouse
handles.output = hObject;
set(hObject, 'WindowButtonMotionFcn', @drag_and_drop);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_mouse wait for user response (see UIRESUME)
% uiwait(handles.fig_mouse);


% --- Outputs from this function are returned to the command line.
function varargout = gui_mouse_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on mouse motion over figure - except title and menu.
function fig_mouse_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to fig_mouse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pos = get(hObject, 'currentpoint'); % get mouse location on figure
x = pos(1); y = pos(2); % assign locations to x and y
set(handles.lbl_x, 'string', ['x loc:' num2str(x)]); % update text for x loc
set(handles.lbl_y, 'string', ['y loc:' num2str(y)]); % update text for y loc 
% get position information of the uicontrol 
bounds = get(handles.lbl_target,'position');
lx = bounds(1); ly = bounds(2);
lw = bounds(3); lh = bounds(4);
% test to see if mouse is within the uicontrol. 
if x >= lx && x <= (lx + lw) && y >= ly && y <= (ly + lh)
    set(handles.lbl_target, 'string', 'IN');
    set(handles.lbl_target, 'backgroundcolor', 'red');
    
    % change to hand shape pointer     
    setfigptr('hand', handles.fig_mouse);
else
    set(handles.lbl_target,'string', 'OUT');
    set(handles.lbl_target, 'backgroundcolor', 'green');    
    
     % change back to normal arrow pointer 
    setfigptr('arrow', handles.fig_mouse);
end

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function fig_mouse_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to fig_mouse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pos = get(hObject, 'currentpoint'); % get mouse location on figure
x = pos(1); y = pos(2); % assign locations to x and y
set(handles.lbl_last_action, 'string', ['Mouse pressed @ X: ', num2str(x), ', Y: ', num2str(y)]);


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function fig_mouse_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to fig_mouse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pos = get(hObject, 'currentpoint'); % get mouse location on figure
x = pos(1); y = pos(2); % assign locations to x and y
set(handles.lbl_last_action, 'string', ['Mouse released @ X: ', num2str(x), ', Y: ', num2str(y)]);
