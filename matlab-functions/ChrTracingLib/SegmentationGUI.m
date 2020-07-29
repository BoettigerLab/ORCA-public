%% Use: [segmentIDs, polygonData] = SegmentationGUI(Image, Pars...)
%
%% --- Return values
% segmentIDs: n x 1 matrix of segment tags values (see 'markerTags' parameter)
% polygonData: Table structure. Pass as 'PolygonData' parameter to reload values.
%
%% --- Parameters
% Image (required first parameter): 
%    - Background image to display during segmentation
%    - alternatively, supply handle to an existing image
% Points: n x 2 matrix of x-y point coordinates to sort
% PolygonData: load ouptut from previous run to load segment data
% MarkerTags: List of bins into which points will be sorted.
% MarkerShape: Shape of plotted points
% MarkerDefaultColor: Color to display untagged points
%
%% --- Updates
% Written by Zack Cinquini, August 2017
% 10/18/17 - Updated to use matlab-figure environment for plotting.
%            This enables zoom-in zoom-out, pan etc. 
%            Also reformatted GUI and renamed as a GUI
%            - Alistair Boettiger 
% 06/13/19 - Updated to also accept existing figure instead of image. 
%            This avoids opening extra figures. 

function varargout = SegmentationGUI(varargin)
% SEGMENTATIONGUI MATLAB code for SegmentationGUI.fig
%      SEGMENTATIONGUI, by itself, creates a new SEGMENTATIONGUI or raises the existing
%      singleton*.
%
%      H = SEGMENTATIONGUI returns the handle to a new SEGMENTATIONGUI or the handle to
%      the existing singleton*.
%
%      SEGMENTATIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEGMENTATIONGUI.M with the given input arguments.
%
%      SEGMENTATIONGUI('Property','Value',...) creates a new SEGMENTATIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SegmentationGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SegmentationGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SegmentationGUI

% Last Modified by GUIDE v2.5 18-Oct-2017 22:47:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SegmentationGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SegmentationGUI_OutputFcn, ...
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


% --- Executes just before SegmentationGUI is made visible.
function SegmentationGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SegmentationGUI (see VARARGIN)
defaults = cell(0,3);
defaults(end+1,:) = {'Points','freeType',[]};
defaults(end+1,:) = {'PolygonData','freeType',{}};  % No option for type 'table'
defaults(end+1,:) = {'MarkerTags','cell',{'1','2','3','4','5','6','7','8','9','10','11','12','13','14'}};
defaults(end+1,:) = {'MarkerShape','string','.'};
defaults(end+1,:) = {'MarkerSize','integer',2};
defaults(end+1,:) = {'MarkerDefaultColor','colormap',[1,1,1]};
defaults(end+1,:) = {'figHandle','freeType',101};
pars = ParseVariableArguments(varargin(2:end),defaults,mfilename);
% -- pass forward into handles
% ------ display parameters
handles.markerTags = pars.MarkerTags;
handles.markerColors = hsv(length(handles.markerTags));
handles.markerShape = pars.MarkerShape;
handles.markerDefaultColor = pars.MarkerDefaultColor;
handles.markerSize = pars.MarkerSize;
handles.figHandle = pars.figHandle;
% ------- data
handles.ptsXY = pars.Points;
handles.npts = size(handles.ptsXY, 1);
handles.pointTags = strings(handles.npts,1);
% -------- initialization
handles.ROI = cell(0,3);  % Columns: {Polygons, Tag, Text}
handles.busyMakingROI = false;
% -- parse required input
if isempty(varargin{1}) % check that something was passed
    error('Please pass Image_Segmentation a 2D or 3-color image: Image_Segmentation(im); or a handle to an existing figure.');
end
im = varargin{1};
if ishandle(im) % if im is a figure handle
    handles.figure100 = im;
    handles.figHandle = handles.figure100; % this seems redundant
    figure(handles.figure100);
    hold on;
else % if im is a image matrix, show the matrix
    handles.figure100 = figure(handles.figHandle); clf; % Alistair added to allow zoom
    % handles.axes2 = gca; % Alistair added [moved down]
    imagesc(im);
    colormap gray;
    hold on;
end
if handles.npts > 0
    s = scatter(handles.ptsXY(:,1) ,handles.ptsXY(:,2), ...
        'Marker', handles.markerShape, ...
        'MarkerEdgeColor', handles.markerDefaultColor,...
        'SizeData',handles.markerSize);
    s.PickableParts = 'none';
end
handles.axes2 = gca; % Alistair added [moved down]
% Set up popup menu
set(handles.ROITypeMenu, 'String', pars.MarkerTags);
% ---- Load polygon data
if ~isempty(pars.PolygonData)
    nPolys = size(pars.PolygonData.Tag,1);
    % Check input tags are valid
    for i = 1:nPolys
        members = ismember(pars.PolygonData.Tag{i}, handles.markerTags);
        z = find(members==0, 1);
        if ~isempty(z)
            error('Error reading input polygon tags: tag not found in list. Make sure it matches with markerTags.');
        end
    end  
    % Load existing polygons
    for i = 1:nPolys
        polyXY = pars.PolygonData.Positions{i};
        h = impoly(handles.axes2, polyXY);  % Alistair converted to axes2 from axes1 
        % h = impoly(handles.axes1, polyXY);  % Alistair converted to axes2 from axes1 
        tag = pars.PolygonData.Tag{i};
        handles = CreateROI(h, tag, handles);
    end
    handles = button_Update_Callback(handles.button_Update, eventdata, handles);
end
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes SegmentationGUI wait for user response (see UIRESUME)
uiwait(handles.figure1);

function color = Tag2Color(tag, handles)
if strcmp(tag, "")
    color = handles.markerDefaultColor;
    return;
end
index = strcmp(tag, handles.markerTags);
color = handles.markerColors(index,:);

function handles = CreateROI(h, tag, handles)
axes(handles.axes2); % Alistair Added
polyXY = getPosition(h);
color = Tag2Color(tag, handles);
setColor(h, color);
t = text(polyXY(1,1),polyXY(1,2),tag);
FormatText(t);
handles.ROI(end+1,:) = {h, tag, t};

function FormatText(t)
t.FontSize = 10;
t.FontWeight = 'bold';
t.HorizontalAlignment = 'center';
t.Color = 'w';
t.PickableParts = 'none';  % Doesn't interact with mouse

% --- Outputs from this function are returned to the command line.
function varargout = SegmentationGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Update points
handles = button_Update_Callback(handles.button_Update, eventdata, handles);
% Generate polygonData output
polygonData = cell(size(handles.ROI,1),2);
for i = 1:size(handles.ROI,1)
    h = handles.ROI{i,1};
    tag = handles.ROI{i,2};
    polygonData(i,:) = {getPosition(h), tag};
end
polygonDataTable = cell2table(polygonData, 'VariableNames', {'Positions', 'Tag'});
varargout{1} = handles.pointTags;
varargout{2} = polygonDataTable;
delete(handles.figure1);  % Close window
% delete(handles.figure100);


% --- Executes on button press in button_New_ROI.
function button_New_ROI_Callback(hObject, eventdata, handles)
% hObject    handle to button_New_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.busyMakingROI)
    return;  % Don't push New ROI button twice
end
handles.busyMakingROI = true;
guidata(hObject, handles);
axes(handles.axes2); % Alistair Added (previously nothing here)
h = impoly();
tag = handles.markerTags{get(handles.ROITypeMenu, 'Value')};
handles = CreateROI(h, tag, handles);
handles = button_Update_Callback(handles.button_Update, eventdata, handles);
handles.busyMakingROI = false;
guidata(hObject, handles);
uiwait(handles.figure1);  % impoly uses a uiresume


% --- Executes on button press in button_Update.
function handles = button_Update_Callback(hObject, eventdata, handles)
% hObject    handle to button_Update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% -- Recalculate pointTags
axes(handles.axes2); % Alistair Added (previously nothing here)
handles.pointTags = strings(handles.npts,1);
ROIsToKeep = true(size(handles.ROI,1),1);
for i = 1:size(handles.ROI,1)   
    if ~isvalid(handles.ROI{i,1})  % Some ROIs may have been deleted.
        ROIsToKeep(i) = false;
        delete(handles.ROI{i,3});
        continue;
    end    
    polyXY = getPosition(handles.ROI{i,1});
    handles.ROI{i,3}.Position = polyXY(1,:);
    if handles.npts > 0
        result = inpolygon(handles.ptsXY(:,1), handles.ptsXY(:,2), polyXY(:,1), polyXY(:,2));
        handles.pointTags(result) = handles.ROI{i,2};
    end
end
handles.ROI = handles.ROI(ROIsToKeep,:);
% --- Recolor all points
for tag = ["", handles.markerTags]
    currPts = strcmp(tag, handles.pointTags);
    color = Tag2Color(tag, handles);
    if sum(currPts) > 0
        s = scatter(handles.ptsXY(currPts,1) ,handles.ptsXY(currPts,2), ...
            'Marker', handles.markerShape,...
            'MarkerEdgeColor', color,...
            'SizeData',handles.markerSize);
        s.PickableParts = 'none';
    end
end
guidata(hObject, handles);


% --- Executes on button press in button_Done.
function handles = button_Done_Callback(hObject, eventdata, handles)
% hObject    handle to button_Done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    uiresume(handles.figure1);  % Triggers SegmentationGUI_OutputFcn
else
    delete(handles.figure1);
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
delete(hObject);
% button_Done_Callback(handles.button_Done, eventdata, handles);


% --- Executes on selection change in ROITypeMenu.
function ROITypeMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ROITypeMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% Hints: contents = cellstr(get(hObject,'String')) returns ROITypeMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ROITypeMenu


% --- Executes during object creation, after setting all properties.
function ROITypeMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ROITypeMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
%
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
