% Use: [segmentIDs, polygonData] = Image_Segmentation(Image, Pars...)

% --- Return values
% segmentIDs: n x 1 matrix of segment tags values (see 'markerTags' parameter)
% polygonData: Table structure. Pass as 'PolygonData' parameter to reload values.

% --- Parameters
% Image (required first parameter): Background image to display during segmentation
% Points: n x 2 matrix of x-y point coordinates to sort
% PolygonData: load ouptut from previous run to load segment data
% MarkerTags: List of bins into which points will be sorted.
% MarkerShape: Shape of plotted points
% MarkerDefaultColor: Color to display untagged points

function varargout = Image_Segmentation(varargin)
% IMAGE_SEGMENTATION MATLAB code for Image_Segmentation.fig
%      IMAGE_SEGMENTATION, by itself, creates a new IMAGE_SEGMENTATION or raises the existing
%      singleton*.
%
%      H = IMAGE_SEGMENTATION returns the handle to a new IMAGE_SEGMENTATION or the handle to
%      the existing singleton*.
%
%      IMAGE_SEGMENTATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGE_SEGMENTATION.M with the given input arguments.
%
%      IMAGE_SEGMENTATION('Property','Value',...) creates a new IMAGE_SEGMENTATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Image_Segmentation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Image_Segmentation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Image_Segmentation

% Last Modified by GUIDE v2.5 25-Aug-2017 15:24:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Image_Segmentation_OpeningFcn, ...
                   'gui_OutputFcn',  @Image_Segmentation_OutputFcn, ...
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


% --- Executes just before Image_Segmentation is made visible.
function Image_Segmentation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Image_Segmentation (see VARARGIN)

error('Obsolete. Please use SegmentationGUI');

% addpath(genpath('C:\Users\Zack\Desktop\Code\matlab-test\ImageSegmentation'));

defaults = cell(0,3);
defaults(end+1,:) = {'Points','freeType',[]};
defaults(end+1,:) = {'PolygonData','freeType',{}};  % No option for type 'table'
defaults(end+1,:) = {'MarkerTags','cell',{'1','2','3','4','5','6','7','8','9','10','11','12','13','14'}};
defaults(end+1,:) = {'MarkerShape','integer','x'};
defaults(end+1,:) = {'MarkerDefaultColor','array',[1,1,1]};
defaults(end+1,:) = {'imHandle','freeType',[]}; % optional, pass an image handle instead
pars = ParseVariableArguments(varargin(2:end),defaults,mfilename);

handles.markerTags = pars.MarkerTags;
handles.markerColors = hsv(length(handles.markerTags));
handles.markerShape = pars.MarkerShape;
handles.markerDefaultColor = pars.MarkerDefaultColor;

handles.ptsXY = pars.Points;
handles.npts = size(handles.ptsXY, 1);
handles.pointTags = strings(handles.npts,1);

handles.ROI = cell(0,3);  % Columns: {Polygons, Tag, Text}
handles.busyMakingROI = false;
handles.pars =pars;


if isempty(pars.imHandle)
    handles.figure101 = figure(101); clf; % Alistair added to allow zoom
else
    handles.figure101 = pars.imHandle;
end
handles.axes2 = gca; % Alistair added 
if isempty(varargin{1}) && isempty(pars.imHandle)
    % delete(handles.figure101); 
    error('Please pass Image_Segmentation a 2D or 3-color image: Image_Segmentation(im); Or pass an existing figure Image_Segmentation([],"imHandles",figHandle);');
elseif ~isempty(varargin{1}) && isempty(pars.imHandle)
    im = varargin{1};
    imagesc(im);
end
colormap gray;
hold on;
if handles.npts > 0
    axes(handles.axes2);
    s = scatter(handles.ptsXY(:,1) ,handles.ptsXY(:,2), ...
        'Marker', handles.markerShape, ...
        'MarkerEdgeColor', handles.markerDefaultColor);
    s.PickableParts = 'none';
end

% Set up popup menu
set(handles.ROITypeMenu, 'String', pars.MarkerTags);

% Load polygon data
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

% UIWAIT makes Image_Segmentation wait for user response (see UIRESUME)
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
function varargout = Image_Segmentation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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
if isempty(handles.pars.imHandle) % if Image_Segmentation created a new figure, delete it when closing.
    delete(handles.figure101);
end

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

% Recalculate pointTags
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

% Recolor all points
for tag = ["", handles.markerTags]
    currPts = strcmp(tag, handles.pointTags);
    color = Tag2Color(tag, handles);
    if sum(currPts) > 0
        s = scatter(handles.ptsXY(currPts,1) ,handles.ptsXY(currPts,2), ...
            'Marker', handles.markerShape, 'MarkerEdgeColor', color);
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
    uiresume(handles.figure1);  % Triggers Image_Segmentation_OutputFcn
else
    delete(handles.figure1);
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);

% button_Done_Callback(handles.button_Done, eventdata, handles);


% --- Executes on selection change in ROITypeMenu.
function ROITypeMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ROITypeMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ROITypeMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ROITypeMenu


% --- Executes during object creation, after setting all properties.
function ROITypeMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ROITypeMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
