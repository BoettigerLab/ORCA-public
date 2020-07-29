function varargout = DaxViewer(varargin)
% DAXVIEWER MATLAB code for DaxViewer.fig
%      DAXVIEWER, by itself, creates a new DAXVIEWER or raises the existing
%      singleton*.
%
%      H = DAXVIEWER returns the handle to a new DAXVIEWER or the handle to
%      the existing singleton*.
%
%      DAXVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DAXVIEWER.M with the given input arguments.
%
%      DAXVIEWER('Property','Value',...) creates a new DAXVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DaxViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DaxViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DaxViewer

% Last Modified by GUIDE v2.5 11-Jun-2017 20:57:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DaxViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @DaxViewer_OutputFcn, ...
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


% --- Executes just before DaxViewer is made visible.
function DaxViewer_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DaxViewer (see VARARGIN)

global DV
if isempty(DV)
    DV = cell(1,1);
else
    DV = [DV;cell(1,1)];
end
id = length(DV);
handles.id = id;

% defaults
DV{handles.id}.verbose = true;
DV{handles.id}.figureOffset = 100; 
DV{handles.id}.frameOffset = 1;  % Hamamatsu camera drops 1st frame
DV{handles.id}.options.frameJump = 50; % number of frames in large frame jump 
DV{handles.id}.currentFrame = 1; 
DV{handles.id}.currentFolder = '';
DV{handles.id}.colormap = colormap(gray);
DV{handles.id}.channelsPerFrame = 1;
DV{handles.id}.currentChannel = 1;

set(handles.instance,'String',['id',num2str(DV{handles.id}.figureOffset+id)]);

% set up figure
figHandle = figure(DV{handles.id}.figureOffset+1); 
figHandle.Name = 'DaxViewer';
DV{handles.id}.figHandle = figHandle;


% Channel Parameters
defaults = cell(0,3);
defaults(end+1,:) = {'minContrast','integer',0};
defaults(end+1,:) = {'maxContrast','integer',0};
parsChn1 = ParseVariableArguments([],defaults,mfilename);

DV{handles.id}.channelPars{1} = parsChn1;

% Choose default command line output for DaxViewer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DaxViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DaxViewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function EditDaxFile_Callback(hObject, eventdata, handles)




% --- Executes on button press in ButtonUpdate.
function ButtonUpdate_Callback(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to ButtonUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




function UploadDax(hObject,eventdata,handles) %#ok<*DEFNU>
global DV;
DV{handles.id}.currentFrame = 1; % reset frame to 1
DV{handles.id}.daxFile  = get(handles.EditDaxFile,'String');

% get number of frames from infFile
infFile = ReadInfoFile(DV{handles.id}.daxFile);
DV{handles.id}.numFrames =  infFile.number_of_frames/DV{handles.id}.channelsPerFrame;
DV{handles.id}.frameSize = infFile.frame_dimensions;
DV{handles.id}.xlim = [1,DV{handles.id}.frameSize(1)];
DV{handles.id}.ylim = [1,DV{handles.id}.frameSize(2)];
UpdateFrameControls(hObject,eventdata,handles);

% get contrast from image
GetDaxFrame(hObject,eventdata,handles);
for c=1:DV{handles.id}.channelsPerFrame
    im = DV{handles.id}.frame(:,:,c);
    DV{handles.id}.channelPars{c}.refPix = max(im(:));
    DV{handles.id}.currentChannel = DV{handles.id}.channelsPerFrame + 1-c;
    ButtonAutoContrast_Callback(hObject, eventdata, handles,'crop',false);
end




% -- Loads dax frame from file
function dax = GetDaxFrame(hObject,eventdata,handles)
global DV;
frameOffset = DV{handles.id}.frameOffset;
currentFrame = DV{handles.id}.currentFrame;
channelsPerFrame = DV{handles.id}.channelsPerFrame;
frameStart = frameOffset +(currentFrame-1)*channelsPerFrame + 1;
frameEnd = frameOffset +(currentFrame)*channelsPerFrame;
dax = ReadDax(DV{handles.id}.daxFile,'startFrame',frameStart,'endFrame',frameEnd,'verbose',false);
DV{handles.id}.frame = dax;


function UpdateImage(hObject,eventdata,handles)
global DV;
figure(DV{handles.id}.figureOffset+1); 
DV{handles.id}.xlim = get(gca,'xlim');
DV{handles.id}.ylim = get(gca,'ylim');
if size(DV{handles.id}.frame,3) == 1
    imagesc(DV{handles.id}.frame); 
    colormap(DV{handles.id}.colormap); 
    caxis([DV{handles.id}.channelPars{1}.minContrast,DV{handles.id}.channelPars{1}.maxContrast]);
else
    im = DV{handles.id}.frame;
    for c = 1:DV{handles.id}.channelsPerFrame
        maxPixel = double(DV{handles.id}.channelPars{c}.refPix);
        maxContrast = double(DV{handles.id}.channelPars{c}.maxContrast);
        minContrast = double(DV{handles.id}.channelPars{c}.minContrast);
        % gain = (maxPixel)/(maxContrast); 
        disp(['channel ',num2str(c),' low: ',num2str(minContrast/maxPixel),' high: ',num2str(maxContrast/maxPixel)]);
        % im(:,:,c) = uint16((im(:,:,c)- minContrast)*gain);
        im(:,:,c) = IncreaseContrast(im(:,:,c),'low',minContrast/maxPixel,'high',maxContrast/maxPixel);
    end
    DV{handles.id}.im = im;
    Ncolor(im); 
end
xlim(DV{handles.id}.xlim);
ylim(DV{handles.id}.ylim);
% 
% % for troubleshooting
% figure(2); clf; 
% subplot(2,1,1); imagesc(DV{handles.id}.im(:,:,1)); colorbar;
% subplot(2,1,2); imagesc(DV{handles.id}.im(:,:,2)); colorbar;

% --------------------------------------------------------------------
function ToolResetZoom_ClickedCallback(hObject, eventdata, handles)
global DV;
DV{handles.id}.xlim = [1,DV{handles.id}.frameSize(1)];
DV{handles.id}.ylim = [1,DV{handles.id}.frameSize(2)];
figure(DV{handles.id}.figureOffset+1);
xlim(DV{handles.id}.xlim);
ylim(DV{handles.id}.ylim);
UpdateImage(hObject,eventdata,handles)


function UpdateContrastControls(hObject,eventdata,handles)
global DV
c = DV{handles.id}.currentChannel;
maxBrightness = double(DV{handles.id}.channelPars{c}.refPix);
set(handles.SliderContrast,'Min',DV{handles.id}.channelPars{c}.minContrast);
set(handles.SliderContrast,'Max',maxBrightness);
set(handles.SliderContrast,'Value',DV{handles.id}.channelPars{c}.maxContrast); 
set(handles.SliderContrast,'SliderStep',[1/maxBrightness,100/maxBrightness]);
set(handles.EditMinContrast,'String',num2str(DV{handles.id}.channelPars{c}.minContrast));
set(handles.EditMaxContrast,'String',num2str(DV{handles.id}.channelPars{c}.maxContrast));


function UpdateFrameControls(hObject,eventdata,handles)
global DV;
set(handles.SliderFrameNum,'Min',1);
set(handles.SliderFrameNum,'Max',DV{handles.id}.numFrames);
set(handles.SliderFrameNum,'Value',DV{handles.id}.currentFrame); 
set(handles.SliderFrameNum,'SliderStep',[1/DV{handles.id}.numFrames,DV{handles.id}.options.frameJump/DV{handles.id}.numFrames]);
get(handles.SliderFrameNum,'Max')


% --- Executes on button press in ButtonAutoContrast.
function ButtonAutoContrast_Callback(hObject, eventdata, handles,varargin)
defaults = cell(0,3);
defaults(end+1,:) = {'crop','boolean',true};
pars = ParseVariableArguments(varargin,defaults,'ButtonAutoContrast');

global DV;

if pars.crop
    figure(DV{handles.id}.figureOffset+1);
    xlims = ceil(get(gca,'xlim'));
    ylims = ceil(get(gca,'ylim'));
    dax = DV{handles.id}.frame(ylims(1):ylims(2)-1,xlims(1):xlims(2)-1,:);
else
    dax = DV{handles.id}.frame;
end

for c=DV{handles.id}.currentChannel
    im = dax(:,:,c); 
    DV{handles.id}.channelPars{c}.minContrast = min(im(:));
    DV{handles.id}.channelPars{c}.maxContrast = max(im(:));
end
UpdateContrastControls(hObject,eventdata,handles);
UpdateImage(hObject,eventdata,handles)


% --- Executes on movement of Frame Slider
function SliderFrameNum_Callback(hObject, eventdata, handles)
% Update slider values 
global DV
DV{handles.id}.currentFrame = round(get(hObject,'Value'));
set(handles.EditFrameNum,'String',num2str(DV{handles.id}.currentFrame));
guidata(hObject, handles);
GetDaxFrame(hObject,eventdata,handles);
UpdateImage(hObject,eventdata,handles);

function EditFrameNum_Callback(hObject, eventdata, handles)
global DV
DV{handles.id}.currentFrame = round(gstr2double(get(hObject,'String')));
set(handles.SliderFrameNum,'Value',DV{handles.id}.currentFrame);
guidata(hObject, handles);
GetDaxFrame(hObject,eventdata,handles);
UpdateImage(hObject,eventdata,handles);

% --- Executes on update of EditMinContrast edit-text box.
function EditMinContrast_Callback(hObject, eventdata, handles)
global DV
c = DV{handles.id}.currentChannel;
DV{handles.id}.channelPars{c}.minContrast = str2double(get(handles.EditMinContrast,'String'));
GetDaxFrame(hObject,eventdata,handles);
UpdateImage(hObject,eventdata,handles);


% --- Executes on slider movement.
function SliderContrast_Callback(hObject, eventdata, handles)
% contrast slider
global DV
c = DV{handles.id}.currentChannel;
DV{handles.id}.channelPars{c}.maxContrast = round(get(handles.SliderContrast,'Value'));
set(handles.EditMaxContrast,'String',num2str(DV{handles.id}.channelPars{c}.maxContrast));
guidata(hObject, handles);
GetDaxFrame(hObject,eventdata,handles);
UpdateImage(hObject,eventdata,handles);

% --- Executes on update of EditMaxContrast edit-text box.
function EditMaxContrast_Callback(hObject, eventdata, handles)
% manual input contrast values
global DV
c = DV{handles.id}.currentChannel;
DV{handles.id}.channelPars{c}.maxContrast = str2double(get(handles.EditMaxContrast,'String'));
if DV{handles.id}.channelPars{c}.maxContrast > DV{handles.id}.channelPars{c}.refPix
    DV{handles.id}.channelPars{c}.refPix = DV{handles.id}.channelPars{c}.maxContrast;
    UpdateContrastControls(hObject, eventdata, handles);
end
set(handles.SliderContrast,'Value',DV{handles.id}.channelPars{c}.maxContrast);
guidata(hObject, handles);
GetDaxFrame(hObject,eventdata,handles);
UpdateImage(hObject,eventdata,handles);



% -- update channel to edit contrast
function EditChannelContrast_Callback(hObject, eventdata, handles)
global DV;
DV{handles.id}.currentChannel = str2double(get(handles.EditChannelContrast,'String'));
UpdateContrastControls(hObject,eventdata,handles);

% -- update channels per frame
function EditChannelsPerFrame_Callback(hObject, eventdata, handles)
global DV;
DV{handles.id}.channelsPerFrame = str2double(get(handles.EditChannelsPerFrame,'String'));

% --------------------------------------------------------------------
function ToolMaxProjectXZ_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ToolMaxProjectXZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ToolMaxProject_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ToolMaxProject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MenuSepFrames_Callback(hObject, eventdata, handles)
% hObject    handle to MenuSepFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DV
DV{handles.id}.channelsPerFrame = 1;

% --------------------------------------------------------------------
function MenuSelectChannels_Callback(hObject, eventdata, handles)
% hObject    handle to MenuSelectChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MenuLoadDax_Callback(hObject, eventdata, handles)
% hObject    handle to MenuLoadDax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DV
if isempty(DV{handles.id}.currentFolder)
    startfolder = pwd;
else
    startfolder = DV{handles.id}.currentFolder;
end

[fileName,pathName,filterIndex] = uigetfile({'*.dax','Dax Movie (*.dax)';...
    '*.*','All Files (*.*)'},'Select Dax File',startfolder);
if filterIndex ~= 0 % loading operation was not canceled 
    try
        if DV{handles.id}.verbose;
            disp(['loading ', pathName,fileName]);
        end
        DV{handles.id}.currentFolder = pathName;
        set(handles.EditDaxFile,'String',[pathName,fileName]); 
    catch
        error(['error reading from ' pathName,fileName]); 
    end
end

guidata(hObject, handles);
UploadDax(hObject, eventdata, handles);

% --------------------------------------------------------------------
function MenuSaveDax_Callback(hObject, eventdata, handles)
% hObject    handle to MenuSaveDax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes during object creation, after setting all properties.
function SliderContrast_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes during object creation, after setting all properties.
function EditMinContrast_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function EditMaxContrast_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function EditDaxFile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function EditFrameNum_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function SliderFrameNum_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --------------------------------------------------------------------
function MenuFile_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MenuDisplay_Callback(hObject, eventdata, handles)
% hObject    handle to MenuDisplay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function EditChannelsPerFrame_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object creation, after setting all properties.
function EditChannelContrast_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
