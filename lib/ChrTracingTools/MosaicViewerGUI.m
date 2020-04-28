function varargout = MosaicViewerGUI(varargin)
% MosaicViewerGUI(imageTiles,tileULs);
% MOSAICVIEWERGUI MATLAB code for MosaicViewerGUI.fig
%      MOSAICVIEWERGUI, by itself, creates a new MOSAICVIEWERGUI or raises the existing
%      singleton*.
%
%      H = MOSAICVIEWERGUI returns the handle to a new MOSAICVIEWERGUI or the handle to
%      the existing singleton*.
%
%      MOSAICVIEWERGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MOSAICVIEWERGUI.M with the given input arguments.
%
%      MOSAICVIEWERGUI('Property','Value',...) creates a new MOSAICVIEWERGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MosaicViewerGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MosaicViewerGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MosaicViewerGUI

% Last Modified by GUIDE v2.5 27-Jun-2019 14:58:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MosaicViewerGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @MosaicViewerGUI_OutputFcn, ...
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


% --- Executes just before MosaicViewerGUI is made visible.
function MosaicViewerGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MosaicViewerGUI (see VARARGIN)

% a global for 
global MV NG
% % allow multiple instances to run
% if isempty(MV) 
%     MV = cell(1,1);
% else
%     MV = [MV;cell(1,1)];
% end
% id = length(MV);
MV = []; % clear MV
NG{1} = []; % clear the NcolorGUI channel 1 for use.  
NG{1}.NcolorGUI_handle=0;
id = 1; % enforce singleton
% set(handles.instance,'String',['inst id',num2str(id)]);
handles.id = id;
handles.busy = false;

% -------- Parse variable arguments
defaults = cell(0,3);
defaults(end+1,:) = {'figHandle','freeType',100};
defaults(end+1,:) = {'mosaicContrastLow','fraction',.1};
defaults(end+1,:) = {'mosaicContrastHigh','fraction',.9999};
defaults(end+1,:) = {'displaySize','integer',2E3};
defaults(end+1,:) = {'displayScale','fraction',0}; % 0 for auto-display size
defaults(end+1,:) = {'names','cell',{}};
defaults(end+1,:) = {'verbose','boolean',true}; % doesn't propegate yet 
defaults(end+1,:) = {'uiwait','boolean',true}; 
defaults(end+1,:) = {'reZero','boolean',true}; 
% MV render parameters
defaults(end+1,:) = {'method',{'mean','sum','edgeBlur','last','first'},'edgeBlur'};
defaults(end+1,:) = {'flatten','boolean',false};
defaults(end+1,:) = {'padMosaic','nonnegative',0};
% mosaic tile adjustments
defaults(end+1,:) = {'transpose','boolean',false};
defaults(end+1,:) = {'fliplr','boolean',false};
defaults(end+1,:) = {'flipud','boolean',false};
defaults(end+1,:) = {'xyShifts','freeType',[]}; % array or empty
pars = ParseVariableArguments(varargin(3:end),defaults,mfilename);
MV{handles.id}.inputPars = pars;

% --------- Parse required arguments
if isempty(varargin{1}) || isempty(varargin{2})
    error('Please pass MosaicViewerGUI a cell-array of images and a vector of upper left coordinates');
end
imageTiles = varargin{1};
ulPositions = round(varargin{2}); % enforce integer values
if pars.reZero
    minXY = min(ulPositions);
    ulPositions = [ulPositions(:,1)-minXY(1)+1,ulPositions(:,2)-minXY(2)+1];
    % MV{handles.id}.reZero = minXY;
% else
%     MV{handles.id}.reZero = [0,0];
end
% full scale, original data
nTiles = length(imageTiles);
if isempty(pars.xyShifts)
    xyShifts = zeros(nTiles,2);
else
    xyShifts = pars.xyShifts;
end

% down-sample image for speed
%   this is now done once, to the global image
MV{handles.id}.fullImageTiles = imageTiles; 
MV{handles.id}.fullUlPositions = ulPositions;
[imageTilesSmall, ulsSmall, downsample] = ...
         ComputeDownsampling(imageTiles,ulPositions,pars.displaySize,'displayScale',pars.displayScale);
% these will be all the tiles
MV{handles.id}.allImageTiles = imageTilesSmall; 
MV{handles.id}.allUlPositions = ulsSmall;
MV{handles.id}.allTileIndices = 1:nTiles;
MV{handles.id}.allXYshifts = xyShifts/downsample;
MV{handles.id}.allManualShifts = zeros(nTiles,2);
% these will be only the tiles in the field of view 
% (originally also all the tiles)
MV{handles.id}.imageTiles = MV{handles.id}.allImageTiles;
MV{handles.id}.tileIndices = MV{handles.id}.allTileIndices;
MV{handles.id}.ulPositions = MV{handles.id}.allUlPositions;
MV{handles.id}.xyShifts = MV{handles.id}.allXYshifts;
MV{handles.id}.manualShifts = MV{handles.id}.allManualShifts;

MV{handles.id}.masterDownsample = downsample;

 
% -----------set up dynamic figure
handles.figureH = figure(pars.figHandle); close; pause(.01); 
handles.figureH = figure(pars.figHandle); pause(.01); 
handles.figureH.Name = 'MosaicViewerGUI';
set(handles.figureH,'KeyPressFcn',{@FigKeyPress,handles});
ResetShifts(hObject,eventdata,handles);

% ------------set up directions
textDir = {'Use the "Select Tile" button to reposition tiles in the mosaic.';
    'Use the "Show Boxes" button to outline tile positions'};
set(handles.TextDir,'String',textDir);

% ----------- Render Mosaic Data
MV{handles.id}.zoom.on = true;
if ~isempty(MV{handles.id}.imageTiles)
    RebuildMosaic(hObject, eventdata, handles);
end
MV{handles.id}.fullMosaicIm = MV{handles.id}.mosaicIm;

% Choose default command line output for MosaicViewerGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% % UIWAIT makes MosaicViewerGUI wait for user response (see UIRESUME)
if pars.uiwait % nargout > 0 
 uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = MosaicViewerGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global MV;
varargout{1} = MV{handles.id}.allManualShifts;
h = figure(handles.figureH);
delete(h);
delete(handles.figure1);  % Close window


function [imageTilesSmall, ulsSmall, downsample, dataSize] = ...
         ComputeDownsampling(imageTiles,uls,displaySize,varargin)
     defaults = cell(0,3);
     defaults(end+1,:) = {'displayScale','fraction',0};
     pars = ParseVariableArguments(varargin,defaults,'ComputeDownsampling');
    xmin = min(uls(:,1));
    xmax = max(uls(:,1));
    ymin = min(uls(:,2));
    ymax = max(uls(:,2));
    [h_i,w_i] = size(imageTiles{1});
    h_m = round(ymax-ymin+1+3*h_i);
    w_m = round(xmax-xmin+1+3*w_i);
    dataSize = round(sqrt(h_m*w_m));
    if pars.displayScale == 0
        downsample = max(ceil(dataSize/displaySize),1);
    else
        downsample = 1/pars.displayScale;
    end
    [nTiles,nChns] = size(imageTiles);
    imageTilesSmall = cell(nTiles,nChns);
    ulsSmall = ceil(uls/downsample);
    for m=1:nTiles
        for n=1:nChns
            imageTilesSmall{m,n} = imresize(imageTiles{m,n},1/downsample);
        end
    end

% ---  Renders the mosaic image from the current imageTiles and uls
function RebuildMosaic(hObject, eventdata, handles)
global MV    
    if ~MV{handles.id}.zoom.on
        MV{handles.id}.mosaicIm = MV{handles.id}.fullMosaicIm;
    else   
        imageTiles = MV{handles.id}.imageTiles;
        uls = MV{handles.id}.ulPositions + MV{handles.id}.xyShifts;
        [nTiles,nChns] = size(imageTiles); %#ok<ASGLU>
        mosaicIms = cell(nChns,1);
        tic
        for n=1:nChns
            mosaicIms{n} = MosaicViewerRender(imageTiles,uls,...
                'parameters',MV{handles.id}.inputPars,...
                'stopOnError',true);
        end
        toc
        MV{handles.id}.mosaicIm = cat(3,mosaicIms{:});
        disp('assembly complete');
    end
    UpdateMosaicImage(handles);

    % ---  updates the figure display 
function UpdateMosaicImage(handles)
    global MV NG
    MV{handles.id}.selectTile = [];
    mosaicIm = MV{handles.id}.mosaicIm;
    % disp(MV{handles.id}.downsample) % print current display resolution
    c_m = size(mosaicIm,3);
    figure(handles.figureH); clf; % note no clear figure
    if c_m == 1
        % 
        dispMosaic = IncreaseContrast(mosaicIm,...
            'high',MV{handles.id}.inputPars.mosaicContrastHigh,...
            'low',MV{handles.id}.inputPars.mosaicContrastLow);
        imagesc(dispMosaic); 
        colormap(gray);
    else 
        if NG{1}.NcolorGUI_handle~=0
            if isfield(NG{1},'contrastMin')
                 MV{handles.id}.inputPars.mosaicContrastLow = NG{1}.contrastMin; 
                 MV{handles.id}.inputPars.mosaicContrastHigh = NG{1}.contrastMax; 
            end
            h = NG{1}.NcolorGUI_handle;
            delete(h);
            NG{1}.NcolorGUI_handle=0;
        end
        contrastMin = MV{handles.id}.inputPars.mosaicContrastLow;
        contrastMax = MV{handles.id}.inputPars.mosaicContrastHigh;
%         dispMosaic= CropEmptyPixels(mosaicIm);
        NcolorGUI(mosaicIm,...
            'names',MV{handles.id}.inputPars.names,...
            'contrastMin',contrastMin,'contrastMax',contrastMax,...
            'figHandle',handles.figureH,'instanceID',1,'restart',false);
         % disp('called');
    end
    
% --- Executes on button press in ButtonSelectTile.
function ButtonSelectTile_Callback(hObject, eventdata, handles) %#ok<*INUSL>
    global MV;
    % select a position on the current mosaic
    textDir = 'Use the cross-hairs to select tile. Then use the arrow keys to move. w/s for up/down, a/d for left right, x to save new position';
    set(handles.TextDir,'String',textDir); guidata(hObject,handles); disp(textDir);
    boxCoords = ULsToBoxCoords(MV{handles.id}.imageTiles,MV{handles.id}.ulPositions);
    UpdateMosaicImage(handles);
    [x,y] = ginput(1);
    isM = x>boxCoords(:,1) & x<boxCoords(:,2) & y>boxCoords(:,3) & y<boxCoords(:,4);
    m = find(isM);
    MV{handles.id}.selectTile = m; % save 
    % remove tile from gray image;
    imTiles = MV{handles.id}.imageTiles;
    ulPos = MV{handles.id}.ulPositions;
    imTiles(m,:) = [];
    ulPos(m,:) = [];
    [~,nChns] = size(imTiles);
    mosaicIms = cell(nChns,1);
    for n=1:nChns
        mosaicIms{n} = MosaicViewerRender(imTiles(:,n),ulPos,...
                        'parameters',MV{handles.id}.inputPars,...
                        'stopOnError',true);
    end
    mosaicIm = cat(3,mosaicIms{:});
    MV{handles.id}.mosaicIm = mosaicIm;
    ResetShifts(hObject,eventdata,handles); % set shifts back to 0 on selection of new tile;
    UpdateOverlay(hObject,eventdata,handles);
    pause(.1);
    
% --- called when figureH is in view and a key is pressed
function FigKeyPress(hObject,eventdata,handles) %#ok<*INUSD>
    global MV
    if ~isempty(MV{handles.id}.selectTile)
        dropImage = false;
        key = eventdata.Key;
        shifts = MV{handles.id}.shifts; % could also have stored this in handles
        switch(key)
            case('w') %want to move top image up
                shifts.upDownOffset = shifts.upDownOffset - 1;
            case('s') %want to move top image down
                shifts.upDownOffset = shifts.upDownOffset + 1;
            case('a') %want to move top image left
                shifts.leftRightOffset = shifts.leftRightOffset - 1;
            case('d') %want to move top image right
                shifts.leftRightOffset = shifts.leftRightOffset + 1;
             case('i') %want to move top image up
                shifts.upDownOffset = shifts.upDownOffset - 10;
            case('k') %want to move top image down
                shifts.upDownOffset = shifts.upDownOffset + 10;
            case('j') %want to move top image left
                shifts.leftRightOffset = shifts.leftRightOffset - 10;
            case('l') %want to move top image right
                shifts.leftRightOffset = shifts.leftRightOffset + 10;
            case('e')
                shifts.angle = shifts.angle + shifts.angleStep;
            case('r')
                shifts.angle = shifts.angle - shifts.angleStep;
            case('x')
                dropImage = true;
            case(28) %For right and left arrows,
            case(29) %will add this soon
        otherwise
        end
        MV{handles.id}.shifts = shifts; % could also stick into handles;
        MV{handles.id}.lastKey = key;
        if ~dropImage
            UpdateOverlay(hObject,eventdata,handles);
        else
            UpdateTilePos(hObject,eventdata,handles);
        end
    end
    
% -- Reset shifts after moving tile is complete    
function ResetShifts(hObject,eventdata,handles)
    global MV
    shifts.upDownOffset = 0;
    shifts.leftRightOffset = 0;
    shifts.angle = 0;
    MV{handles.id}.shifts = shifts;

% --- Show a two color overlay of the mosaic when moving a tile
function UpdateOverlay(hObject,eventdata,handles)
    global MV   
    im2 = nanmean(MV{handles.id}.mosaicIm,3);
    im1 = zeros(size(im2),'uint16');
    imageTiles = MV{handles.id}.imageTiles;
    boxCoords = ULsToBoxCoords(MV{handles.id}.imageTiles,MV{handles.id}.ulPositions);
    shifts = MV{handles.id}.shifts;
    m = MV{handles.id}.selectTile;
    im1((boxCoords(m,3):boxCoords(m,4))+shifts.upDownOffset,...
        (boxCoords(m,1):boxCoords(m,2))+shifts.leftRightOffset) =  imageTiles{m};
    im1 = IncreaseContrast(im1,'low',MV{handles.id}.inputPars.mosaicContrastLow,'high',MV{handles.id}.inputPars.mosaicContrastHigh);
    im2 = IncreaseContrast(im2,'low',MV{handles.id}.inputPars.mosaicContrastLow,'high',MV{handles.id}.inputPars.mosaicContrastHigh);
    mosaicOverlay = Ncolor(cat(3,im1,im2)); 
    figure(handles.figureH); clf;
    imagesc(mosaicOverlay); 

function boxCoords = ULsToBoxCoords(im,uls)
    if iscell(im)
        im = im{1};
    end
    [h_i,w_i] = size(im); 
    ulz = uls;%  [uls(:,1) - min(uls(:,1))+1, uls(:,2)-min(uls(:,2))+1];
    boxCoords = round([ulz(:,1),ulz(:,1)+w_i-1,ulz(:,2),ulz(:,2)+h_i-1]);

    
% Save the new tile position into the the data
function UpdateTilePos(hObject,eventdata,handles)
    global MV;
    shifts = MV{handles.id}.shifts;
    tileIndices = MV{handles.id}.tileIndices;
    t = MV{handles.id}.selectTile;
    m = tileIndices(t);
    MV{handles.id}.manualShifts(t,1) = MV{handles.id}.manualShifts(t,1)+shifts.leftRightOffset*MV{handles.id}.masterDownsample;
    MV{handles.id}.manualShifts(t,2) = MV{handles.id}.manualShifts(t,2)+shifts.upDownOffset*MV{handles.id}.masterDownsample;     
    MV{handles.id}.allManualShifts(m,1) = MV{handles.id}.allManualShifts(m,1)+shifts.leftRightOffset*MV{handles.id}.masterDownsample;
    MV{handles.id}.allManualShifts(m,2) = MV{handles.id}.allManualShifts(m,2)+shifts.upDownOffset*MV{handles.id}.masterDownsample;
    RebuildMosaic(hObject,eventdata,handles);
 
    
% --- Executes on button press in ButtonZoomIn.
function ButtonZoomIn_Callback(hObject, eventdata, handles)
    global MV
    %% new approach
    % zoom in by removing tiles outside current FOV -- this should make
    % replotting the image (changing color, sliding a tile, etc) faster. 
    MV{handles.id}.zoom.on = true;
    figure(handles.figureH);
    xlims = get(gca,'xlim');
    ylims = get(gca,'ylim');
    ulPositions = MV{handles.id}.ulPositions;
    inView = ulPositions(:,1) > xlims(1) & ulPositions(:,1) < xlims(2) & ...
             ulPositions(:,2) > ylims(1) & ulPositions(:,2) < ylims(2);
    if sum(inView)>0
        ulPositions = ulPositions(inView,:);
        rezero = min(ulPositions);
        ulPositions = [ulPositions(:,1)-rezero(1)+1,ulPositions(:,2)-rezero(2)+1];
        MV{handles.id}.ulPositions = ulPositions;
        MV{handles.id}.imageTiles = MV{handles.id}.imageTiles(inView,:);  % just the local parameters
        MV{handles.id}.tileIndices = MV{handles.id}.tileIndices(inView);
        MV{handles.id}.xyShifts = MV{handles.id}.xyShifts(inView,:);
        MV{handles.id}.manualShifts = MV{handles.id}.manualShifts(inView,:);
        RebuildMosaic(hObject, eventdata, handles);
    else
        warning('no tiles in view! try to zoom out a bit more'); 
    end


% --- Executes on button press in ButtonZoomOut.
function ButtonZoomOut_Callback(hObject, eventdata, handles)
    global MV
    % new approach
    MV{handles.id}.zoom.on = false;
    MV{handles.id}.imageTiles = MV{handles.id}.allImageTiles;  % just the local parameters
    MV{handles.id}.ulPositions = MV{handles.id}.allUlPositions;
    MV{handles.id}.tileIndices = MV{handles.id}.allTileIndices;
    MV{handles.id}.xyShifts = MV{handles.id}.allXYshifts;
    MV{handles.id}.manualShifts = MV{handles.id}.allManualShifts;    
    RebuildMosaic(hObject, eventdata, handles);



% --- Executes on button press in ButtonShowBoxes.
function ButtonShowBoxes_Callback(hObject, eventdata, handles)
    global MV;
    % UpdateMosaicImage(handles);
    figure(handles.figureH);
    lineHandle = findobj(gca,'Type','Rectangle');
    textHandle = findobj(gca,'Type','Text');
    if ~isempty(lineHandle) || ~isempty(textHandle)
        delete(lineHandle);
        delete(textHandle);
    else
        xy = MV{handles.id}.ulPositions;
        imTiles = MV{handles.id}.imageTiles;
        nFOV = length(imTiles);
        for f=1:nFOV
            [h,w] = size(imTiles{f});
            hold on; text(xy(f,1),xy(f,2),num2str(f),'color','c');
            r = rectangle('position',[xy(f,1),xy(f,2),w,h]);
            r.EdgeColor = [1 0 1];
        end
    end

% --- Executes on button press in ButtonDone.
function ButtonDone_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonDone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    uiresume(handles.figure1);  % Triggers SegmentationGUI_OutputFcn
else
    delete(handles.figure1);
end

% --- Executes during object deletion, before destroying properties.
function ButtonShowBoxes_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to ButtonShowBoxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% required to exist for ButtonDone now


% --------------------------------------------------------------------
function MenuEdit_Callback(hObject, eventdata, handles)
% hObject    handle to MenuEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MenuEditDisplayOptions_Callback(hObject, eventdata, handles)
% hObject    handle to MenuEditDisplayOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global MV
editPars = MV{handles.id}.inputPars;
% remove non-editable parameters 
if isfield(editPars,'uiwait')
    editPars = rmfield(editPars,{'uiwait','reZero','figHandle','xyShifts'});
end
editPars = SimpleParameterGUI(editPars);
MV{handles.id}.inputPars = editPars;
