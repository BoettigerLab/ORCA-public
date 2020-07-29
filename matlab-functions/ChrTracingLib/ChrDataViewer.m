function varargout = ChrDataViewer(varargin)
% CHRDATAVIEWER MATLAB code for ChrDataViewer.fig
%      CHRDATAVIEWER, by itself, creates a new CHRDATAVIEWER or raises the existing
%      singleton*.
%
%      H = CHRDATAVIEWER returns the handle to a new CHRDATAVIEWER or the handle to
%      the existing singleton*.
%
%      CHRDATAVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHRDATAVIEWER.M with the given input arguments.
%
%      CHRDATAVIEWER('Property','Value',...) creates a new CHRDATAVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ChrDataViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ChrDataViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ChrDataViewer

% Last Modified by GUIDE v2.5 10-Nov-2017 19:35:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ChrDataViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @ChrDataViewer_OutputFcn, ...
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


% --- Executes just before ChrDataViewer is made visible.
function ChrDataViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ChrDataViewer (see VARARGIN)

    global DV
    if isempty(DV) || isempty(DV{1})
        DV = cell(1,1);
    else
        DV = [DV;cell(1,1)];
    end
    id = length(DV);
    set(handles.DVinstance,'String',['inst id',num2str(id)]);
    handles.id = id;

    DV{handles.id}.dataFolder = '';
    DV{handles.id}.currentSpot = 1;
    DV{handles.id}.tableKeepReject = table();

    % defaults
    defaults = cell(0,3);
    % Fiducial alignment defaults
    defaults(end+1,:) = {'upsample','positive',4};  % 8 for accuracy 2 for speed
    defaults(end+1,:) = {'maxShiftXY','positive',4};
    defaults(end+1,:) = {'maxShiftZ','positive',6};
    % Fiducial fitting defaults
    defaults(end+1,:) = {'fidMinPeakHeight', 'positive', 200};
    defaults(end+1,:) = {'fidMinSep', 'nonnegative', 5};  % Min separation between peaks in pixels.  Closer than this will be averaged
    defaults(end+1,:) = {'fidKeepBrightest','integer',2};  % Max number of peaks to allow
    % Data fitting defaults
    defaults(end+1,:) = {'datBoxXY','positive',5}; % box radius in pixels
    defaults(end+1,:) = {'datBoxZ','positive',7}; % box radius in pixels
    defaults(end+1,:) = {'datMinPeakHeight', 'positive', 500};
    defaults(end+1,:) = {'datMaxFitWidth', 'positive', 6};
    defaults(end+1,:) = {'datMinSep', 'nonnegative', 3};
    defaults(end+1,:) = {'datMinHBratio','nonnegative',1.2}; % peak value over background value
    defaults(end+1,:) = {'datMinAHratio','nonnegative',.25}; % fitted height over background vs peak value
    defaults(end+1,:) = {'datMaxUncert','nonnegative',2}; % pixels
    pars = ParseVariableArguments([],defaults,'ChrTracer2_FitSpots');

    DV{handles.id}.fitPars = pars;

% Choose default command line output for ChrDataViewer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ChrDataViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ChrDataViewer_OutputFcn(hObject, eventdata, handles)  %#ok<*INUSL>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DV;
% Get default command line output from handles structure
varargout{1} = handles.output;



function EditDataFolder_Callback(hObject, eventdata, handles) %#ok<*DEFNU,*INUSD>
    global DV

    DV{handles.id}.dataFolder = get(handles.EditDataFolder,'String');
    if isempty(DV{handles.id}.dataFolder) 
        error('Please specify a data folder containing i4d files');
    end
    % fovFolder = 'J:\2017-09-19_IMR90_chr21_CT2out\fov001\'
    fovFolder = DV{handles.id}.dataFolder; 
    i4ds = cellstr(ls([fovFolder,'*AlignedData.i4d']));
    if isempty(i4ds)
        error('Please specify a data folder containing i4d files');
    end
    i4ds = strcat(fovFolder,i4ds);
    DV{handles.id}.numCells = length(i4ds); 
    DV{handles.id}.i4ds = i4ds;
    
    s = strfind(fovFolder,'fov');
    DV{handles.id}.fov = str2double(fovFolder(s+3:s+5));

    set(handles.SliderCellNum,'Max',DV{handles.id}.numCells);
    set(handles.SliderCellNum,'Min',1);
    set(handles.SliderCellNum,'Value',1);
    set(handles.SliderCellNum,'SliderStep',[1/DV{handles.id}.numCells,5/DV{handles.id}.numCells]);
    guidata(hObject, handles);

function EditCellNum_Callback(hObject, eventdata, handles)
    global DV
    DV{handles.id}.currentSpot = str2double(get(handles.EditCellNum,'String'));


% --- Executes on slider movement.
function SliderCellNum_Callback(hObject, eventdata, handles)
    global DV
    currentSpot = round(get(handles.SliderCellNum,'Value'));
    set(handles.EditCellNum,'String',num2str(currentSpot));
    DV{handles.id}.currentSpot = currentSpot;
    guidata(hObject, handles);
    ButtonLoadData_Callback(hObject, eventdata, handles);


%=========================================================================%
% --- Executes on button press in ButtonLoadData.
function ButtonLoadData_Callback(hObject, eventdata, handles)

    global DV
    
    
    % read in parameters
    showXY = get(handles.CheckDisplayXY,'Value');
    showXZ = get(handles.CheckDisplayXZ,'Value');
    loadDat = get(handles.CheckLoadData,'Value');
    loadFid = get(handles.CheckLoadFid,'Value');
    loadFits = get(handles.CheckLoadFits,'Value');
    showLevels = get(handles.CheckDisplayLevels,'Value');
    showMap = get(handles.CheckDisplayMap,'Value');
    showPolymer = get(handles.CheckDisplayPolymer,'Value');
    showOverlay = get(handles.CheckDisplayOverlays,'Value');

    i4ds = DV{handles.id}.i4ds;
    c = DV{handles.id}.currentSpot;
    [spotFolder,spotName] = fileparts(regexprep(i4ds{c},'_AlignedData.i4d',''));
    nameX = strfind(spotName,'(');
    nameY = strfind(spotName,',');
    DV{handles.id}.locusX = str2double( spotName(nameX+1:nameY-1) );
    DV{handles.id}.locusY = str2double( spotName(nameY+1:end-1) );
    DV{handles.id}.spotName = spotName;

    %---------------------- load stuff --------------------------------% 
    if loadDat
        datMat = ReadImage4D(i4ds{c},'showPlots',false);
        nIms = size(datMat,4);
    else
        datMat = [];
    end
    DV{handles.id}.datMat = datMat;
    if loadFits
        fits = readtable( regexprep(i4ds{c},'AlignedData.i4d','fits.csv') );
    else
        fits = [];
    end
    DV{handles.id}.fits = fits;
    if loadFid
        fidMat = ReadImage4D(regexprep(i4ds{c},'AlignedData','AlignedFid'),'showPlots',false);
        nIms = size(datMat,4);
    else
        fidMat = [];
    end
    DV{handles.id}.fidMat = fidMat;
    tileLabels = cellstr( num2str((1:nIms)') );
    %---------------------- end load stuff --------------------------------% 
    %---------------------- display stuff-----------------------------%
    if showXY
        figure(101); clf;
        if loadFid && loadDat
            subplot(2,1,1); PlotProjection4D(fidMat,'fits',[],'projection','xy','tileLabels',tileLabels); 
            subplot(2,1,2); PlotProjection4D(datMat,'fits',fits,'projection','xy','tileLabels',tileLabels); 
        elseif loadFid
            PlotProjection4D(fidMat,'fits',[],'projection','xy','tileLabels',tileLabels);  
        elseif loadDat
            PlotProjection4D(datMat,'fits',fits,'projection','xy','tileLabels',tileLabels);  
        end
    end
    if showXZ
        figure(102); clf;
        if loadFid && loadDat
            subplot(2,1,1); PlotProjection4D(fidMat,'fits',[],'projection','xz','tileLabels',tileLabels); 
            subplot(2,1,2); PlotProjection4D(datMat,'fits',fits,'projection','xz','tileLabels',tileLabels); 
        elseif loadFid
            PlotProjection4D(fidMat,'fits',[],'projection','xz','tileLabels',tileLabels); 
        elseif loadDat
            PlotProjection4D(datMat,'fits',fits,'projection','xz','tileLabels',tileLabels); 
        end
    end
    
    if showOverlay
        if showXY
            figure(103); clf;
            if loadFid && loadDat
                subplot(1,2,1); Ncolor((2/nIms)*IncreaseContrast(squeeze(max(fidMat,[],3))));
                subplot(1,2,2); Ncolor((2/nIms)*IncreaseContrast(squeeze(max(datMat,[],3))));
            elseif loadFid
                Ncolor((2/nIms)*IncreaseContrast(squeeze(max(fidMat,[],3))));
            elseif loadDat
                Ncolor((2/nIms)*IncreaseContrast(squeeze(max(datMat,[],3))));
            end
        end
        if showXZ
           figure(104); clf;
            if loadFid && loadDat
                subplot(1,2,1); Ncolor((2/nIms)*IncreaseContrast(squeeze(max(permute( fidMat,[1,3,2,4]),[],3))));
                subplot(1,2,2); Ncolor((2/nIms)*IncreaseContrast(squeeze(max(permute( datMat,[1,3,2,4]),[],3))));
            elseif loadFid
                Ncolor((2/nIms)*IncreaseContrast(squeeze(max(permute( fidMat,[1,3,2,4]),[],3))));
            elseif loadDat
                Ncolor((2/nIms)*IncreaseContrast(squeeze(max(permute( datMat,[1,3,2,4]),[],3))));
            end
        end
         if showXY && showXZ
               figure(103); clf;
            if loadFid && loadDat
                subplot(2,2,1); Ncolor((2/nIms)*IncreaseContrast(squeeze(max(fidMat,[],3))));
                subplot(2,2,2); Ncolor((2/nIms)*IncreaseContrast(squeeze(max(datMat,[],3))));
                subplot(2,2,3); Ncolor((2/nIms)*IncreaseContrast(squeeze(max(permute( fidMat,[1,3,2,4]),[],3))));
                subplot(2,2,4); Ncolor((2/nIms)*IncreaseContrast(squeeze(max(permute( datMat,[1,3,2,4]),[],3))));
            elseif loadFid
                subplot(1,2,1); Ncolor((2/nIms)*IncreaseContrast(squeeze(max(fidMat,[],3))))
                subplot(1,2,2); Ncolor((2/nIms)*IncreaseContrast(squeeze(max(permute( fidMat,[1,3,2,4]),[],3))));
            elseif loadDat
                subplot(1,2,1); Ncolor((2/nIms)*IncreaseContrast(squeeze(max(fidMat,[],3))));
                subplot(1,2,2); Ncolor((2/nIms)*IncreaseContrast(squeeze(max(permute( datMat,[1,3,2,4]),[],3))));
            end
         end
    end

    


    if showLevels
        pixMax = zeros(nIms,1);
        for i=1:nIms
            datIm = datMat(:,:,:,i);
            pixMax(i) = max(datIm(:));
        end
        figure(105); clf; 
        bar(pixMax);
    end

    if showMap && loadFits && ~isempty(fits)
        if max(fits.idx) == 1
            map = squareform(pdist([fits.x,fits.y,fits.z]));
            figure(106); clf;
            imagesc(map); colorbar;
        else
            i = fits.idx==1;
            k = fits.idx==2;
            map1 = squareform(pdist([fits.x(i),fits.y(i),fits.z(i)]));
            map2 = squareform(pdist([fits.x(k),fits.y(k),fits.z(k)]));
            figure(106); clf;
            subplot(1,2,1); imagesc(map1); colorbar;
            subplot(1,2,2); imagesc(map2); colorbar;
        end
        
    end

    if showPolymer && loadFits && ~isempty(fits)
        if max(fits.idx) == 1
            figure(107); clf;
            PlotPolymerTube([fits.x,fits.y,fits.z],...
                'sphereRadius',60,...
                'maxJump',600,...
                'number',true);
        else
            figure(107); clf;
            i = fits.idx==1;
            k = fits.idx==2;
            subplot(1,2,1);
            PlotPolymerTube([fits.x(i),fits.y(i),fits.z(i)],...
                'sphereRadius',60,...
                'maxJump',600,...
                'number',true);
            subplot(1,2,2);
            PlotPolymerTube([fits.x(k),fits.y(k),fits.z(k)],...
                'sphereRadius',60,...
                'maxJump',600,...
                'number',true);
        end
    end
    
        showFidFOV = get(handles.CheckLoadFOVfid,'Value');
    if showFidFOV
        fovstr = num2str(DV{handles.id}.fov,'%03d');
        fidDax = [spotFolder,'\..\fov',fovstr,'_h001_fid.dax'];
        dax = ReadDax(fidDax);
        fidIm = max(dax,[],3);
        figure(108); clf; imagesc(fidIm); colormap(gray); hold on;
        plot(DV{handles.id}.locusX,DV{handles.id}.locusY,'yo');
    end
    
    
    %----------------------end display stuff-----------------------------%
%========================================================================%





% --- Executes on button press in ButtonRefitAll.
function ButtonRefitAll_Callback(hObject, eventdata, handles)
% refit data
    global DV;
    global scratchPath;

    datMat = DV{handles.id}.datMat;
    fidMat = DV{handles.id}.fidMat;
    fits = DV{handles.id}.fits;
    c = DV{handles.id}.currentSpot;

    pars = DV{handles.id}.fitPars;
    pars.fov = DV{handles.id}.fov; 
    pars.saveFolder = scratchPath;
    pars.overwrite = true;
    pars.showPlots = true; 

    % need to convert 4D mat to cell
    nChns = size(fidMat,4);
    fidSpts = cell(nChns,1);
    for n=1:nChns
        fidSpts{n} = fidMat(:,:,:,n);
    end

    % need to convert 4D mat to cell
    nChns = size(datMat,4);
    datSpts = cell(nChns,1);
    for n=1:nChns
        datSpts{n} = datMat(:,:,:,n);
    end
    locusXY = zeros(c,2);
    locusXY(c,:) = [DV{handles.id}.locusX,DV{handles.id}.locusY];
    newFits = ChrTracer2_FitSpots(fidSpts,datSpts,c,'parameters',pars,'lociXY',locusXY);
    DV{handles.id}.newFits = newFits; 

% --- Executes on button press in ButtonFitParameters.
function ButtonFitParameters_Callback(hObject, eventdata, handles)
% create a popup menu to specify new fit parameters
global DV
DV{handles.id}.fitPars = ChrTracer_ParameterGUI(DV{handles.id}.fitPars);

% --- Executes on button press in ButtonSaveFits.
function ButtonSaveFits_Callback(hObject, eventdata, handles)
% save table to disk
global DV
tableName = [DV{handles.id}.spotName,'_fits.csv'];
writetable(DV{handles.id}.newFits,tableName);

% --- Executes during object creation, after setting all properties.
function EditDataFolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditDataFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function EditCellNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditCellNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes during object creation, after setting all properties.
function SliderCellNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SliderCellNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in CheckDisplayXZ.
function CheckDisplayXZ_Callback(hObject, eventdata, handles)
% hObject    handle to CheckDisplayXZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckDisplayXZ


% --- Executes on button press in CheckDisplayXY.
function CheckDisplayXY_Callback(hObject, eventdata, handles)
% hObject    handle to CheckDisplayXY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckDisplayXY


% --- Executes on button press in CheckDisplayOverlays.
function CheckDisplayOverlays_Callback(hObject, eventdata, handles)
% hObject    handle to CheckDisplayOverlays (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckDisplayOverlays


% --- Executes on button press in CheckLoadFid.
function CheckLoadFid_Callback(hObject, eventdata, handles)
% hObject    handle to CheckLoadFid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckLoadFid


% --- Executes on button press in CheckLoadData.
function CheckLoadData_Callback(hObject, eventdata, handles)
% hObject    handle to CheckLoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckLoadData


% --- Executes on button press in CheckLoadFits.
function CheckLoadFits_Callback(hObject, eventdata, handles)
% hObject    handle to CheckLoadFits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckLoadFits


% --- Executes on button press in ButtonKeep.
function ButtonKeep_Callback(hObject, eventdata, handles)
% record a good fit
global DV
addTable = table({DV{handles.id}.spotName},...
                DV{handles.id}.fov,...
                DV{handles.id}.currentSpot,...
                true,...
                'VariableNames',{'name','fov','spot','keep'});
DV{handles.id}.tableKeepReject = cat(1,DV{handles.id}.tableKeepReject,addTable);
disp(['added ',DV{handles.id}.spotName,' to keep list']);


% --- Executes on button press in ButtonReject.
function ButtonReject_Callback(hObject, eventdata, handles)
% record a bad fit
global DV
addTable = table({DV{handles.id}.spotName},...
                DV{handles.id}.fov,...
                DV{handles.id}.currentSpot,...
                false,...
                'VariableNames',{'name','fov','spot','keep'});
DV{handles.id}.tableKeepReject = cat(1,DV{handles.id}.tableKeepReject,addTable);
disp(['added ',DV{handles.id}.spotName,' to reject list']);

% --- Executes on button press in ButtonRefitData.
function ButtonRefitData_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonRefitData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in CheckDisplayLevels.
function CheckDisplayLevels_Callback(hObject, eventdata, handles)
% hObject    handle to CheckDisplayLevels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckDisplayLevels


% --- Executes on button press in CheckDisplayMap.
function CheckDisplayMap_Callback(hObject, eventdata, handles)
% hObject    handle to CheckDisplayMap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckDisplayMap


% --- Executes on button press in CheckDisplayPolymer.
function CheckDisplayPolymer_Callback(hObject, eventdata, handles)
% hObject    handle to CheckDisplayPolymer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckDisplayPolymer






% --- Executes on button press in CheckLoadFOVdata.
function CheckLoadFOVdata_Callback(hObject, eventdata, handles)
% hObject    handle to CheckLoadFOVdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckLoadFOVdata


% --- Executes on button press in CheckLoadFOVfid.
function CheckLoadFOVfid_Callback(hObject, eventdata, handles)
% hObject    handle to CheckLoadFOVfid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckLoadFOVfid


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DV;
% Hint: delete(hObject) closes the figure
try
    DV(handles.id) = [];
catch
end
delete(hObject);
