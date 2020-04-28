function varargout = MosaicBuilderGUI(varargin)
% MOSAICBUILDERGUI MATLAB code for MosaicBuilderGUI.fig
%      MOSAICBUILDERGUI, by itself, creates a new MOSAICBUILDERGUI or raises the existing
%      singleton*.
%
%      H = MOSAICBUILDERGUI returns the handle to a new MOSAICBUILDERGUI or the handle to
%      the existing singleton*.
%
%      MOSAICBUILDERGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MOSAICBUILDERGUI.M with the given input arguments.
%
%      MOSAICBUILDERGUI('Property','Value',...) creates a new MOSAICBUILDERGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MosaicBuilderGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MosaicBuilderGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MosaicBuilderGUI

% Last Modified by GUIDE v2.5 04-Jun-2018 12:40:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MosaicBuilderGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @MosaicBuilderGUI_OutputFcn, ...
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


% --- Executes just before MosaicBuilderGUI is made visible.
function MosaicBuilderGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MosaicBuilderGUI (see VARARGIN)


global MB
% if isempty(MB)
%     MB = cell(1,1);
% else
%     MB = [MB;cell(1,1)];
% end
% id = length(MB);
MB = cell(1,1); id = 1;  % enforce singleton for now
% set(handles.instance,'String',['inst id',num2str(id)]);
handles.id = id;

% Choose default command line output for MosaicBuilderGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

DefaultParameters(hObject, eventdata, handles) 


% UIWAIT makes MosaicBuilderGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MosaicBuilderGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function DefaultParameters(hObject, eventdata, handles) 
global MB


% Find Files parameters ----------------------------------------------%
defaults = cell(0,3);
defaults(end+1,:) = {'sampleName', 'string', 'Hyb_100\Conv*.dax'};  % 'Hyb_102\Conv*.dax'
defaults(end+1,:) = {'chnDataFolder','string','S:\2018-05-15_WT_2-4hrs(HoxRNA)_(3kbDNA)\RNA_Expt\'};
defaults(end+1,:) = {'chnSaveFolder','string','S:\2018-05-15_WT_2-4hrs(HoxRNA)_(3kbDNA)_CT2out\RNA\'};
defaults(end+1,:) = {'chnTableFile','string','S:\2018-05-15_WT_2-4hrs(HoxRNA)_(3kbDNA)\RNA_Expt\rnaTableFile.xlsx'};
defaults(end+1,:) = {'verbose', 'boolean', true}; 
defaults(end+1,:) = {'veryverbose','boolean',false};
defaults(end+1,:) = {'overwrite', 'boolean', false}; 
MB{handles.id}.parsFindFiles  = ParseVariableArguments([],defaults,'MosaicBuilder FindFiles');


% Align Hybes  parameters ----------------------------------------------%
defaults = cell(0,3);
defaults(end+1,:) = {'showFigs', 'boolean', true}; 
defaults(end+1,:) = {'overwrite','boolean',false};
defaults(end+1,:) = {'selectFOV','freeType',[]};
defaults(end+1,:) = {'showCorrAlign','boolean',true};
defaults(end+1,:) = {'alignContrastHigh', 'fraction', .999}; 
defaults(end+1,:) = {'alignContrastLow', 'fraction', 0.92}; 
defaults(end+1,:) = {'angles', 'float', 0}; 
MB{handles.id}.parsAlignHybes = ParseVariableArguments([],defaults,'MosaicBuilder AlignHybes');

% Align Tiles (Dax to Mosaic) parameters ----------------------------------------------%
defaults = cell(0,3);
defaults(end+1,:) = {'corrAlign', 'boolean', true}; 
defaults(end+1,:) = {'interactive', 'boolean', true};
defaults(end+1,:) = {'trimBorder', 'integer', 10}; 
defaults(end+1,:) = {'flatten', 'boolean', true};  % use data to flatten field (sometimes bad for embryos if there is never data in the edges of the image)
defaults(end+1,:) = {'minOverlap', 'fraction', .01}; % min overlap to align
defaults(end+1,:) = {'maxShift', 'integer', 200};  % max shift in pixels
defaults(end+1,:) = {'corrAlignLow', 'fraction', .9};  % contrast for correlation alignment 
defaults(end+1,:) = {'corrAlignHigh', 'fraction', .9995}; % contrast for correlation alignment to ID objects 
defaults(end+1,:) = {'minGrad', 'float', 50}; % min gradient change to be counted as a valid shift. Make larger if you want it not to attempt to move bad matches so much
defaults(end+1,:) = {'maxSize', 'positive', 300}; % controls speed, coarse alignment resolution
defaults(end+1,:) = {'showFinal', 'boolean', true}; % just a display 
defaults(end+1,:) = {'downsampleFinal', 'integer', 10};
% more parameters
defaults(end+1,:) = {'transpose','boolean',true}; % see the mosaic parameters in the Hal parameter file used to collect the movie 
defaults(end+1,:) = {'fliplr','boolean',true};  % see the mosaic parameters in the Hal parameter file used to collect the movie
defaults(end+1,:) = {'flipud','boolean',false};  % see the mosaic parameters in the Hal parameter file used to collect the movie
% less common options, possibly remove these.
defaults(end+1,:) = {'showplot', 'boolean', true}; % need to check if this is still used
defaults(end+1,:) = {'pix_to_mm', 'positive', 6.55}; % not recommended to change this.  OLD: 6.55; % scope 1 6.45; % scope 2  
defaults(end+1,:) = {'projectAll', 'boolean', true};  % not recommended to change this
defaults(end+1,:) = {'method',{'mean','sum','edgeBlur'},'edgeBlur'}; % not recommended to change this
defaults(end+1,:) = {'sortData', 'boolean', true}; % not recommended to change this
defaults(end+1,:) = {'showProgress', 'boolean', false}; % not recommended to change this
defaults(end+1,:) = {'downsample', 'integer', 1}; % probably should remove this option
defaults(end+1,:) = {'saveCorr', 'boolean', false}; % not used probably should remove this
MB{handles.id}.parsAlignTiles = ParseVariableArguments([],defaults,'MosaicBuilder AlignTiles');

% Adjust Tiles ----------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'figHandle','freeType',100};
defaults(end+1,:) = {'mosaicContrastLow','fraction',.1};
defaults(end+1,:) = {'mosaicContrastHigh','fraction',.9999};
defaults(end+1,:) = {'displaySize','integer',1.5E3};
MB{handles.id}.parsAdjustTiles = ParseVariableArguments([],defaults,'MosaicBuilder AdjustTiles');

% Tile All Hybes ----------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'displayMosaic','boolean',true};
defaults(end+1,:) = {'downsample','integer',10};
defaults(end+1,:) = {'justNames','boolean',false};
defaults(end+1,:) = {'flatten','boolean',false}; % for display only
MB{handles.id}.parsTileAllHybes= ParseVariableArguments([],defaults,'MosaicBuilder TileAllHybes');

% Display Mosaic -----------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'figHandle','freeType',100};
defaults(end+1,:) = {'mosaicContrastLow','fraction',.1};
defaults(end+1,:) = {'mosaicContrastHigh','fraction',.9999};
defaults(end+1,:) = {'displaySize','integer',1.5E3};
MB{handles.id}.parsDisplayMosaic = ParseVariableArguments([],defaults,'MosaicBuilder DisplayMosaic');

% Save mosaic ---------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'saveFolder','string',''};
defaults(end+1,:) = {'saveName','string','mosaicTileData'};
MB{handles.id}.parsSaveMosaic = ParseVariableArguments([],defaults,'MosaicBuilder SaveMosaic');

% step names and directions ---------------------------
MB{handles.id}.currStepName = 'FindFiles';
MB{handles.id}.stepNames = {'FindFiles';
                            'AlignHybes';
                            'AlignTiles';
                            'AdjustTiles';
                            'TileAllHybes';
                            'DisplayMosaic';
                            'SaveMosaic'};
% Directions
MB{handles.id}.Directions = {'Find files: select a search folder and a file base name and select "Run Step" when ready';
	 'Align Hybes: Correct for inter-hybe drift and save aligned image data. This may take a long time for large data files.';  % It should be possible to skip the someday for speed. For now not recommended to skip
	 'Align Tiles: Correct for inter-tile drift and save new stage coordinates. To ensure an accurate reconstruction, it is recommended that the "interactive" parameter is set to true. For RNA+DNA datasets it is recommended to use only the DNA data for stage alignment. For RNA data, set "corrAlign" to false and "interactive" to false to load the data using the recorded stage coordinates.';
     'Adjust Tiles: OPTIONAL. Use MosaicGUI to select individual tiles in the current mosaic and manually adjust their positions.';
     'Tile all hybes: Construct mosaic from all hybes using the inter-tile drift determined in the previous steps. Recommended for RNA data only.  Set the Parameter "justNames" to "true" to speed up execution if you do not intend to view the mosaic now.  This will save the names for the mosaic analyzer, but you wont be able see them in the viewer in the next step.';
     'Display mosaic: OPTIONAL. Display Only. View multicolor mosaic. Use Parameters to select display channels.';
     'Save mosaic: Save new tile postions for all aligned data'};

                        
% --- Executes on button press in ButtonPars.
function ButtonPars_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonPars (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global MB
currStep = MB{handles.id}.currStepName;
if strcmp(currStep,'FindFiles')
    MB{handles.id}.parsFindFiles = SimpleParameterGUI(MB{handles.id}.parsFindFiles);
elseif strcmp(currStep,'AlignHybes')
    MB{handles.id}.parsAlignHybes = SimpleParameterGUI(MB{handles.id}.parsAlignHybes);
elseif strcmp(currStep,'AlignTiles')
    MB{handles.id}.parsAlignTiles = SimpleParameterGUI(MB{handles.id}.parsAlignTiles);
elseif strcmp(currStep,'AdjustTiles')
    MB{handles.id}.parsAdjustTiles = SimpleParameterGUI(MB{handles.id}.parsAdjustTiles);
elseif strcmp(currStep,'TileAllHybes')
    MB{handles.id}.parsTileAllHybes = SimpleParameterGUI(MB{handles.id}.parsTileAllHybes);
elseif strcmp(currStep,'DisplayMosaic')
    MB{handles.id}.parsDisplayMosaic  = SimpleParameterGUI(MB{handles.id}.parsDisplayMosaic);
elseif strcmp(currStep,'SaveMosaic')
    MB{handles.id}.parsSaveMosaic = SimpleParameterGUI(MB{handles.id}.parsSaveMosaic);
end

% --- Executes on button press in ButtonRun.
function ButtonRun_Callback(hObject, eventdata, handles)
% hObject    handle to ButtonRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global MB;

% shorthand for a few common parameters
currStep = MB{handles.id}.currStepName;
chnDataFolder = MB{handles.id}.parsFindFiles.chnDataFolder;
chnSaveFolder = MB{handles.id}.parsFindFiles.chnSaveFolder;
chnTable = readtable( MB{handles.id}.parsFindFiles.chnTableFile);
nameField = find(contains(chnTable.Properties.VariableNames,'name','IgnoreCase',true));
if length(nameField)>1
    nameField = nameField(end);
end
chnNames = chnTable{:,nameField};
allNames = cellfun(@(x) strsplit(x,','),chnNames,'UniformOutput',false);
chnNameList = cat(2,allNames{:})';

if isempty(MB{handles.id}.parsSaveMosaic.saveFolder)
    MB{handles.id}.parsSaveMosaic.saveFolder = chnSaveFolder;
end

%% Step 1: check if data is found
if strcmp(currStep,'FindFiles')
    % search for data.  Report findings.
    % if files already exist no need to create more. 
    hyb1Tag = 'h001';
    datas = cellstr(ls([chnDataFolder,filesep,MB{handles.id}.parsFindFiles.sampleName]));
    nFOV = sum(~cellfun(@isempty,datas)); 
    MB{handles.id}.nFOV = nFOV;
    result1 = ['found ',num2str(nFOV),' FOVs in raw data folder ',chnDataFolder,'    '];
    disp(result1);
    if nFOV == 0  % NO FOV
        result2 = 'Data folder is empty! check your file paths';
        disp(result2);
        set(handles.TextDir,'String',cat(2,result1,result2)); guidata(hObject, handles);
    else % Yes found some data
        procDaxFiles = cellstr(ls([chnSaveFolder,filesep,'fov*_',hyb1Tag,'_data1.dax']));
        nProcessed = length(procDaxFiles);
        result2 = ['found ',num2str(nProcessed),' processed FOVs already in target save folder ' chnSaveFolder,'    '];
        disp(result2);
        if nProcessed >= nFOV && MB{handles.id}.parsFindFiles.overwrite 
            % Data is already written but overwrite is requested  
            result3 = ['overwrite is on, existing files in ',chnSaveFolder ,'will be overwritten.']; 
            set(handles.TextDir,'String',cat(2,result1,result2,result3));
        elseif nProcessed >= nFOV && ~MB{handles.id}.parsFindFiles.overwrite
            % Data is already written, no need to process. Skip ahead.
            result3 = 'found existing data, skipping aligment.  Please Wait...';
            set(handles.TextDir,'String',cat(2,result1,result2,result3));
            pause(5); 
            ButtonNext_Callback(hObject, eventdata, handles); pause(1);
            ButtonNext_Callback(hObject, eventdata, handles); 
        elseif nProcessed < nFOV && nProcessed > 0 && ~MB{handles.id}.parsFindFiles.overwrite
            % Data was partially written, resume from here. 
            result3 = 'found some data'; 
            set(handles.TextDir,'String',cat(2,result1,result2,result3));
        else
            % no prior data existss, start fresh
            result3 = 'starting fresh';
            set(handles.TextDir,'String',cat(2,result1,result2,result3));
        end
        disp(result3);    
    end
    
%% Step 2: register and save images of each channel data
elseif strcmp(currStep,'AlignHybes')
    hyb1Tag = 'h001';
    nFOV = MB{handles.id}.nFOV;
    showFigs = MB{handles.id}.parsAlignHybes.showFigs;
    selectFOV = MB{handles.id}.parsAlignHybes.selectFOV;
    overwrite = MB{handles.id}.parsAlignHybes.overwrite;
    SetFigureSavePath(chnSaveFolder,'makeDir',true,'verbose',false);
    fidMaps = cell(nFOV,1);
    datMaps = cell(nFOV,1);
    if isempty(selectFOV) || selectFOV==0
        selectFOV = 1:nFOV;
    end
  
    for f=selectFOV         % f = 11
        try
            daxFiles = ls([chnSaveFolder,'fov',num2str(f,'%03d'),'_',hyb1Tag,'_data1.dax']);
            if isempty(daxFiles) || overwrite 
            [dataOut,pars] = ChrTracer3_LoadData(MB{handles.id}.parsFindFiles.chnTableFile,...
                'fov',f,'refHybe',1,'saveFolder',chnSaveFolder,...
                'overwrite',MB{handles.id}.parsAlignHybes.overwrite,...
                'alignContrastHigh',MB{handles.id}.parsAlignHybes.alignContrastHigh,...
                'alignContrastLow',MB{handles.id}.parsAlignHybes.alignContrastLow,...
                'angles',MB{handles.id}.parsAlignHybes.angles,...
                'showCorrAlign',MB{handles.id}.parsAlignHybes.showCorrAlign);
            % 'alignmentBoxWidth',MB{handles.id}.parsAlignHybes.alignmentBoxWidth,...
            % ,... 
            fidFrames = cat(3,dataOut.fiducialFrames{:});
            if showFigs
                figure(1); clf; Ncolor(3/size(fidFrames,3)*IncreaseContrast(fidFrames,'high',.999));
                pause(.1); title(['Fid images fov ',num2str(f)]); 
                figTile = figure(2); clf; TileImageStack(IncreaseContrast(fidFrames,'high',.999),...
                    'tileLabels',chnNames);
                pause(.1);
                SaveFigure(figTile,'name',['TileFigFOV',num2str(f)],'formats',{'png'},'overwrite',true);
            end
            
            disp(['processing FOV ',num2str(f)]); 
            [fidMaps{f},datMaps{f},fidFlatFrames,datFlatFrames] = ...
                ChrTracer3_SaveAlignedData(chnTable,dataOut.regData,...
                'memMapData',showFigs,...
                'dataFolder',MB{handles.id}.parsFindFiles.chnDataFolder,...
                'saveFolder',MB{handles.id}.parsFindFiles.chnSaveFolder,...
                'overwrite',MB{handles.id}.parsAlignHybes.overwrite,...
                'parameters',pars);
            if showFigs
             fidTileFig = figure(2); clf; TileImageStack(IncreaseContrast(fidFlatFrames,'high',.999),...
                'tileLabels',chnNames);
             SaveFigure(fidTileFig,'name',['FidTile_FOV',num2str(f,'%02d')],'formats',{'png'},'overwrite',true);
             datTileFig = figure(3); clf; TileImageStack(IncreaseContrast(datFlatFrames,'high',.999),...
                'tileLabels',chnNameList);
             SaveFigure(datTileFig,'name',['DatTile_FOV',num2str(f,'%02d')],'formats',{'png'},'overwrite',true);
            end
            end
        catch er
            warning(er.getReport);
        end
        
    end
    disp('file saving complete');

elseif strcmp(currStep,'AlignTiles')
    %% Step 3: Mosaic of fiducial data
    hyb1Tag = 'h001'; % should become an input parameter
    chnSaveFolder = [MB{handles.id}.parsFindFiles.chnSaveFolder,filesep]; %  (short hand); 
    MB{handles.id}.parsAlignTiles.saveCorrFolder = ...
        SetFigureSavePath([chnSaveFolder,'CorrAlign\'],'makeDir',true);
    showFinal = MB{handles.id}.parsAlignTiles.showFinal;
    downsampleFinal = MB{handles.id}.parsAlignTiles.showFinal;
    hyb1Files = cellstr(ls( [chnSaveFolder,'fov*_',hyb1Tag,'_fid.dax']));
    daxFiles = strcat(chnSaveFolder,hyb1Files);
    % load the daxfiles.  Any transpose rotate issues are corrected here as well  
    [imageTiles,stageXY] = DaxToImageTiles(daxFiles,'parameters',MB{handles.id}.parsAlignTiles);
    % register the images   
    [mosaicImage,boxCoords,xyShifts] = MosaicViewerRender(imageTiles,stageXY,'parameters',MB{handles.id}.parsAlignTiles);
    % render the mosiac
    if showFinal
        if MB{handles.id}.parsFindFiles.verbose
            disp('rendering final mosaic image...');
        end
        figure(2); clf;
        mosaicIm = imresize(mosaicImage,1/downsampleFinal);
        imagesc(IncreaseContrast(mosaicIm,'low',.6,'high',.999));
        colormap(gray);
    end
    % mosaic1 holds the original shift data
    MB{handles.id}.mosaic1.xyShifts = xyShifts;
    MB{handles.id}.mosaic1.mosaicImage = mosaicImage;
    % the export data is stored in the top level structure
    MB{handles.id}.imageTiles = imageTiles;
    MB{handles.id}.boxCoords = boxCoords;
    MB{handles.id}.tileNames = daxFiles; % keep the names
elseif strcmp(currStep,'AdjustTiles')
    imageTiles = MB{handles.id}.imageTiles;
    stageXY = MB{handles.id}.boxCoords(:,[1,3]);
    newStageXY = MosaicViewerGUI(imageTiles(:,end),stageXY,'uiwait',true,'names',chnNames);
    [h_i,w_i] = size(imageTiles{1});
    boxCoords = [newStageXY(:,1),newStageXY(:,1)+w_i,newStageXY(:,2),newStageXY(:,2)+h_i];
    MB{handles.id}.boxCoords = boxCoords;
    
elseif strcmp(currStep,'TileAllHybes')
    %% all channel mosaic
    % this should use the reference frame and shifts computed from the
    % first channel.
    % because mpars_all now contains non-empty pars.shiftsXY, correlation drift
    % corrections will not be recomputed.
    justNames = MB{handles.id}.parsTileAllHybes.justNames;
    displayMosaic = MB{handles.id}.parsTileAllHybes.displayMosaic;
    downsample = MB{handles.id}.parsTileAllHybes.downsample;
    nFOV = size(MB{handles.id}.boxCoords,1);
    nDatasPerFid = length(cellstr(ls([chnSaveFolder,'fov001','_h',num2str(1,'%03d'),'_data*.dax'])));
    nFids = size(allNames,1);
    % size(allNames{1},2)
    % nChns = nDatasPerFid * nFids; 
    nChns = length(chnNameList);
    mosaics_all  = cell(1,nChns);
    tiles_all = cell(nFOV,nChns);
    boxCoords_all = cell(1,nChns);
    daxNames_all = cell(nFOV,nChns);
    r=0;
    for h=1:nFids
        for d=1:nDatasPerFid
            r=r+1;
            disp(['processing mosaic for ',chnNameList{r}]);
            daxFiles_r = strcat(chnSaveFolder,cellstr(ls([chnSaveFolder,'fov*','_h',num2str(h,'%03d'),'_data',num2str(d),'.dax'])));
            daxNames_all(:,r) = daxFiles_r;
            if ~justNames
                [tiles_all(:,r),stageXY] = DaxToImageTiles(daxFiles_r,'parameters',MB{handles.id}.parsAlignTiles);
                 % this is just for display purposes, the full resolution data is in the tiles_all variable 
                 % in future should be able to turn off this request. 
                 mosaics_all{r} = MosaicViewerRender(tiles_all(:,r),stageXY,...
                     'parameters',MB{handles.id}.parsAlignTiles,'corrAlign',false,...
                     'xyShifts',MB{handles.id}.mosaic1.xyShifts,'downsample',downsample,...
                     'reZeroUL',true,'interactive',false,'showplot',false,...
                     'flatten',MB{handles.id}.parsTileAllHybes.flatten);
            end
        end
    end
    
    % -----add fiducial image from Hyb 1 at end of stack---------------
    daxFiles_r = strcat(chnSaveFolder,cellstr(ls([chnSaveFolder,'fov*','_h',num2str(1,'%03d'),'_fid.dax'])));
    daxNames_all(:,end+1) = daxFiles_r;
    if ~justNames
        [tiles_all(:,end+1),stageXY] = DaxToImageTiles(daxFiles_r,'parameters',MB{handles.id}.parsAlignTiles);
         % this is just for display purposes, the full resolution data is in the tiles_all variable 
         % in future should be able to turn off this request. 
         mosaics_all{end+1} = MosaicViewerRender(tiles_all(:,r),stageXY,...
             'parameters',MB{handles.id}.parsAlignTiles,'corrAlign',false,...
             'xyShifts',MB{handles.id}.mosaic1.xyShifts,'downsample',downsample,...
             'reZeroUL',true,'interactive',false,'showplot',false);
    end
       
    disp('finished loading RNA mosaics');
    % Display Mosaic
    mosaicImages = cat(3,mosaics_all{:}); % a downsampled view
    % set background to equal the mode
    if displayMosaic && ~justNames
        for r=1:nChns
            temp = mosaicImages(:,:,r);
            bkd = mode(nonzeros(temp(:)));
            temp(temp<=bkd) = bkd;
            mosaicImages(:,:,r) = temp;
        end
        sel = 2:min(5,nChns);% [28,6,23,10,4,18,29,8];
        disp(chnNameList(sel));
        contrastImA = IncreaseContrast(.2*mosaicImages(:,:,sel),'low',.85,'high',.99995);
        figure(2); clf; Ncolor(contrastImA);
        axis image;
    else
        mosaicImages = [];
    end
    if ~justNames
        MB{handles.id}.mosaicImages = mosaicImages; % don't overwrite if we didn't load new images
        MB{handles.id}.imageTiles = tiles_all; % not sure we need to save these 
    end
    MB{handles.id}.tileNames =daxNames_all;
       
elseif strcmp(currStep,'DisplayMosaic')
    %% Just viewing fun of Display Mosaic

    stageXY = MB{handles.id}.boxCoords(:,[1,3]);
    MosaicViewerGUI(MB{handles.id}.imageTiles,stageXY,...
         'names',chnNameList,'uiwait',true);
    
elseif strcmp(currStep,'SaveMosaic')
    %% Save Mosaic 
    % new approach:
    % just save a small text table of the updated tile locations, names
    % and paths to original and or compressed data
    saveFolder = MB{handles.id}.parsSaveMosaic.saveFolder; 
    saveRoot = MB{handles.id}.parsSaveMosaic.saveName;
    boxCoords = MB{handles.id}.boxCoords;
    tileNames = MB{handles.id}.tileNames;
    tileFile = [saveFolder,saveRoot,'_tileFileNames.txt'];
    chnFile = [saveFolder,saveRoot,'_chnNames.txt'];
    tilePositionsFile = [saveFolder,saveRoot,'_tilePositions.csv'];
    tileMatFile = [saveFolder,saveRoot,'.mat'];
    % write tile names
    writeFile = OverwriteDlgBox(tileFile);
    if writeFile
        fid = fopen(tileFile,'w+'); 
        fprintf(fid,'%s\r\n',tileNames{:}); 
        fclose(fid);
        result1 = ['saving ',tileFile]; 
        disp(result1);  set(handles.TextDir,'String',result1); guidata(hObject, handles);
    end
    % write table
    xMin = boxCoords(:,1); xMax=boxCoords(:,2); yMin = boxCoords(:,3); yMax = boxCoords(:,4);
    tilePosTable = table(xMin,xMax,yMin,yMax);
    writeFile = OverwriteDlgBox(tilePositionsFile);
    if writeFile
        writetable(tilePosTable,tilePositionsFile);
        result1 = ['saving ',tilePositionsFile]; 
        disp(result1);  set(handles.TextDir,'String',result1); guidata(hObject, handles);
    end  
    % write channel names
    writeFile = OverwriteDlgBox(chnFile);
    if writeFile
        fid = fopen(chnFile,'w+'); 
        fprintf(fid,'%s\r\n',chnNameList{:}); 
        fclose(fid);
        result1 = ['saving ',chnFile]; 
        disp(result1);  set(handles.TextDir,'String',result1); guidata(hObject, handles);
    end
    % write mat file
    mosaicData.chnNames = chnNameList;
    mosaicData.tileNames = tileNames;
    mosaicData.boxCoords = boxCoords; %#ok<STRNU>
    writeFile = OverwriteDlgBox(tileMatFile);
    if writeFile
        result1 = ['saving ',tileMatFile]; 
        disp(result1);  set(handles.TextDir,'String',result1); guidata(hObject, handles);
        save(tileMatFile,'mosaicData');
    end
end






% --- Executes on button press in ButtonNext.
function ButtonNext_Callback(hObject, eventdata, handles)
global MB;
step = find(strcmp(MB{handles.id}.stepNames, MB{handles.id}.currStepName));
directions = MB{handles.id}.Directions;
numSteps = length(MB{handles.id}.stepNames);
if step < numSteps
     step = step + 1;
     MB{handles.id}.currStepName = MB{handles.id}.stepNames{step};
     set(handles.TextDir,'String',directions{step}); guidata(hObject, handles);
else
    warning('Already reached last step'); 
end


% --- Executes on button press in ButtonBack.
function ButtonBack_Callback(hObject, eventdata, handles) %#ok<*INUSL>
global MB
step = find(strcmp(MB{handles.id}.stepNames, MB{handles.id}.currStepName));
directions = MB{handles.id}.Directions;
if step > 1
     step = step - 1;
     MB{handles.id}.currStepName = MB{handles.id}.stepNames{step};
     set(handles.TextDir,'String',directions{step}); guidata(hObject, handles);
else
    warning('Currently on step 1'); 
end
