function [datMat,datInfo,datTable] = ReadImage4D(fileName, varargin)
% ----------------------------------------------------------------------- %
% Read in a 4D binary image file
% 
% see also SaveImage4D, PlotProjection4D
% ----------------------------------------------------------------------- %
% Alistair Boettiger (boettiger@stanford.edu)
% CC BY Aug 8 2017
% ----------------------------------------------------------------------- %

% default parameters
defaults = cell(0,3);
% defaults for ReadImage4D
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'showPlots','boolean',false};
defaults(end+1,:) = {'loadTable','boolean',false};
% defaults for PlotProjection4D
defaults(end+1,:) = {'autoContrast', 'boolean', true}; 
% defaults for image tile
defaults(end+1,:) = {'mode',{'subplots','single'},'single'};  % render as subplots (slow and flexible) or as a single image  
defaults(end+1,:) = {'numRows', 'positive', 4}; 
defaults(end+1,:) = {'gap', 'nonnegative', 1}; 
defaults(end+1,:) = {'colorTiles', 'boolean', true}; 
defaults(end+1,:) = {'colormap', 'colormap', 'hsv'}; 
defaults(end+1,:) = {'showPanelNumber', 'boolean', true}; 
defaults(end+1,:) = {'fontSize', 'positive', 10}; 
defaults(end+1,:) = {'tileLabels','freeType',{}};
defaults(end+1,:) = {'fixBlack','boolean',true};
% defaults for plotting fits
defaults(end+1,:) = {'fits','freeType',[]};
defaults(end+1,:) = {'xyzNames','cell',{}}; % {'xc','yc','zc'}
defaults(end+1,:) = {'nmXYpix','positive',154};
defaults(end+1,:) = {'nmZpix','positive',100};
defaults(end+1,:) = {'figHandle','freeType',0};
defaults(end+1,:) = {'MarkerSize','positive',20};
pars = ParseVariableArguments(varargin,defaults,mfilename);

% Read info file
datInfo = ReadTextToStruct(fileName); 

% Set filepointer
fid = fopen(fileName);
if fid < 0
    error(['Invalid file: ' infoFile.localPath fileName]);
end

% Load data 
sizeData = datInfo.total_rows*datInfo.total_columns*datInfo.height_zstack*datInfo.numColors;
binaryFormat = 'l';
fseek(fid,0,'bof'); % bits/(bytes per bit) 
movie = fread(fid, sizeData, '*uint16', binaryFormat);
fclose(fid);

% reshape data into correct form indicated in info file
datMat = reshape(movie,[datInfo.total_rows, datInfo.total_columns, datInfo.height_zstack, datInfo.numColors]);

% Optional display progress
if pars.verbose
    disp(['loaded ',fileName]);
end


if pars.loadTable
    tableName = regexprep(fileName,'_AlignedData.i4d','_fits.csv');
    try
        datTable =  readtable(tableName); 
    catch
        datTable = [];
    end
else
    datTable = [];
end
if pars.showPlots 
    if pars.figHandle == 0
        f1 = figure();
        f2 = figure();
    else
        f1 = pars.figHandle(1);
        f2 = pars.figHandle(2);
    end
    
    %% new version
    figure(f1); PlotProjection4D(datMat,'projection','xy','fits',datTable,'parameters',pars);
    figure(f2); PlotProjection4D(datMat,'projection','xz','fits',datTable,'parameters',pars);
%     %% old version
%     figure(); 
%     datMatXY = squeeze(max(datMat,[],3));
%     TileImageStack(datMatXY);
%     figure();
%     datMatXZ = squeeze(max(permute(datMat,[3,2,1,4]),[],3));
%     TileImageStack(datMatXZ);
end