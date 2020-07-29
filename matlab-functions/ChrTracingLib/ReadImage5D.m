function [datMat,datInfo,datTable] = ReadImage5D(fileName, varargin)
% ----------------------------------------------------------------------- %
% Read in a 5D binary image file
% 
% see also WriteImage5D, ReadImage4D, SaveImage4D, PlotProjection4D
% ----------------------------------------------------------------------- %
% Alistair Boettiger (boettiger@stanford.edu)
% CC BY Aug 8 2017
% 
% ----------------------------------------------------------------------- %

% default parameters
defaults = cell(0,3);
% defaults for ReadImage4D
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'showPlots','boolean',false};
defaults(end+1,:) = {'loadTable','boolean',false};
defaults(end+1,:) = {'spots','integer',inf}; % list of spots in ascending order
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
pars = ParseVariableArguments(varargin,defaults,mfilename);


% Read info file
datInfo = ReadTextToStruct(fileName); 

% Set filepointer
binaryFormat = 'l'; % hard coded;
bytes_per_unit = 2; % for uint16 (1 for uint8, 4 for single etc)
fid = fopen(fileName,'r',binaryFormat);
if fid < 0
    error(['Invalid file: ' infoFile.localPath fileName]);
end



if isinf(pars.spots)
    % Load data 
    sizeData = datInfo.total_rows*datInfo.total_columns*datInfo.height_zstack*datInfo.numClrs*datInfo.numSpts;% *bytes_per_unit;
    fseek(fid,0,'bof'); % bits/(bytes per bit) 
    movie = fread(fid, sizeData, '*uint16',0, binaryFormat);
    fclose(fid);
    % reshape data into correct form indicated in info file
    datMat = reshape(movie,[datInfo.total_rows, datInfo.total_columns, datInfo.height_zstack, datInfo.numClrs, datInfo.numSpts]);
    % Optional display progress
    if pars.verbose
        disp(['loaded ',fileName]);
    end
else
     sizeSpot = datInfo.total_rows*datInfo.total_columns*datInfo.height_zstack*datInfo.numClrs*bytes_per_unit;
     numSpots = length(pars.spots);
     movies = cell(numSpots,1);
     fseek(fid,sizeSpot*(pars.spots(1)-1),'bof');
     ftell(fid)
     temp = fread(fid, sizeSpot/bytes_per_unit, '*uint16',0, binaryFormat); % *uint16
     movies{1} = reshape(temp,[datInfo.total_rows, datInfo.total_columns, datInfo.height_zstack, datInfo.numClrs]); 
     for s=2:numSpots
          % if the spot is not consecutive, we fseek forward from current
          % postion the appropriate number of steps.  
          stp = pars.spots(s) - pars.spots(s-1) -1; 
          fseek(fid,sizeSpot*stp,'cof');  
          temp = fread(fid, sizeSpot/bytes_per_unit, '*uint16',0, binaryFormat); 
          movies{s} = reshape(temp,[datInfo.total_rows, datInfo.total_columns, datInfo.height_zstack, datInfo.numClrs]);    
     end
     fclose(fid);
     datMat = cat(5,movies{:});
     % datInfo.numSpts
    
end


if pars.loadTable
    % needs to be updated to i5d style tables
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
    %% new version
    figure(); PlotProjection4D(datMat,'projection','xy','fits',datTable,'parameters',pars);
    figure(); PlotProjection4D(datMat,'projection','xz','fits',datTable,'parameters',pars);
%     %% old version
%     figure(); 
%     datMatXY = squeeze(max(datMat,[],3));
%     TileImageStack(datMatXY);
%     figure();
%     datMatXZ = squeeze(max(permute(datMat,[3,2,1,4]),[],3));
%     TileImageStack(datMatXZ);
end