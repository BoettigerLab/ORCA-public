function [datMat,datInfo] = ReadImage4D(fileName, varargin)
% ----------------------------------------------------------------------- %
% Read in a 4D binary image file
% 
% ----------------------------------------------------------------------- %
% Alistair Boettiger (boettiger@stanford.edu)
% CC BY Aug 8 2017
% ----------------------------------------------------------------------- %

% default parameters
defaults = cell(0,3);
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'showPlots','boolean',false};
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

if pars.showPlots
    figure(); 
    datMatXY = squeeze(max(datMat,[],3));
    TileImageStack(datMatXY);
    figure();
    datMatXZ = squeeze(max(permute(datMat,[3,2,1,4]),[],3));
    TileImageStack(datMatXZ);
end