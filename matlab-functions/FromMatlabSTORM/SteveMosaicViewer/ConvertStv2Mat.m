function ConvertStv2Mat(varargin)
% convert all .stv files to .mat files so they can be read.
% if files have already been converted, do nothing. 

global pythonPath matlabStormPath stvfile

if length(varargin) == 0
    mosaicFile = stvfile;
else
    mosaicFile = varargin{1};
end

disp(['loading ',mosaicFile]);

mosaic_to_matlab =[matlabStormPath,'\GUIs\Library\STORMrender\mosaic_to_matlab.py'];
runLocation = ''; %  Local
% runLocation = ' &'; % External

mosaicFolder = [fileparts(mosaicFile),'/'];
mosaicFolder = regexprep(mosaicFolder,'\','/'); % python prefers linux slashes 
mosaicName = dir([mosaicFolder,'/','*.msc']);

stvFiles = dir([mosaicFolder,'*.stv']);
matFiles = dir([mosaicFolder,'*.mat']); 

if isempty(stvFiles)
    warning('no .stv files found in folder');
    disp(mosaicFolder);
end

if length(matFiles) < length(stvFiles);
    stv1 =mosaicName.name; %  'mosaic_352.stv';
    dir(mosaic_to_matlab)

    % 
    ccall = ['C:\Anaconda2\','python.exe ',...
        mosaic_to_matlab,...
        ' ',mosaicFolder,'/',stv1,runLocation];

    disp(ccall);
    system(ccall);
end

