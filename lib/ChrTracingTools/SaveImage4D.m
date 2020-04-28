function SaveImage4D(datMat,fullSaveName,varargin)
% ----------------------------------------------------------------------- %
% save a 4D image as a flat binary 16 bit file 
% with accompany info file
% 
% ----------------------------------------------------------------------- %
% Alistair Boettiger (boettiger@stanford.edu)
% CC BY Aug 8 2017
% ----------------------------------------------------------------------- %

% default parameters
defaults = cell(0,3);
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'description','string',''};
pars = ParseVariableArguments(varargin,defaults,mfilename);

if ~contains(fullSaveName,'.i4d')
    error('SaveImage4D file name must be of type ".i4d"');
end

% setup info file
[nRows,nCols,nStks,numColors] = size(datMat);
datInfo.total_rows = nRows;
datInfo.total_columns = nCols;
datInfo.height_zstack = nStks;
datInfo.numColors = numColors; 
datInfo.description = pars.description;

% a little information needed for saving 
sizeData = datInfo.total_rows*datInfo.total_columns*datInfo.height_zstack*datInfo.numColors;
binaryFormat = 'l';

% write binary file
fid = fopen(fullSaveName,'w+');
fwrite(fid, uint16(reshape(datMat,sizeData,1)), 'int16',binaryFormat);
fclose(fid);

% write info file
infoName = regexprep(fullSaveName,'.i4d','.inf');
SaveStructAsText(datInfo,infoName,'verbose',false); 

if pars.verbose
   disp(['wrote ',fullSaveName]);  
end

