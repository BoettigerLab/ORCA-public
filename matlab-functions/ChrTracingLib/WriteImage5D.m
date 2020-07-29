function WriteImage5D(im5D,fullSaveName,varargin)
% ----------------------------------------------------------------------- %
% save a 5D image as a flat binary 16 bit file 
% with accompany info file
% 
% ----------------------------------------------------------------------- %
% Alistair Boettiger (boettiger@stanford.edu)
% CC BY April 2020
% ----------------------------------------------------------------------- %


%%

% default parameters
defaults = cell(0,3);
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'description','string',''};
defaults(end+1,:) = {'fov','integer',0};
defaults(end+1,:) = {'lociX','positive',0};
defaults(end+1,:) = {'lociY','positive',0};
%pars = ParseVariableArguments([],defaults,mfilename);
pars = ParseVariableArguments(varargin,defaults,mfilename);

if ~contains(fullSaveName,'.i5d')
    error('WriteImage5D file name must be of type ".i5d"');
end

% setup info file
[nRows,nCols,nStks,numClrs,numSpts] = size(im5D);
datInfo.encoding = 'little endian'; % these are hardcoded options, we merely report them for portability  
datInfo.dataType = 'uint16'; % these are hardcoded options, we merely report them for portability
datInfo.total_rows = nRows;
datInfo.total_columns = nCols;
datInfo.height_zstack = nStks;
datInfo.numClrs = numClrs; 
datInfo.numSpts = numSpts; 
datInfo.fov = pars.fov;
datInfo.lociX = pars.lociX;
datInfo.lociY = pars.lociY;
datInfo.description = pars.description;

% a little information needed for saving 
sizeData = datInfo.total_rows*datInfo.total_columns*datInfo.height_zstack*datInfo.numClrs*datInfo.numSpts;
binaryFormat = 'l';

 % write binary file
 fid = fopen(fullSaveName,'w+');
 fwrite(fid, uint16(reshape(im5D,sizeData,1)), 'int16',binaryFormat);
 fclose(fid);

% write info file
infoName = regexprep(fullSaveName,'.i5d','.inf');
SaveStructAsText(datInfo,infoName,'verbose',false); 

if pars.verbose
   disp(['wrote ',fullSaveName]);  
end

