function [movieN,infoPerCell,movieFid,parameters] = LoadMERFISHdax(foundFiles,varargin)
% inputs - foundFiles
% outputs
% 

defaults = cell(0,3);
defaults(end+1,:) = {'imageTag', 'string', 'STORM'}; % Base tag for all images
defaults(end+1,:) = {'saveDataPath', 'string', ''}; % Base tag for all images
defaults(end+1,:) = {'bitFrames', 'integer', 1}; %   vector of the frames in the movie containing images of RNA bits  
defaults(end+1,:) = {'fidFrames', 'integer', []}; % 
defaults(end+1,:) = {'overwriteDax', 'boolean', false}; % 
defaults(end+1,:) = {'verbose', 'boolean', true}; % 
defaults(end+1,:) = {'maxCells', 'integer', inf}; % 
parameters = ParseVariableArguments(varargin, defaults, mfilename);

%% Main Function

% x and y dim of raw dax movie (could be determined automatically): 
infoFile = ReadInfoFile(foundFiles(1).filePath);
yDim = infoFile.frame_dimensions(2);
xDim = infoFile.frame_dimensions(1); 
numFrames = infoFile.number_of_frames; % not used?

dataPath = [fileparts(foundFiles(1).filePath),filesep];
% bitMovies = strcmp({foundFiles.movieType},parameters.imageTag);
numHybes = length(unique([ foundFiles.hybNum]));
numCells = min( length(unique([ foundFiles.cellNum])), parameters.maxCells);
numChns = length(parameters.bitFrames);
numBits = numHybes*numChns;
if isfield(foundFiles,'isFiducial');
    isFiducial = [foundFiles.isFiducial];
else
    isFiducial = false(size(foundFiles));
end

hybeNamesAll = {foundFiles.name}';

parameters.numHybes = numHybes;
parameters.numBits = numBits;
parameters.numChns = numChns;

if ~isempty(parameters.saveDataPath) && ~parameters.overwriteDax
    hasData = dir([parameters.saveDataPath,'*Balanced*dax']);
    startCell = length(hasData)+1;
    if parameters.verbose && startCell > 1
        disp(['found existing data, starting from cell ',num2str(startCell)]);
    end
else
    startCell = 1;
end

% one dax file per cell (this is the most natural)
movieN = zeros(yDim,xDim,numCells,numBits,'uint16');  % raw movie per cell
movieFid = zeros(yDim,xDim,numCells,numHybes,'uint16'); % fiducial per cell 
infoPerCell = repmat(infoFile,1,numCells);
for c=startCell:numCells   
    try
        if parameters.verbose
            disp(['loading data for cell ',num2str(c),'...']);
        end
        for h=1:numHybes
            fileName = [dataPath,...
                 hybeNamesAll{~isFiducial &...
                 [foundFiles.cellNum]==c-1 &...
                 [foundFiles.hybNum]==h-1}];
            [daxMov,infoN] = ReadDax(fileName,'verbose',false);   %  07/02/15, removed from list of logicals on names: bitMovies &   
            if h==1; infoPerCell(c) = infoN; end % save info file

            % store movie in a rationally indexed array
            f = numChns*(h-1)+1:numChns*h;
            movieN(:,:,c,f) = daxMov(:,:,parameters.bitFrames);

            % read corresponding fiducial movie (if necessary).  This could be gotten from a specific frame for the other movie 
            if ~isempty(parameters.fidFrames)
                 movieFid(:,:,c,h) = daxMov(:,:,parameters.fidFrames);
            end
            if sum(isFiducial)>0
                fidName = [dataPath,...
                                  hybeNamesAll{isFiducial &...
                                  [foundFiles.cellNum]==c-1 &...
                                  [foundFiles.hybNum]==h-1}];                         
                movieFid(:,:,c,h) = ReadDax(fidName,'verbose',false);   %  07/02/15, removed from list of logicals on names: bitMovies &
            end
        end
    catch er
        warning(er.getReport);
        warning(['error reading data from cell ',num2str(c),' skipping this cell.']); 
    end
end

junkOnLens = uint16(mean(reshape(movieFid,yDim,xDim,numCells*numHybes),3));
junkOnLens = junkOnLens - min(junkOnLens(:));
% figure(1); clf; imagesc(junkOnLens); colorbar; colormap(gray);
movieFid = movieFid-repmat(junkOnLens,1,1,numCells,numHybes);