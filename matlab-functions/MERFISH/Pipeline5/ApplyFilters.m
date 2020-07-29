function [moviesOut,infoOut,parameters] = ApplyFilters(movieNW,infoNW,varargin)

% parameters
psf = fspecial('gaussian',10,2);

defaults = cell(0,3);
defaults(end+1,:) = {'upsample', 'positive', 3}; % for deconvolution
defaults(end+1,:) = {'psf', 'array', psf}; % for deconvolution
defaults(end+1,:) = {'numIters', 'positive', 5}; % for deconvolution
defaults(end+1,:) = {'numChns', 'positive', 1}; % for deconvolution
defaults(end+1,:) = {'figVis', 'string', 'off'}; % 
defaults(end+1,:) = {'closeOnSave', 'boolean', true}; % 
defaults(end+1,:) = {'saveFigures', 'boolean', true}; % 
defaults(end+1,:) = {'showImages', 'boolean', false}; % 
defaults(end+1,:) = {'overwriteDax', 'boolean', false}; % 
defaults(end+1,:) = {'laplaceFilter', 'boolean', true}; % 
defaults(end+1,:) = {'deconvolve', 'boolean', true}; % 
defaults(end+1,:) = {'verbose', 'boolean', true}; 
parameters = ParseVariableArguments(varargin, defaults, mfilename);

%% Main Function

diskFilt2 = fspecial('laplacian'); % (hard coded intentionally);

[nRows,nCols,numCells,numBits] = size(movieNW);

movieNWL = movieNW; % warped laplace
movieNWLD = imresize(movieNW,parameters.upsample); % warped, laplace, deconvovled
infoNWL = infoNW;
infoNWLD = infoNW;

origName = {infoNW.localName}';
daxRawName = regexprep(origName,{'ImagesCell','.inf'},{'ImagesRawCell','.dax'}); 
daxLaplaceName = regexprep(origName,{'ImagesCell','.inf'},{'ImagesLaplaceCell','.dax'}); 
daxFinalName = regexprep(origName,{'ImagesCell','.inf'},{'ImagesBalancedCell','.dax'}); 
newDataPath = infoNW.localPath;

for c=1:numCells 
    file1Exists = ~isempty( dir([newDataPath,daxRawName{c}]) );
    file2Exists = ~isempty( dir([newDataPath,daxLaplaceName{c}]) );
    file3Exists = ~isempty( dir([newDataPath,daxFinalName{c}]) );
    if parameters.verbose 
        disp(['filtering and saving data for cell ',num2str(c)]);
    end
    testData = movieNW(:,:,c,1);
    notEmpty = sum(testData(:)) > 0;
    
    if (~(file1Exists && file2Exists && file3Exists) || parameters.overwriteDax) && notEmpty    
        for b=1:numBits         
            % Apply Laplacae filter
            if parameters.laplaceFilter
                temp = imfilter(double(movieNW(:,:,c,b)),diskFilt2,'symmetric');
                temp = 2^16-temp; temp = temp-.9999*mode(temp(:)); temp = 2^16*temp/max(temp(:));  % invert, set background at zero;
            else
                temp = double(movieNW(:,:,c,b));
            end
            movieNWL(:,:,c,b) = uint16(temp); 
            
            % Upsample and deconvolve
            if parameters.deconvolve
                temp = imresize(temp,parameters.upsample);
                temp = deconvlucy(temp,parameters.psf,parameters.numIters); % deconvolution to sharpen
                temp = temp - min(temp(:)); temp = 2^16*temp/quantile(temp(:),.9999);  % histogram equalize
                movieNWLD(:,:,c,b) = uint16(temp);
            else
                movieNWLD(:,:,c,b) = uint16(temp);
            end
        end
        
         daxOut = squeeze(movieNW(:,:,c,:)); 
         infoNWR = infoNW(c); % info file for raw data file
         infoNWR.localName = regexprep(daxRawName{c},'.dax','.inf'); 
         infoNWR.number_of_frames = size(daxOut,3); 
         WriteDAXFiles(daxOut,infoNWR,'verbose',false); 
        
        % Write dax files  
        if parameters.laplaceFilter
            daxOutLaplace = squeeze(movieNWL(:,:,c,:)); 
            infoNWL(c) = infoNW(c); % info file for raw data file
            infoNWL(c).localName = regexprep(daxLaplaceName{c},'.dax','.inf'); 
            infoNWL(c).number_of_frames  = size(daxOutLaplace,3); 
            WriteDAXFiles(daxOutLaplace,infoNWL(c),'verbose',false); 
        end
        
        if parameters.deconvolve
            daxOutNormalized = squeeze(movieNWLD(:,:,c,:)); 
            infoNWLD(c) = infoNW(c); % info file for warped data  
            [yDim,xDim,~] = size(daxOutNormalized);
            infoNWLD(c).frame_dimensions = [xDim,yDim];
            infoNWLD(c).frame_size = infoNWLD(c).frame_dimensions(1)*infoNWLD(c).frame_dimensions(2);
            infoNWLD(c).number_of_frames  = size(daxOutNormalized,3); 
            infoNWLD(c).notes = ['upsampled by ',num2str(parameters.upsample)];
            infoNWLD(c).localName = regexprep(daxFinalName{c},'.dax','.inf'); 
            WriteDAXFiles(daxOutNormalized,infoNWLD(c),'verbose',false); 
        end
                   
    else
        if parameters.verbose
            disp(['Skipping warp for cell ',num2str(c),' Using existing file...']);
        end
    end
end


moviesOut{1} = cell(0,1);
infoOut{1} = cell(0,1);
if parameters.laplaceFilter
    moviesOut{end+1} = movieNWL;
    infoOut{end+1} = infoNWL;
end
if parameters.deconvolve
   moviesOut{end+1} = movieNWLD;
   infoOut{end+1} = infoNWLD;
end


