function images = CombineHybImages(imageData, cellIDs, varargin)
% ------------------------------------------------------------------------
% images = CombineHybImages(imageData, cellIDs, varargin)
%--------------------------------------------------------------------------
% Necessary Inputs
% imageData/A structure array with an element for each image used to create
%   the elements in words.  See CreateWordsStructure for information on
%   field names. 
%--------------------------------------------------------------------------
% Outputs
% images/A NxMXLstructure containing information from the report
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% Jeffrey Moffitt 
% jeffmoffitt@gmail.com
% October 23, 2014
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'colormap', 'function', @jet};
defaults(end+1,:) = {'verbose', 'boolean', true};
% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);        

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if (nargin < 2 || ...
        ~isstruct(imageData) || ...
        ~isempty(setdiff(fields(CreateImageDataStructure(0)), fields(imageData)) ))
    error('matlabFunctions:invalidArguments', 'Invalid imageData structures.');
end

% -------------------------------------------------------------------------
% Determine cellIDs to use and sort in order from lowest to highest
% -------------------------------------------------------------------------
foundIDs = unique([imageData.cellNum]);
if isempty(cellIDs)
    cellIDs = foundIDs;
else
    cellIDs = intersect(foundIDs, cellIDs);
end
cellIDs = sort(cellIDs);

% -------------------------------------------------------------------------
% Determine cellIDs to use and sort in order from lowest to highest
% -------------------------------------------------------------------------
images = uint16(zeros(length(cellIDs), imageData(1).imageW, imageData(1).imageH, 3));
for i=1:length(cellIDs)
    if parameters.verbose
        display(['Creating combined hyb image for cell ' num2str(cellIDs(i))]);
    end
    
    tic;
    % Extract local image data and sort by hyb number
    localImageData = imageData([imageData.cellNum]==cellIDs(i));
    [~, sind] = sort([localImageData.hybNum]);
    localImageData = localImageData(sind);
    
    alignedIm = zeros(localImageData(1).imageH, localImageData(1).imageW, length(localImageData));
    for h=1:length(localImageData)
        dax = max(ReadDax(localImageData(h).infFilePath,'endFrame',1, 'verbose', false),[],3); % max project the first frame
        [H,W] = size(dax); 
        tformInv = fliptform(localImageData(h).tform); 
        alignedDax = imtransform(dax,tformInv,...
                        'XYScale',1,'XData',[1 W],'YData',[1 H]);
        alignedIm(:,:,h) = alignedDax;
    end
    
    % Combine image and add color
    images(i,:,:,:) = Ncolor(uint16(alignedIm), parameters.colormap(length(localImageData)));
    if parameters.verbose
        display(['...finished in ' num2str(toc) ' s']);
    end
end