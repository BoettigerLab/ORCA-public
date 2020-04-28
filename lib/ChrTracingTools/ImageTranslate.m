function [B,RB] = ImageTranslate(varargin)
%
% NOTE: accelerated version of imtranslate
% 
%IMTRANSLATE Translate image.
%   B = IMTRANSLATE(A,TRANSLATION) translates image A by a translation
%   vector TRANSLATION. TRANSLATION is of the form [TX TY] for 2-D inputs,
%   and [TX TY TZ] for 3-D inputs. If TRANSLATION is a two-element vector
%   and A has more than two dimensions, a 2-D translation is applied to A
%   one plane at a time. TRANSLATION can be fractional.
%
%   [B,RB] = IMTRANSLATE(A,RA,TRANSLATION) translates the spatially
%   referenced image defined by A and RA by a translation vector
%   TRANSLATION. TRANSLATION is in the world coordinate system. The output
%   is a translated spatially referenced image defined by B and RB.
%
%   B = IMTRANSLATE(A,TRANSLATION,METHOD) translates image A, using the interpolation
%   method specified by METHOD. METHOD is a string that can have one of the
%   following values. The default value is enclosed in braces ({}).
%
%        'nearest'      Nearest neighbor interpolation
%
%        {'linear'}     Linear interpolation (Default)
%
%        'cubic'        Cubic interpolation. Note: This interpolation
%                       method can produce pixel values outside the original
%                       range.
%
%   [___] = IMTRANSLATE(___,Name, Value,...) translates the input image
%   using name-value pairs to control various aspects of the translation.
%
%   Parameters include:
%
%      'OutputView' - A string that defines the output world limits. Valid
%                     strings are 'same' and 'full'. When 'OutputView' is
%                     'same', the output world limits are the same as the
%                     input image. When 'OutputView' is 'full', the world
%                     limits are the bounding rectangle that includes both the
%                     input image and the translated output image.
%
%                     Default value: 'same'
%
%      'FillValues'   An array containing one or several fill values.
%                     Fill values are used for output pixels when the
%                     corresponding inverse transformed location in the
%                     input image is completely outside the input image
%                     boundaries.
%
%                     If A is 2-D then 'FillValues' must be a scalar. If A
%                     is 3-D and TRANSLATION is a 3-element vector, then
%                     'FillValues' must be a scalar. If A is N-D and
%                     TRANSLATION is a 2-element vector, then 'FillValues'
%                     may be either scalar or an array whose size matches
%                     dimensions 3 to N of A. For example, if A is a uint8
%                     RGB image that is 200-by-200-by-3, then 'FillValues'
%                     can be a scalar or a 3-by-1 array. In this RGB image
%                     example, possibilities for 'FillValues' include:
%
%                           0                 - fill with black
%                           [0;0;0]           - also fill with black
%                           255               - fill with white
%                           [255;255;255]     - also fill with white
%                           [0;0;255]         - fill with blue
%                           [255;255;0]       - fill with yellow
%
%                     If A is 4-D with size 200-by-200-by-3-by-10, then
%                     FillValues' can be a scalar or a 3-by-10 array.
%
%                     Default value: 0
%
%    Notes
%    -----
%
%    1) IMTRANSLATE is optimized for integrally valued TRANSLATION vectors.
%
%    2) When 'OutputView' is 'full' and translation is a fractional number
%    of pixels, the world limits of the output spatial referencing object
%    Rout are expanded to the nearest full pixel increment such that Rout
%    contains both the original and translated images at the same
%    resolution as the input image A. The additional image extent in each
%    is added on one side of the image, in the direction that the
%    translation vector points. For example, when translation is
%    fractional and positive in both X and Y, then the maximum of
%    XWorldLimits and YWorldLimits is expanded to enclose the 'full'
%    bounding rectangle at the resolution of the input image.
%
%    Class Support
%    -------------
%    A can be of any nonsparse, numeric class except uint64 and int64. A
%    can also be logical. TRANSLATION is a nonsparse, real valued numeric
%    vector.  The class of B is the same as the class of A. RA and RB are
%    spatial referencing objects of class imref2d or imref3d.
%
%   Example 1
%   ---------
%   Translate image I by 5.3 pixels in X and -10.1 pixels in Y.
%
%       I = imread('pout.tif');
%       J = imtranslate(I,[5.3, -10.1],'FillValues',255);
%       figure, imshow(I);
%       figure, imshow(J);
%
%   Example 2
%   ---------
%   Translate a 3-D MRI dataset. Use 'OutputView' to obtain full
%   translated image volume without clipping.
%
%       s = load('mri');
%       mriVolume = squeeze(s.D);
%       sizeIn = size(mriVolume);
%       hFigOriginal = figure;
%       hAxOriginal  = axes;
%       slice(double(mriVolume),sizeIn(2)/2,sizeIn(1)/2,sizeIn(3)/2);
%       grid on, shading interp, colormap gray
%
%       % Apply a translation in the X,Y direction
%       mriVolumeTranslated = imtranslate(mriVolume,[40,30,0],'OutputView','full');
%
%       % Visualize axial slice plane taken through center of volume
%       sliceIndex = round(sizeIn(3)/2);
%       axialSliceOriginal   = mriVolume(:,:,sliceIndex);
%       axialSliceTranslated = mriVolumeTranslated(:,:,sliceIndex);
%
%       imshowpair(axialSliceOriginal,axialSliceTranslated,'montage');
%
%   See also IMRESIZE, IMROTATE, IMWARP.

%   Copyright 2013 The MathWorks, Inc.

narginchk(2,inf);

[R_A, varargin] = preparseSpatialReferencingObjects(varargin{:});

[A,translation,method,outputView,fillValues] = parseInputs(varargin{:});

inputSpatialReferencingNotSpecified = isempty(R_A);
if inputSpatialReferencingNotSpecified
    if isa(R_A,'imref3d')
        R_A = imref3d(size(A));
    else
        R_A = imref2d(size(A));
    end
else
    % Check agreement of input spatial referencing object with input image.
    checkSpatialRefAgreementWithInputImage(A,R_A);
end

is2DProblem = isequal(numel(translation),2);

integrallyValuedTranslation = inputSpatialReferencingNotSpecified && all(mod(translation,1) == 0);

if integrallyValuedTranslation && isreal(A)
    % As a performance optimization, we treat non-spatially referenced,
    % real valued problems with integral translations in all dimensions as
    % a special case.
    if is2DProblem
        [B,RB] = translateIntegerShift2D(A,R_A,translation,method,outputView,fillValues);
    else
        [B,RB] = translateIntegerShift3D(A,R_A,translation,method,outputView,fillValues);
    end
    
else
    
    if is2DProblem
        [B,RB] = translate2D(A,R_A,translation,method,outputView,fillValues);
    else
        [B,RB] = translate3D(A,R_A,translation,method,outputView,fillValues);
    end
    
end

end

function [out,Rout] = translateIntegerShift2D(A,RA,translation,~,outputView,fillValues)

% This code path is a special case for non-spatially referenced cases in
% which the translation is integrally valued in all dimensions. We can
% avoid the computational cost of interpolation in these cases and form the
% output image with simple indexing.

Rout = computeOutputSpatialRef(RA,translation,outputView);

% Determine size of output image
inputSize = size(A);
if length(inputSize) < 3
    outputSize = Rout.ImageSize;
    numPlanes = 1;
else
    outputSize = [Rout.ImageSize, inputSize(3:end)];
    numPlanes = sum(inputSize(3:end));
end

% Pre-allocate output image to FillValue.
fillValueCastToOutputType = cast(fillValues,class(A));
if isscalar(fillValues)
    % This pre-allocation has to work with logical values as well as
    % numeric types.
    % out = cast(zeros(outputSize),class(A));  % Cast can be supper slow!
    out = zeros(outputSize,class(A)); % this is substantially faster
    out(:) = fillValueCastToOutputType;
else
    % out = cast(zeros(outputSize),class(A)); % Cast can be supper slow!
    out = zeros(outputSize,class(A));
    for i = 1:length(fillValues)
        out(:,:,i) = fillValueCastToOutputType(i);
    end
end

[XWorldBoundingSubscripts,YWorldBoundingSubscripts] = Rout.intrinsicToWorld([1 Rout.ImageSize(2)],...
    [1 Rout.ImageSize(1)]);

UWorld = XWorldBoundingSubscripts - translation(1);
VWorld = YWorldBoundingSubscripts - translation(2);

% Figure out whether bounding rectangle of reverse mapped pixel centers
% includes any in bounds locations in the source image.
intrinsicSourceBoundingRectangleU = [UWorld(1) UWorld(2) UWorld(2) UWorld(1)];
intrinsicSourceBoundingRectangleV = [VWorld(1) VWorld(1) VWorld(2) VWorld(2)];
% Contains is true inside the world limits and the boundary of the world
% limits.
locationsInSourceMapToDestination = any(RA.contains(intrinsicSourceBoundingRectangleU,...
    intrinsicSourceBoundingRectangleV));

% If there are locations in the source image that map into the destination,
% use indexing to form the output image. Otherwise, return all fill values.
if locationsInSourceMapToDestination
    
    % Clip reverse mapped boundaries to boundaries that live entirely
    % within A.
    UWorldClippedToBounds = [max(1,UWorld(1)), min(RA.ImageSize(2),UWorld(2))];
    VWorldClippedToBounds = [max(1,VWorld(1)), min(RA.ImageSize(1),VWorld(2))];
    
    % At this point we know the locations in source that map into valid
    % locations in the destination image. We want to forward map these into
    % corresponding subscripts in our output image.
    NonFillOutputLocX = UWorldClippedToBounds + translation(1);
    NonFillOutputLocY = VWorldClippedToBounds + translation(2);
    [outputR,outputC] = Rout.worldToSubscript(NonFillOutputLocX,NonFillOutputLocY);
    
    % Where the output locations map into valid, in-bounds locations in A,
    % assign the output values by simple indexing. No interpolation is
    % required since the translation is integrally valued.
    for i = 1:numPlanes
        out(outputR(1):outputR(2),outputC(1):outputC(2),i) = A(VWorldClippedToBounds(1):VWorldClippedToBounds(2),...
            UWorldClippedToBounds(1):UWorldClippedToBounds(2),i);
    end
    
end


end

function [out,Rout] = translateIntegerShift3D(A,RA,translation,~,outputView,fillValues)

% This code path is a special case for non-spatially referenced cases in
% which the translation is integrally valued in all dimensions. We can
% avoid the computational cost of interpolation in these cases and form the
% output image with simple indexing.

Rout = computeOutputSpatialRef(RA,translation,outputView);

fillValueCastToOutputType = cast(fillValues,class(A));
% This pre-allocation has to work with logical values as well as
% numeric types.
out = cast(zeros(Rout.ImageSize),class(A));
out(:) = fillValueCastToOutputType;

[XWorldBoundingSubscripts,YWorldBoundingSubscripts,ZWorldBoundingSubscripts] = Rout.intrinsicToWorld([1 Rout.ImageSize(2)],...
    [1 Rout.ImageSize(1)],[1 Rout.ImageSize(3)]);

UWorld = XWorldBoundingSubscripts - translation(1);
VWorld = YWorldBoundingSubscripts - translation(2);
WWorld = ZWorldBoundingSubscripts - translation(3);

% Figure out whether bounding rectangle of reverse mapped pixel centers
% includes any in bounds locations in the source image. In specific the
% bounding cube, wind CCW in the UV plane on the lower W face of the cube,
% then repeat winding CCW in the UV plane on the upper W face of the cube.
intrinsicSourceBoundingCubeU = [UWorld(1), UWorld(2), UWorld(2), UWorld(1),...
                                UWorld(1), UWorld(2), UWorld(2), UWorld(1)];
                            
intrinsicSourceBoundingCubeV = [VWorld(1), VWorld(1), VWorld(2), VWorld(2),...
                                VWorld(1), VWorld(1), VWorld(2), VWorld(2)];
                            
intrinsicSourceBoundingCubeW = [WWorld(1), WWorld(1), WWorld(1), WWorld(1),...
                                WWorld(2), WWorld(2), WWorld(2), WWorld(2)];
                            
% Contains is true inside the world limits and the boundary of the world
% limits.
locationsInSourceMapToDestination = any(RA.contains(intrinsicSourceBoundingCubeU,...
                                                    intrinsicSourceBoundingCubeV,...
                                                    intrinsicSourceBoundingCubeW));
                                                    
% If there are locations in the source image that map into the destination,
% use indexing to form the output image. Otherwise, return all fill values.
if locationsInSourceMapToDestination
    
    % Clip reverse mapped boundaries to boundaries that live entirely
    % within A.
    UWorldClippedToBounds = [max(1,UWorld(1)), min(RA.ImageSize(2),UWorld(2))];
    VWorldClippedToBounds = [max(1,VWorld(1)), min(RA.ImageSize(1),VWorld(2))];
    WWorldClippedToBounds = [max(1,WWorld(1)), min(RA.ImageSize(3),WWorld(2))];

    % At this point we know the locations in source that map into valid
    % locations in the destination image. We want to forward map these into
    % corresponding subscripts in our output image.
    NonFillOutputLocX = UWorldClippedToBounds + translation(1);
    NonFillOutputLocY = VWorldClippedToBounds + translation(2);
    NonFillOutputLocZ = WWorldClippedToBounds + translation(3);

    [outputR,outputC,outputP] = Rout.worldToSubscript(NonFillOutputLocX,NonFillOutputLocY,NonFillOutputLocZ);
    
    % Where the output locations map into valid, in-bounds locations in A,
    % assign the output values by simple indexing. No interpolation is
    % required since the translation is integrally valued.
    out(outputR(1):outputR(2),outputC(1):outputC(2),outputP(1):outputP(2)) = A(VWorldClippedToBounds(1):VWorldClippedToBounds(2),...
        UWorldClippedToBounds(1):UWorldClippedToBounds(2),WWorldClippedToBounds(1):WWorldClippedToBounds(2));
    
end

end

function [out,Rout] = translate2D(A,RA,translation,method,outputView,fillValues)

Rout = computeOutputSpatialRef(RA,translation,outputView);

% Compute spatially referenced case as a general 2-D affine transformation.
tform = affine2d([1 0 0; 0 1 0; translation(1:2) 1]);
[out,Rout] = imwarp(A,RA,tform,method,'OutputView',Rout,'fillValues',fillValues,...
    'SmoothEdges', true);

end


function [out,Rout] = translate3D(A,RA,translation,method,outputView,fillValues)

Rout = computeOutputSpatialRef(RA,translation,outputView);

% Compute spatially referenced case as a general 3-D affine transformation.
tform = affine3d([1 0 0 0; 0 1 0 0; 0 0 1 0; translation(1:3) 1]);
[out,Rout] = imwarp(A,RA,tform,method,'OutputView',Rout,'fillValues',fillValues,...
    'SmoothEdges', true);

end


function Rout = computeOutputSpatialRef(RA,translation,outputView)

if strcmp(outputView,'same')
    Rout = RA;
else
    
    % imtranslate(___,'OutputView','full');
    [XWorldLimitsOut,numColsOutput] = computeFullExtentAndGridSizePerDimension(RA.XWorldLimits,...
        RA.PixelExtentInWorldX,translation(1));
    
    [YWorldLimitsOut,numRowsOutput] = computeFullExtentAndGridSizePerDimension(RA.YWorldLimits,...
        RA.PixelExtentInWorldY,translation(2));
    
    if ~isa(RA,'imref3d')
        
        Rout = imref2d([numRowsOutput numColsOutput],XWorldLimitsOut,YWorldLimitsOut);
             
    else
           
        [ZWorldLimitsOut,numPlanesOutput] = computeFullExtentAndGridSizePerDimension(RA.ZWorldLimits,...
            RA.PixelExtentInWorldZ,translation(3));
        
        Rout = imref3d([numRowsOutput numColsOutput numPlanesOutput],XWorldLimitsOut,YWorldLimitsOut,ZWorldLimitsOut);        
    end
    
end

end

function [worldLimitsOut,numPixelsInDimOutput] = computeFullExtentAndGridSizePerDimension(inputWorldLimits,inputWorldPixelExtentInDim,translationInDim)

% The full bounding rectangle is the bounding rectangle that
% includes the original and translated images
worldLimitsTranslated = inputWorldLimits+translationInDim;
minInDim = min(inputWorldLimits(1),worldLimitsTranslated(1));
maxInDim = max(inputWorldLimits(2),worldLimitsTranslated(2));

worldLimitsFullIdeal = [minInDim maxInDim];
idealFullExtent = diff(worldLimitsFullIdeal);

% Compute the number of pixels necessary to capture the entire full
% bounding box at the input image resolution. If the full extent is
% not evenly divisible by the input image resolution, use ceil to
% guarantee that we completely capture the full bounding box at the
% input image resolution.
numPixelsInDimOutput = ceil(idealFullExtent ./ inputWorldPixelExtentInDim);

% Compute the extent in world units of the output image, determined
% by the input image resolution and the number of pixels in the output
% image along each dimension.
outputImageExtentInDim = numPixelsInDimOutput*inputWorldPixelExtentInDim;

% If the ideal full image extent is not evenly divisible by the
% input image resolution, then the ceil will have added additional
% image extent. Compute the additional image extent.
addedImageExtentInDim  = outputImageExtentInDim - idealFullExtent;

% Add the additional image extent in each dimension on one side of
% the output image. Increase the full bounding box on the side that
% the translation vector points toward. This allows for a gradual
% transition as the translated image moves in sub-pixel increments.
if translationInDim >=0
    worldLimitsOut = worldLimitsFullIdeal + [0 addedImageExtentInDim];
else
    worldLimitsOut = worldLimitsFullIdeal + [-addedImageExtentInDim 0];
end

end


function [A,translation,method,outputView,fillValues] = parseInputs(varargin)

cachedInterpolationMethod = 'linear';
cachedOutputView = 'same';

supportedNumericClasses = {'uint8','uint16','uint32','int8','int16',...
    'int32','single','double','logical'};

p = inputParser();
p.addRequired('A',@validateInputImage);
p.addRequired('TRANSLATION',@validateTranslation);
p.addOptional('METHOD',cachedInterpolationMethod,@validateInterp);
p.addParameter('OutputView',true,@validateOutputView)
p.addParameter('FillValues',0,@validateFillValues)

p.parse(varargin{:});

A               = p.Results.A;
translation     = double(p.Results.TRANSLATION);
method          = cachedInterpolationMethod;
outputView      = cachedOutputView;
fillValues      = p.Results.FillValues;

if any(strcmp(method,{'linear','cubic'}));
    method = strcat('bi',method);
end

postValidateTranslation(A,translation);
checkFillValues(fillValues,A,translation);

%----------------------------------
    function TF = validateInputImage(A)
        
        supportedImageAttributes = {'nonsparse','finite','nonempty'};
        
        validateattributes(A,supportedNumericClasses,supportedImageAttributes,mfilename,'A');
        
        TF = true;
        
    end

%-----------------------------------------
    function TF = validateInterp(interpMethod)
        
        cachedInterpolationMethod = validatestring(interpMethod,{'nearest','linear','cubic','bilinear','bicubic'},...
            mfilename,'METHOD');
        
        TF = true;
        
    end

%---------------------------------------------
    function TF = validateOutputView(boundingBox)
        
        cachedOutputView = validatestring(boundingBox,{'same','full'},...
            mfilename,'OutputView');
        
        TF = true;
        
    end

end

function TF = validateFillValues(fillVal)

validateattributes(fillVal,{'numeric'},...
    {'nonempty','nonsparse'},'imtranslate','FillValues');

TF = true;

end

%---------------------------------------------
function TF = validateTranslation(translation)

supportedNumericClasses = {'uint8','uint16','uint32','int8','int16',...
    'int32','single','double'};

validateattributes(translation,supportedNumericClasses,{'nonempty','vector','real','nonsparse','finite'},...
    mfilename,'TRANSLATION');

if ~any(numel(translation) == [2 3])
    error(message('images:imtranslate:invalidTranslationLength','TRANSLATION'));
end

TF = true;

end

%----------------------------------------------
function postValidateTranslation(A,translation)

if numel(translation) == 3
    if ndims(A) ~= 3
        error(message('images:imtranslate:threeDimensionalTranslationImageMismatch','A','TRANSLATION'));
    end
end

end

function [R_A,varargin] = preparseSpatialReferencingObjects(varargin)
% This should be abstracted into a package

if (nargin > 1) && isa(varargin{2},'imref2d')
    validateattributes(varargin{2},{'imref2d'},{'scalar','nonempty'},'imwarp','RA');
    R_A = varargin{2};
    varargin(2) = [];
else
    % We don't want to actually assign the default spatial referencing
    % object until the rest of the input arguments have been validated.
    % Assign empty spatial referencing arguments as a flag that we need to
    % assign the identity spatial referencing object after input
    % parsing/validation has finished.
    translation = varargin{2};
    validateTranslation(translation);
    postValidateTranslation(varargin{1},translation);
    if (numel(translation) == 2)
        R_A = imref2d.empty();
    else
        R_A = imref3d.empty();
    end
    
end

end

function checkSpatialRefAgreementWithInputImage(A,RA)
% This should be abstracted into a package

if ~sizesMatch(RA,A)
    error(message('images:imwarp:spatialRefDimsDisagreeWithInputImage','ImageSize','RA','A'));
end

end

function checkFillValues(fillValues,inputImage,translation)

planeAtATimeProblem = numel(translation)==2  && ~ismatrix(inputImage);

scalarFillValuesRequired = ~planeAtATimeProblem;
if scalarFillValuesRequired && ~isscalar(fillValues)
    error(message('images:imtranslate:scalarFillValueRequired','''FillValues'''));
end

if planeAtATimeProblem && ~isscalar(fillValues)
    sizeImage = size(inputImage);
    
    % MxNxP input image is treated as a special case. We allow [1xP] or
    % [Px1] fillValues vector in this case.
    validFillValues = isequal(sizeImage(3:end),size(fillValues)) ||...
        (isequal(ndims(inputImage),3) && isvector(fillValues)...
        && isequal(length(fillValues),sizeImage(3)));
    
    if ~validFillValues
        error(message('images:imwarp:fillValueDimMismatch','''FillValues''','''FillValues''','A'));
    end
end

end





