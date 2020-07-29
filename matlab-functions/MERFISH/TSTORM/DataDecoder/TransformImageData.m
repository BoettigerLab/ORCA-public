function [imageData, parameters] = TransformImageData(imageData,fiducialData,varargin)
%--------------------------------------------------------------------------
% Necessary Inputs
% imageData/A structure array with elements equal to the number of
%   images to align. This structure must contain an mList field.
% fiducialData/A structure array with elements equal to the number of
%   images to align. This structure must contain an mList field. 
%--------------------------------------------------------------------------
% Outputs
% imageData/A structure array with elements equal to the number of
%   images to align. This structure must contain an mList field.
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%
%--------------------------------------------------------------------------
% Alistair Boettiger
% boettiger@fas.harvard.edu
% Jeffrey Moffitt (revisions and reorganization)
% jeffmoffitt@gmail.com
% September 9, 2014
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Creative Commons License CC BY NC SA
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'verbose', 'boolean', true};  % 
defaults(end+1,:) = {'printedUpdates', 'boolean', true}; 
    
% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2
    error('matlabSTORM:invalidArguments', 'requires imageData and fiducialData');
end


% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

if parameters.printedUpdates && parameters.verbose
    display('--------------------------------------------------------------');
    display('Shifting images');
end
fieldsToTransferFrom = {'tform', 'warpErrors', 'hasFiducialError', 'fiducialErrorMessage', 'uID'};
fieldsToTransferTo = {'tform', 'warpErrors', 'hasFiducialError', 'fiducialErrorMessage', 'fidUID'};
for j=1:length(imageData)
    switch class(fiducialData(j).tform)
        case 'affine2d'
            [imageData(j).mList.xc, imageData(j).mList.yc] = ...
                transformPointsForward(fiducialData(j).tform, double(imageData(j).mList.x), ...
                double(imageData(j).mList.y));
        otherwise
            [imageData(j).mList.xc, imageData(j).mList.yc] = tforminv(...
                fiducialData(j).tform, ...
                double(imageData(j).mList.x), ...
                double(imageData(j).mList.y));
    end
    for k=1:length(fieldsToTransferFrom)
        imageData(j).(fieldsToTransferTo{k}) = fiducialData(j).(fieldsToTransferFrom{k});
    end        
end
