function mapsOut = CatMaps(cellMaps,varargin)
% force concatinate maps of different xy sizes along the z-dimension
% this version keeps the same upper left coordinate, tosses out data from
% columns >b and pads with NaNs if data is smaller than b.
% 
% Updated version - if b is not passed, b is autocomputed as the maximum
%   number of barcodes in any of the matrices in the cell array. 
%   Note, b is passed as a second variable, not a flag name pair, which
%   presevers backwards compatibility with earlier scripts and is compact.
%
% see also CatPolys

%% parse inputs
% % handle two different formats of input:
% CatMaps(cellMaps,b,'Name',value)
% or
% CatMaps(cellMaps,'Name',value);

b=0;
if rem(length(varargin),2)==1 % if there is an odd number of entries in varargin 
    b = varargin{1};
    varin = varargin(2:end);
else
    varin = varargin;
end

% now parse the optional flag variables
defaults = cell(0,3);
defaults(end+1,:) = {'nRpts','integer',0};
pars = ParseVariableArguments(varin,defaults,mfilename);


if b==0
    arraySizes = cellfun(@size,cellMaps,'UniformOutput',false);
    isEmpty = cellfun(@length,arraySizes)~=3;
    arraySizes(isEmpty) = [];
    arraySizes = cat(1,arraySizes{:});
    b1 = max(arraySizes(:,1));
    b2 = max(arraySizes(:,2));
    b = max([b1,b2]);
end

%% Main function
ds = size(cellMaps(:,1),1);
catMap = cell(ds,1);
for d=1:ds
    currMaps = cellMaps{d};
    currMaps = currMaps(1:end-pars.nRpts,1:end-pars.nRpts,:);
    [h,w,z] = size(currMaps);
    mapBlank = nan(b,b,z);
    h1 = min(h,b);
    w1 = min(w,b);
    mapBlank(1:h1,1:w1,:) = currMaps(1:h1,1:w1,:);
    catMap{d} = mapBlank;
end
mapsOut = cat(3,catMap{:});