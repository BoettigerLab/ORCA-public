function polysOut = CatPolys(cellPolys,varargin)
% Input is a cell array of polymers
% force concatinate polymers of different xy sizes along the z-dimension
% this version keeps the same upper left coordinate, tosses out data from
% columns >b and pads with NaNs if data is smaller than b.
% 
% Updated version - if b is not passed, b is autocomputed as the maximum
%   number of barcodes in any of the matrices in the cell array. 
%   Note, b is passed as a second variable, not a flag name pair, which
%   presevers backwards compatibility with earlier scripts and is compact.

% % handle two different formats of input:
% CatPolys(cellPolys,b,'Name',value)
% or
% CatPolys(cellPolys,'Name',value);

b=0;
if rem(length(varargin),2)==1
    b = varargin{1};
    varin = varargin(2:end);
else
    varin = varargin;
end

% now parse the inputs
defaults = cell(0,3);
defaults(end+1,:) = {'nRpts','integer',0};
pars = ParseVariableArguments(varin,defaults,mfilename);


if b==0
    arraySizes = cellfun(@size,cellPolys,'UniformOutput',false);
    isEmpty = cellfun(@length,arraySizes)~=3;
    arraySizes(isEmpty) = [];
    arraySizes = cat(1,arraySizes{:});
    b = max(arraySizes(1));
end

% main function
ds = size(cellPolys(:,1),1);
catPoly = cell(ds,1);
for d=1:ds
    currPolys = cellPolys{d};
    currPolys = currPolys(1:end-pars.nRpts,:,:);
    [h,dims,z] = size(currPolys);
    if dims==0
       continue
    end
    matBlank = nan(b,dims,z);
    h1 = min(h,b);
    matBlank(1:h1,:,:) = currPolys(1:h1,:,:);
    catPoly{d} = matBlank;
end
polysOut = cat(3,catPoly{:});