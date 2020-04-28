function [triMap,xRatio] = SquareToTriMap(squareMap,varargin)
% convert a square Hi-C map to a triangular Hi-C map

defaults = cell(0,3);
defaults(end+1,:) = {'top','fraction',1/4};
defaults(end+1,:) = {'zeroToInf','boolean',true};
defaults(end+1,:) = {'cleanDiag','boolean',true};

pars = ParseVariableArguments(varargin,defaults,mfilename);


squareMap = InterpMapNans(squareMap); % imrotate can't handle nans;
if pars.cleanDiag 
    diagValue = nanmean(diag(squareMap));
    rs = size(squareMap,1);
    squareMap(diag(true(rs,1))) = 1E-10;
end
triMap = imrotate(imresize(squareMap,4,'nearest'),45,'nearest');
w1 = size(squareMap,1);
[h,w2] = size(triMap);
triMap = triMap(round(h*pars.top:h/2+1),:); 
if pars.zeroToInf
    triMap(triMap==0) = inf;
end
if pars.cleanDiag
    triMap(triMap==1E-10) = diagValue;
end
xRatio = w1/w2;