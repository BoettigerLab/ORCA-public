function [imOut,norm] = NormMap(imIn,varargin)
% imOut - normalized relative to median distance of each n-diag

defaults = cell(0,3);
defaults(end+1,:) = {'start','integer',4}; % not used
defaults(end+1,:) = {'stop','integer',8}; % not used
defaults(end+1,:) = {'method',{'square','diag'},'diag'};
defaults(end+1,:) = {'max','positive',inf};
pars = ParseVariableArguments(varargin,defaults,mfilename);

if strcmp(pars.method,'diag')
    nHybes= size(imIn,1);
    % compute distance normalization
    norm = ones(nHybes-1,1);
    normMap = zeros(size(imIn));
    for r=1:min(nHybes-1, pars.max)
        norm(r)  = nanmedian(diag(imIn,r));
        selMap = boolean(diag(true(nHybes-r,1),r) + diag(true(nHybes-r,1),-r));
        normMap =normMap+ norm(r)*double(selMap);
    end
    imOut = normMap; 

else
    % original behavior, no longer default
    [nRows,nCols] = size(imIn);
    im = imIn;
    norm = nanmean(im,1);
    im = im./repmat(norm,nRows,1);
    norm2 = nanmean(im,2);
    im = im./repmat(norm2,1,nCols);
    imOut = im;

end