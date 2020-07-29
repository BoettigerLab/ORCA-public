function [normByDist,normMap] = ComputeNormMap(refMap,varargin)
% normByDist - power law fit
% normMap - normalized relative to median distance of each n-diag

defaults = cell(0,3);
defaults(end+1,:) = {'start','integer',4};
defaults(end+1,:) = {'stop','integer',8};
defaults(end+1,:) = {'method',{'powerlaw','median'},'powerlaw'};
defaults(end+1,:) = {'power','float',0};
pars = ParseVariableArguments(varargin,defaults,mfilename);


    nHybes= size(refMap,1);
    % compute distance normalization
    
    norm = ones(nHybes-1,1);
    normMap = zeros(size(refMap));
    for r=1:nHybes-1
        norm(r)  = nanmedian(diag(refMap,r));
        selMap = boolean(diag(true(nHybes-r,1),r) + diag(true(nHybes-r,1),-r));
        normMap =normMap+ norm(r)*double(selMap);
    end
    medMap = normMap; 
    
    x=(1:nHybes-1)';
    if pars.power == 0
        ftype = fittype('c*x+b','coeff',{'c','b'},'ind','x');
        xx = log10(x(pars.start:end-pars.stop));
        yy = log10(norm(pars.start:end-pars.stop));
        skip = isinf(yy) | yy==0 | isnan(yy);
        xx = xx(~skip); 
        yy = yy(~skip);
        linfit = fit(xx,yy,ftype,'StartPoint',[1,0]);
        y = x.^linfit.c*10^linfit.b;
        % figure(4); clf; loglog(norm); hold on; plot(x,y,'r.');
    else
        y = x.^pars.power;
    end

    % compute norm by distance
    normMap = zeros(size(refMap));
    normByDist = normMap;
    for r=1:nHybes-1
        selMap = boolean(diag(true(nHybes-r,1),r) + diag(true(nHybes-r,1),-r));
        normMap = normMap+ y(r)*double(selMap);
        normByDist = normByDist + norm(r)*double(normMap);
    end
if strcmp(pars.method,'median') && nargout == 1
    normByDist = medMap;
end
% figure(4); clf; imagesc(normMap); colorbar;