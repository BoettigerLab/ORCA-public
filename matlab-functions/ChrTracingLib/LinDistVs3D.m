function ps = LinDistVs3D(imIn,varargin)
% imOut - normalized relative to median distance of each n-diag

defaults = cell(0,3);
defaults(end+1,:) = {'method',{'square','diag','powerlaw'},'diag'};
defaults(end+1,:) = {'max','positive',inf};
defaults(end+1,:) = {'power','double',-1}; % only used in powerlaw
defaults(end+1,:) = {'stat',{'median','mean'},'median'}; % only used in 'diag'

pars = ParseVariableArguments(varargin,defaults,mfilename);

nHybes = size(imIn,1);
ps = nan(nHybes,1); 
 for r=1:min(nHybes-1, pars.max)
    if strcmp(pars.stat,'median')
        ps(r)  = nanmedian(diag(imIn,r));
    elseif strcmp(pars.stat,'mean')
        ps(r)  = nanmean(diag(imIn,r));
    end
 end