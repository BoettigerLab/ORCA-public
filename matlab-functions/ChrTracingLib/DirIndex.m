function [dirI] = DirIndex(map,varargin)
% returns difference in upstream vs. downstream contacts/distances
% See InsulationScore

defaults = cell(0,3);
defaults(end+1,:) = {'w','positive',10};
defaults(end+1,:) = {'showPlot','boolean',false};
defaults(end+1,:) = {'threshold','float',100};

pars = ParseVariableArguments(varargin,defaults,mfilename);

nReads = size(map,1);
w = pars.w; % window size in bins
dirI = zeros(1,nReads);
for r=1:nReads % r = 20
    st = max(1,r-w);
    en = min(nReads,r+w);
    downstream = map(r,st:r-1);
    upstream = map(r,r+1:en);

    B = nanmedian(downstream);
    A = nanmedian(upstream);
    E= (A+B)/2;
    dirI(r) = (B-A)./abs(B-A)* (  (A-E)^2/E + (B-E)^2/E);
end

