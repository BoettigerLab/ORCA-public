
function [pmap,hmap] = DistPvalue(dmap1,dmap2,varargin)
defaults = cell(0,3);
defaults(end+1,:) = {'stat',{'ranksum','kstest'},'ranksum'};
pars = ParseVariableArguments(varargin,defaults,mfilename);

N = size(dmap1,1);
pmap = zeros(N,N);
hmap = zeros(N,N);
for a=1:N
    for b=a+1:N
        v1 = squeeze(dmap1(a,b,:));
        v2 = squeeze(dmap2(a,b,:));
        v1 = v1(~isnan(v1));
        v2 = v2(~isnan(v2));
        if isempty(v1) || isempty(v2)
            pmap(a,b) = 1;
            hmap(a,b) = 0;
        else
            if strcmp(pars.stat,'ranksum')
                [pmap(a,b),hmap(a,b)] = ranksum(v1,v2);
            else
                [hmap(a,b),pmap(a,b)] = kstest2(v1,v2);
            end
        end
    end
end
pmap = pmap + pmap';
hmap = hmap + hmap';