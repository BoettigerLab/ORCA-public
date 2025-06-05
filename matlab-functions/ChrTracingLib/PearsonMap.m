function [pmap,rmap] = PearsonMap(mmap,varargin)
% input - contact map or average distance map
% 2nd input optional, a normalization map; 

defaults = cell(0,3);
defaults(end+1,:) = {'method',{'square','diag','powerlaw'},'diag'};
defaults(end+1,:) = {'max','positive',inf};
defaults(end+1,:) = {'power','float',-1}; % only used in powerlaw
defaults(end+1,:) = {'stat',{'median','mean'},'median'}; % only used in 'diag'


if rem(length(varargin),2)==0
    varin = varargin;
    rmap =[];
else
    varin = varargin(2:end);
    rmap = mmap./varargin{1};
end
pars = ParseVariableArguments(varin,defaults,mfilename);

if isempty(rmap)
    normM = NormMap(mmap,'parameters',pars);
    rmap = mmap./normM;
end

nB = size(rmap,1);
pmap = nan(nB,nB);
for b=1:nB
    for c=1:nB
        try
        goodData = ~isinf(rmap(b,:)) & ~isnan(rmap(b,:)) & ~isinf(rmap(:,c)') & ~isnan(rmap(:,c)');
        pmap(b,c) = corr(rmap(b,goodData)',rmap(goodData,c));
        catch
        end
    end
end
rmap = log2(rmap);

% figure(1); clf; imagesc(log2(rmap)); caxis([-3 3]); GetColorMap('BlueWhiteRed');
% figure(2); clf; imagesc(pmap); caxis([-1 1]); GetColorMap('BlueWhiteRed');
