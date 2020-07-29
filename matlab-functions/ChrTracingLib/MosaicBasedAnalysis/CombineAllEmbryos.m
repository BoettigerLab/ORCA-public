function [eMapsAbd,eTablesAbd] = CombineAllEmbryos(embDataAbd,varargin)
% [eMapsAbd,eTablesAbd] = CombineAllEmbryos(embDataAbd)

defaults = cell(0,3);
defaults(end+1,:) = {'rs','integer',[]};
defaults(end+1,:) = {'minMedDist','positive',150};
pars = ParseVariableArguments(varargin,defaults,mfilename);


numEmbs = length(embDataAbd);
eTablesAbd = cell(numEmbs,1);
eMapsAbd = cell(numEmbs,1);
for e=1:numEmbs
    tableTemp = embDataAbd(e).table;
    tableTemp.emb = e*ones(height(tableTemp),1);
    eTablesAbd{e} = tableTemp;
    eMapsAbd{e} = embDataAbd(e).maps;
end
eTablesAbd = cat(1,eTablesAbd{:});
eMapsAbd = cat(3,eMapsAbd{:});



if isempty(pars.rs)
    rs = 1:size(eMapsAbd,1);
else
    rs = pars.rs;
end

% Get rid of fixed point artefacts, trim to rs
map2 = eMapsAbd(rs,rs,:);
mapStack = reshape(map2,size(map2,1)*size(map2,2),size(map2,3));
% figure(3); clf; hist(nanmedian(mapStack),100);
c1 = nanmedian(mapStack) < pars.minMedDist;
eMapsAbd(:,:,c1) = nan(size(eMapsAbd,1),size(eMapsAbd,1),sum(c1));