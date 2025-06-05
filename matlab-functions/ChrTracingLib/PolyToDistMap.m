function distMap = PolyToDistMap(poly,varargin)
% return the distance map, NxNxC given a polymer, Nx3xC
% a simple/lazy fucntion turning 5 lines of code into 1

defaults = cell(0,3);
defaults(end+1,:) = {'dims','integer',1:3};
pars = ParseVariableArguments(varargin,defaults,mfilename); 


[nB,~,nC] = size(poly);
distMap = nan(nB,nB,nC);
for c=1:nC
    distMap(:,:,c) = squareform(pdist(poly(:,pars.dims,c)));
end