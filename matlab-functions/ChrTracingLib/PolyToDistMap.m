function distMap = PolyToDistMap(poly)
% return the distance map, NxNxC given a polymer, Nx3xC
% a simple/lazy fucntion turning 5 lines of code into 1

[nB,~,nC] = size(poly);
distMap = nan(nB,nB,nC);
for c=1:nC
    distMap(:,:,c) = squareform(pdist(poly(:,:,c)));
end