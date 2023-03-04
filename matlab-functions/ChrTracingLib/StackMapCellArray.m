function mapMat = StackMapCellArray(maps)

sz = cellfun(@size,maps,'UniformOutput',false);
sz = cat(1,sz{:});
nb = max(sz(:,1));
nd = max(sz(:,2));

mapCell = maps;
needPad = find(sz(:,1)' < nb);
for n=needPad
    currMap = maps{n};
    padDim1 = nb - sz(n,1);
    padDim2 = nd - sz(n,2);
    newMap = padarray(currMap,[padDim1,padDim2],'post');
    mapCell{n} = newMap;
end
mapMat = cat(3,mapCell{:});
