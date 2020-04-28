function map = InsertGaps(map,gaps)
% note gaps must be added in order from low to high

for i=1:length(gaps)
    [nHybs_1,~,nCells] = size(map);
    g = gaps(i);
    map = [map(1:g-1,:,:); nan(1,nHybs_1,nCells); map(g:end,:,:)];
    map = [map(:,1:g-1,:), nan(nHybs_1+1,1,nCells), map(:,g:end,:)];
end