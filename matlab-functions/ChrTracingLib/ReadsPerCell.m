function readsPerCell = ReadsPerCell(inputData)
% return a count of the number of reads per cell if passed a distance map
% nReads x nReads x nCells 
% or if passed a polymer map nReads x nDims x nCells

[n1,n2,n3] = size(inputData);
if n1==n2 % what was passed is a distmap
    readsPerCell = sqrt(sum(~isnan(reshape(permute(inputData,[3,2,1]),n3,n1^2)),2));
else % was a poly
    readsPerCell = squeeze(sum(~isnan(inputData(:,1,:))));
end
