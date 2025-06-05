function readsPerBar = ReadsPerBarcode(inputData)
% Returns the Read efficiency per barcode

[n1,n2,n3] = size(inputData);
if n1==n2 % what was passed is a distmap
    error('please pass a polymer')
%     % reshape into cells x nB^2
%      readsPerBar  = sqrt(sum(~isnan(reshape(permute(inputData,[3,2,1]),n3,n1^2)),2));
% readsPerBar  = sqrt(sum(~isnan(reshape(permute(inputData,[3,2,1]),n3,n1^2)),2));
% 
%     readsPerBar = sum(sum(~isnan(inputData),3)/n3,2)/n2;
else % was a poly
    readsPerBar = squeeze(sum(~isnan(inputData(:,1,:)),3))/n3;
    % squeeze(sum(~isnan(inputData(:,1,:))));
end
