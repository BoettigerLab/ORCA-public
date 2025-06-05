function [d22,d12,d11] = NormMedianDistance(currMap,normMap,r1,r2,bH)

nB = size(currMap,1);

if ~isempty(normMap)
    m0 = currMap./normMap;
else
    m0 = currMap;
end
if ~isempty(bH)
    m0(bH,:) = NaN;
    m0(:,bH) = NaN;
end
m0( eye(nB)>0  ) = NaN;
m11 = m0(r1,r1);
m12 = m0(r1,r2);
m22 = m0(r2,r2);
d12 =  nanmedian(m12(:));
d11 = nanmedian(m11(:));
d22 = nanmedian(m22(:));