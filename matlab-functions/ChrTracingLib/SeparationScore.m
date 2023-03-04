function [sepMedR2,sepMedR1] = SeparationScore(currMap,normMap,r1,r2,bH)

nB = size(currMap,1);

m0 = currMap./normMap;
m0(bH,:) = NaN;
m0(:,bH) = NaN;
m0( eye(nB)>0  ) = NaN;
m11 = m0(r1,r1);
m12 = m0(r1,r2);
m22 = m0(r2,r2);
sepMedR1 =  nanmedian(m12(:))./nanmedian(m11(:));
sepMedR2 = nanmedian(m12(:))./nanmedian(m22(:));