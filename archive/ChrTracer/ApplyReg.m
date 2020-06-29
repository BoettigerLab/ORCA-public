function imOut = ApplyReg(imIn,regData)
% imOut = ApplyReg(imIn,regData)
% 

imOut = ImageTranslate(imIn,[regData.xshift,regData.yshift]);
imOut = imrotate(imOut, regData.angle,'bilinear','crop');
if ~isempty(regData.tform)
    imOut = imwarp(imOut,regData.tform.tobj,'OutputView',imref2d(size(imOut)));
end
