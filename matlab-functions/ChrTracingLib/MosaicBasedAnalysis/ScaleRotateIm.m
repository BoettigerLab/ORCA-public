function im = ScaleRotateIm(im,alignValues)

scale1 = alignValues.rescale;
scale2 = alignValues.rescale2;
theta1 = alignValues.theta;
theta2 = alignValues.theta2;

 % scaling matrix
S1 = [scale1   0     0;
     0       scale1  0;
     0         0      1];
 % scaling matrix
S2 = [scale2   0     0;
     0       scale2  0;
     0         0      1];



if scale1 ~= 1
   im = imwarp(im,affine2d(S1),'OutputView',imref2d(size(im))  ); 
end
if theta1 ~= 0
    im = imrotate(im,theta1,'bilinear','crop');
end
if scale2 ~=1
    im = imwarp(im,affine2d(S2),'OutputView',imref2d(size(im))  ); 
end
if theta2 ~=0
    im = imrotate(im,theta2,'bilinear','crop');
end