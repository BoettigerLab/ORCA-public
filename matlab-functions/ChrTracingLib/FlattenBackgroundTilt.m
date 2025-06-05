function [imOut,bkd] = FlattenBackgroundTilt(imIn)


[h,w,z] = size(imIn);
imIn2D = max(imIn,[],3);

x=1:h; 
y1 = mean(imIn2D,2) ;
linFit_1 = polyfit(x,y1,1);
yy = polyval(linFit_1,x);
x=1:w; 
y2 = mean(imIn2D,1) ;
linFit_2 = polyfit(x,y2,1);
xx = polyval(linFit_2,x);
[XX,YY] = meshgrid(xx,yy);
% YY = flipud(YY); % images start 0 at top
bkd = double(XX.*YY); 

% figure(6); clf; imagesc(XX.*YY); colormap('parula'); colorbar;
% imN = makeuint(double(imIn2D)./double(XX.*YY),16);
% figure(6); clf; imagesc(imN); colormap('parula'); colorbar;
% figure(7); clf; subplot(2,1,1); bar(y1);%  title()
% subplot(2,1,2); bar(y2);%


% imOut = double(imIn)./ repmat(bkd,1,1,z);
imOut = makeuint(double(imIn)./ repmat(bkd,1,1,z),16);
