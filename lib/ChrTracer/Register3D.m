function [outputImage,shifts] = Register3D(refImage,inputImage,varargin)
% [outputImage,shifts] = Register3D(refImage,inputImage)
% 
% i3 = Register3D(i1,i2,'upsample',4,'showplots',false)
% 
% Alistair Boettiger
% Aug 4, 2017
% CC BY NC

defaults = cell(0,3);
defaults(end+1,:) = {'upsample','positive',4};
defaults(end+1,:) = {'showplots','boolean',false};
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'maxShiftXY','positive',inf};
defaults(end+1,:) = {'maxShiftZ','positive',inf};
defaults(end+1,:) = {'center','nonnegative',0};
defaults(end+1,:) = {'threshold','nonnegative',.6};
pars = ParseVariableArguments(varargin,defaults,mfilename);


[nRows,nCols,nStks] = size(refImage);
if pars.center==0
    pars.center = [nRows/2,nCols/2,nStks/2]+.5;
end

xi = max(round(pars.center(1)-pars.maxShiftXY),1);
xe = min(round(pars.center(1)+pars.maxShiftXY),nCols);
yi = max(round(pars.center(2)-pars.maxShiftXY),1);
ye = min(round(pars.center(2)+pars.maxShiftXY),nRows);
zi = max(round(pars.center(3)-pars.maxShiftZ), 1);
ze = min(round(pars.center(3)+pars.maxShiftZ), nStks);
    
i1 = refImage(yi:ye,xi:xe,zi:ze);
i2 = inputImage(yi:ye,xi:xe,zi:ze);

im1= imresize(i1,pars.upsample); % try imresize3, change all 3 dimensions 

bx =  [i1(1,:,:), i1(end,:,:)];
by =  [i1(:,1,:), i1(:,end,:)];
bz =  [i1(:,:,1), i1(:,:,end)];
edge1 = quantile( [bx(:) ; by(:); bz(:)],.9); % .8

bx =  [i2(1,:,:), i2(end,:,:)];
by =  [i2(:,1,:), i2(:,end,:)];
bz =  [i2(:,:,1), i2(:,:,end)];
edge2 = quantile( [bx(:) ; by(:); bz(:)],.9); % .8

% Get a small volume as template
im2= imresize(i2,pars.upsample);  % try imresize3, change all 3 dimensions 
im2 = ImageTranslate(im2,[0,0]);

% KEY: send edge values to zero, because we will pad with zeros when
% doing the correlation based alignment.  
if pars.threshold ~= 0
    im1 = im1 - 1*edge1;
    im2 = im2 - 1*edge2;
else
    theta = quantile(im1(:), pars.threshold);
    im1(im1 < theta) = 0;
    theta = quantile(im2(:), pars.threshold);
    im2(im2<theta) = 0;
end

% Calculate SDD between template and image
I_SSD=template_matching(im2,im1);
% Find maximum correspondence
[~,indmax] = max(I_SSD(:));
[y,x,z]=ind2sub(size(I_SSD),indmax);
% [y,x,z]=ind2sub(size(I_SSD),find(I_SSD==max(I_SSD(:))));
[h,w,d] = size(im1);
xshift = x-w/2;
yshift = y-h/2;
zshift = z-d/2;
shifts.xshift = xshift/pars.upsample;
shifts.yshift = yshift/pars.upsample;
shifts.zshift = zshift; % /pars.upsample; % updated with imresize3

if pars.verbose
    disp(['xshift ', num2str(shifts.xshift),...
        ' yshift ',num2str(shifts.yshift),...
        ' zshift ',num2str(shifts.zshift)]);
end
temp = TranslateImage(imresize(inputImage,pars.upsample),xshift,yshift ,'zshift',zshift);
i3 = imresize(temp,1/pars.upsample); 
i3(i3==0) = edge2;
outputImage = i3;
   
if pars.showplots
    figure(10); clf;
    subplot(2,2,1); Ncolor(IncreaseContrast( cat(3, max(refImage,[],3), max(inputImage,[],3)) ) );
    subplot(2,2,2); Ncolor(IncreaseContrast( cat(3, squeeze(max(refImage,[],2)), squeeze(max(inputImage,[],2))) ));

    subplot(2,2,3); Ncolor(IncreaseContrast( cat(3, max(refImage,[],3), max(outputImage,[],3)) ));
    subplot(2,2,4); Ncolor(IncreaseContrast( cat(3, squeeze(max(refImage,[],2)), squeeze(max(outputImage,[],2))) ));

    figure(11); clf;
    subplot(2,2,1); Ncolor(IncreaseContrast( cat(3, max(im1,[],3), max(im2,[],3)) ) );
    subplot(2,2,2); Ncolor(IncreaseContrast( cat(3, squeeze(max(im1,[],2)), squeeze(max(im2,[],2))) ));

    im3 = TranslateImage(im2,xshift,yshift ,'zshift',zshift);
    subplot(2,2,3); Ncolor(IncreaseContrast( cat(3, max(im1,[],3), max(im3,[],3)) ));
    subplot(2,2,4); Ncolor(IncreaseContrast( cat(3, squeeze(max(im1,[],2)), squeeze(max(im3,[],2))) ));
end


% 