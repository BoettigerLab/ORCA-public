function [xyz,f] = MergeXyXzLists(xy_list,xz_list,s,numHybes,varargin)
% defaults = cell(0,3);
% defaults(end+1,:) = {'zStep', 'integer', 100}; % nm  per zstep
% defaults(end+1,:) = {'npp', 'positive', 150}; % nm per pixel
%
%-------------------------------------------------------------------------%
% Alistair Boettiger (boettiger@stanford.edu)
% February 10, 2017
% Copyright: CC BY NC
%-------------------------------------------------------------------------%

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'zStep', 'integer', 100}; % nm  per zstep
defaults(end+1,:) = {'npp', 'positive', 150}; % nm per pixel

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);



spotS = xy_list.frame >= (s-1)*(numHybes+1)+1 & xy_list.frame < s*(numHybes+1);
x = double(xy_list.x(spotS));
y = double(xy_list.y(spotS));
f = xy_list.frame(spotS)-(s-1)*(numHybes+1);

spotS = xz_list.frame >= (s-1)*(numHybes+1)+1 & xz_list.frame < s*(numHybes+1);
xz = double(xz_list.x(spotS));
z = double(xz_list.y(spotS))*parameters.zStep/parameters.npp;
fz = xz_list.frame(spotS)-(s-1)*(numHybes+1);

% match xy data with xz data based on x-positions
xyz = zeros(length(x),3); 
k = 0;
for h=1:numHybes
    xs = x(f==h);
    ys = y(f==h);
    xzs = xz(fz==h);
    zs = z(fz==h);
    if isempty(xzs)
        xzs=NaN;
        zs =NaN;
    end
    for i=1:length(xs)
        k=k+1;
        [~,idx] = min((xs(i)-xzs).^2);
        xyz(k,:) = [xs(i),ys(i),zs(idx)];
    end
end

%-------------------------------------------------------------------------%

% plotting to match spots for z-fit and x-fit
troubleshoot =false;
if troubleshoot
    h=1;
    im_xy = permute(squeeze(max(cy3Spots{s}(:,:,ss:end-es,:),[],3)),[2,1,3]);  % 
    im_xz = permute(squeeze(max(cy3Spots{s}(:,:,ss:end-es,:),[],2)),[2,1,3]);  % 
    figure(1); clf; 
    subplot(1,2,1); imagesc(im_xy(:,:,h)); title('xy');
    hold on; plot(x(f==h),y(f==h),'ro');
    subplot(1,2,2); imagesc(im_xz(:,:,h)); title('xz');
    hold on; plot(xz(fz==h),z(fz==h),'ro');
end