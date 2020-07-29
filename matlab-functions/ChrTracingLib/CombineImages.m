function [imOut,blend] = CombineImages(im1,im2,varargin)
% imOut = CombineImages(im1,im2)
% mean takes the average of the two images, ignoring nans;
% sum takes the sum of the two images, ignoring nans.
% edgeBlur will remove the appearance of seems upon joining two images by
% detecting the zero/nan valued part of each image and use alpha feathering
% of each image in approaching this edge. This is useful in combining
% images that don't fully overlap.

defaults = cell(0,3);
defaults(end+1,:) = {'method',{'mean','sum','edgeBlur','last','first'},'mean'};
defaults(end+1,:) = {'blendFast','boolean',true};
defaults(end+1,:) = {'checkNaNs','boolean',false}; % should be made true for max backwards compatability 
defaults(end+1,:) = {'blend','freeType',[]};
pars = ParseVariableArguments(varargin,defaults,mfilename);

if pars.checkNaNs
    im1(isnan(im1)) = 0;
    im2(isnan(im2)) = 0;
end


switch(pars.method)
    case 'mean'
        imOut = nanmean(cat(3,im1,im2),3);
    case 'sum'
        imOut = nansum(cat(3,im1,im2),3);
    case 'last'
        imOut = im2;
    case 'first'
        imOut = im1;
    case 'edgeBlur'
        origClass = class(im1);
        im1 = double(im1);
        im2 = double(im2);
        %--------------new version---------------
        overlap = im1>0 & im2>0;
        % for speed, change the smallest number of elements to/from 0
        [h_i,w_i] =size(im1);
        hasSmallOverlap = sum(overlap(:)) < h_i*w_i/2; % small overlap
            if hasSmallOverlap
                im1_noOverlap = im1;
                im1_noOverlap(overlap)=0;
                im2_noOverlap = im2;
                im2_noOverlap(overlap)=0;
            end
        if isempty(pars.blend)
            if pars.blendFast
                imE = im1>0;    
                imE = imresize(imE,1/5);
                blend1 = GetBlendMask(imE,'symmetric',true);
                blend1 = imresize(blend1,[h_i,w_i]); % 'nearest'/'bilinear'  doesn't speed up much 
            else
                blend1 = GetBlendMask(im1,'symmetric',true);
            end
            blend2= 1-blend1;   
        else
           blend1 = pars.blend{1}; 
           blend2 = pars.blend{2};
        end
        
        imOut = zeros(h_i,w_i);
        imOut(overlap) = im1(overlap).*blend1(overlap)+im2(overlap).*blend2(overlap);
        if hasSmallOverlap
            imOut = imOut + im1_noOverlap + im2_noOverlap;
        else
            imOut(~overlap) = im1(~overlap)+im2(~overlap);
        end        
        % this is much faster than cat, nansum(data,3);
        imOut = cast(imOut,origClass);
        %  % --------- old version -----------
        % blend1 = GetBlendMask(im1); 
        % blend2= 1-blend1;           
        % imOut = nansum( cat(3, im1.*blend1, im2.*blend2 ), 3);
        % imOut = cast(imOut,origClass);
end
if nargout > 1
    blend = {blend1,blend2};
end
% 
%  figure(8); clf; 
%  Ncolor(IncreaseContrast(cat(3,im1,im2),'high',.9995));
%  figure(9); clf; 
%  Ncolor(IncreaseContrast(cat(3,im1_overlap,im2_overlap),'high',.9995));
%  figure(10); clf; 
%  subplot(1,3,1); imagesc(im1_overlap.*blend1);
%  subplot(1,3,2); imagesc(im2_overlap.*blend2);
%  subplot(1,3,3); Ncolor(IncreaseContrast(cat(3,im1_overlap.*blend1,im2_overlap.*blend2),'high',.9995));
%  colorbar;