function [im,pads] = MergeImages(im1,im2,varargin)
% imOut = MergeImages(im1,im2)  - returns a imOut, a combination of image 1
%           and image 2 (im1 and im2),
% imOut = MergeImages(im1,im2,'autocontrast',false) sets autocontrast
%           option to false; 
% combine two images, even if of different size. Optionally allows im2 to
% be transltaed or rotated.  By default, images are auto-contrasted to
% equalize intensities without data loss (default range [0,1]). 

defaults = cell(0,3);
defaults(end+1,:) = {'xshift','float',0};
defaults(end+1,:) = {'yshift','float',0};
defaults(end+1,:) = {'padValue','float',0};
defaults(end+1,:) = {'rotateCW','float',0};
defaults(end+1,:) = {'contrast','boolean',true};
defaults(end+1,:) = {'high1','positive',1}; % contrast parameters
defaults(end+1,:) = {'high2','positive',1};
defaults(end+1,:) = {'low1','nonnegative',0};
defaults(end+1,:) = {'low2','nonnegative',0};
defaults(end+1,:) = {'shiftFirst','boolean',true};
defaults(end+1,:) = {'align',{'upperLeft','lowerRight','center','matchFirst'},'upperLeft'};

pars = ParseVariableArguments(varargin,defaults,mfilename); 
pars.xshift = round(pars.xshift); % required to be integer
pars.yshift = round(pars.yshift);  % required to be integer
if pars.contrast
    try
        im1 = IncreaseContrast(im1,'low',pars.low1,'high',pars.high1);
        im2 = IncreaseContrast(im2,'low',pars.low2,'high',pars.high2);
    catch er
        warning(er.message);
    end
end

im1_origin = [0,0]; % depreciated. will remove later;

im1o = im1;
im2o = im2;
[h1,w1,c1] = size(im1);
[h2,w2,c2] = size(im2);
%%
pads = [h2-h1, w2-w1];
if h2-h1 > 0
    hPad1 = true;
else
    hPad1 = false;
end
if w2-w1 > 0
    wPad1 = true;
else
    wPad1 = false;
end
if strcmp(pars.align,'upperLeft')      
    if hPad1
        im1 = padarray(im1,[h2-h1,0],pars.padValue,'post');
    else
        im2 = padarray(im2,[h1-h2,0],pars.padValue,'post');
    end
    if wPad1
        im1 = padarray(im1,[0,w2-w1],pars.padValue,'post');
    else
        im2 = padarray(im2,[0,w1-w2],pars.padValue,'post');
    end
elseif strcmp(pars.align,'lowerRight')
    if hPad1
        im1 = padarray(im1,[h2-h1,0],pars.padValue,'pre');
    else
        im2 = padarray(im2,[h1-h2,0],pars.padValue,'pre');
    end
    if wPad1
        im1 = padarray(im1,[0,w2-w1],pars.padValue,'pre');
    else
        im2 = padarray(im2,[0,w1-w2],pars.padValue,'pre');
    end
elseif strcmp(pars.align,'center')
    hpre = floor((h2-h1)/2);
    hpst = ceil((h2-h1)/2);
    wpre = floor((w2-w1)/2);
    wpst = ceil((w2-w1)/2);
    if hPad1
        im1 = padarray(im1,[hpre,0],pars.padValue,'pre');
        im1 = padarray(im1,[hpst,0],pars.padValue,'post');
    else
        im2 = padarray(im2,[-hpre,0],pars.padValue,'pre');
        im2 = padarray(im2,[-hpst,0],pars.padValue,'post');
    end
    if wPad1
        im1 = padarray(im1,[0,wpre],pars.padValue,'pre');
        im1 = padarray(im1,[0,wpst],pars.padValue,'post');
    else
        im2 = padarray(im2,[0,-wpre],pars.padValue,'pre');
        im2 = padarray(im2,[0,-wpst],pars.padValue,'post');
    end
elseif strcmp(pars.align,'matchFirst')
    if hPad1
       im2 = im2(1:h1,:,:);
    else
        im2 = padarray(im2,[h1-h2,0],pars.padValue,'post');
    end
    if wPad1
        im2 = im2(:,1:w1,:);
    else
        im2 = padarray(im2,[0,w1-w2],pars.padValue,'post');
    end
end
    
im2s = im2;
if pars.shiftFirst
    im2s = TranslateImage(im2s,pars.xshift,pars.yshift);
end
im2s = imrotate(im2s,pars.rotateCW,'bilinear','crop');
if ~pars.shiftFirst
    im2s = TranslateImage(im2s,pars.xshift,pars.yshift);
end
im = cat(3,im1,im2s);

if nargout == 0
    Ncolor(im);
end