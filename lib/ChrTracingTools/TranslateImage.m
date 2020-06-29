function imageOut = TranslateImage(imageIn,xshift,yshift,varargin)
% imageOut = TranslateImage(imageIn,xshift,yshift)
%
% Note on function:
%   Can be called on 2D multi-color or 3D images.
%   For 3D images, only integer pixel values are allowed.  
%   For 2D images, decimal pixels are allowed. The parameter 'upsample'
%   determines the degree of 2D interpolation, which determines the
%   accuracy at which decimal pixels are interpreted.  Set 'upsample' to 2
%   to allow half pixel accuracy, 4 to allow quarter pixel, etc.
%
% Alistair Boettiger
% Updated 4/8/17 to use faster padarray 
% Updated 8/5/17 to clarify use of upsample on 3D images  

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'padValue', 'freeType', 0};
defaults(end+1,:) = {'verbose', 'boolean', false};
defaults(end+1,:) = {'upsample','positive',1};
defaults(end+1,:) = {'zshift','float',0};
defaults(end+1,:) = {'usePad','boolean',true}; 
defaults(end+1,:) = {'trim','boolean',true}; 
defaults(end+1,:) = {'resize3D','boolean',false}; % sometimes 3rd dim is color  
% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);


if parameters.upsample ~= 1
    xshift = round(xshift*parameters.upsample);
    yshift = round(yshift*parameters.upsample);
    if parameters.resize3D % 3rd dim is z-position
        zshift =  parameters.zshift*parameters.upsample; % round(parameters.zshift*parameters.upsample);
        imageIn = imresize3(imageIn,parameters.upsample); % note, imresize does not rescale zshift
    else % 3rd dim is color
        zshift =   parameters.zshift; % round(parameters.zshift*parameters.upsample);
        imageIn = imresize(imageIn,parameters.upsample); % note, imresize does not rescale zshift
    end
else
    zshift = parameters.zshift;
end

% enforce 
xshift = round(xshift);
yshift = round(yshift);
zshift = round(zshift); 


dType = class(imageIn);
[h,w,c] = size(imageIn);
imageOut = imageIn;
if parameters.usePad   % padarray is faster than concatination 
    if xshift >= 0 && yshift >= 0
        imageOut = padarray(imageOut,[yshift,xshift],parameters.padValue,'pre');
        if parameters.trim; imageOut = imageOut(1:h,1:w,:); end
    elseif xshift < 0 && yshift < 0
        imageOut = padarray(imageOut,[-yshift,-xshift],parameters.padValue,'post');
        if parameters.trim; imageOut = imageOut(end-h+1:end,end-w+1:end,:); end
    elseif xshift >= 0 && yshift < 0
        imageOut = padarray(imageOut,[0,xshift],parameters.padValue,'pre');
        imageOut = padarray(imageOut,[-yshift,0],parameters.padValue,'post');
        if parameters.trim; imageOut = imageOut(end-h+1:end,1:w,:); end
    elseif xshift < 0 && yshift >= 0 
        imageOut = padarray(imageOut,[0,-xshift],parameters.padValue,'post');
        imageOut = padarray(imageOut,[yshift,0],parameters.padValue,'pre');
        if parameters.trim; imageOut = imageOut(1:h,end-w+1:end,:); end
    end

    if zshift > 0
        imageOut = padarray(imageOut,[0,0,zshift],parameters.padValue,'pre');
        if parameters.trim; imageOut = imageOut(:,:,1:c); end
    elseif zshift < 0
        imageOut = padarray(imageOut,[0,0,-zshift],parameters.padValue,'post'); 
        if parameters.trim; imageOut = imageOut(:,:,end-c+1:end); end
    end

else

% % old version ( this is slower than padarray)
% 
% xshift 
if xshift > 0
    % imageOut =[parameters.padValue*ones(h,xshift), imageOut(:,1:end-xshift+1)];
    imageOut =[parameters.padValue*ones(h,xshift,c,dType), imageOut(:,1:end-xshift,:)];
elseif xshift < 0
    imageOut = [imageOut(:,-xshift+1:end,:), parameters.padValue*ones(h,-xshift,c,dType)];
end

if yshift > 0
    imageOut =[parameters.padValue*ones(yshift,w,c,dType); imageOut(1:end-yshift,:,:)];
elseif yshift < 0
    imageOut = [imageOut(-yshift+1:end,:,:); parameters.padValue*ones(-yshift,w,c,dType)];
end

if parameters.zshift > 0
    imageOut = cat(3,parameters.padValue*ones(h,w,zshift,dType),imageOut(:,:,1:end-zshift));
elseif parameters.zshift < 0 
    imageOut = cat(3,imageOut(:,:,-zshift+1:end),parameters.padValue*ones(h,w,-zshift,dType));
end

end

if parameters.verbose
    disp(['x-shifted ',num2str(xshift/parameters.upsample)]);
    disp(['y-shifted ',num2str(yshift/parameters.upsample)]);
    if zshift ~= 0
        disp(['z-shifted ',num2str(zshift/parameters.upsample)]);
    end
end

if parameters.upsample ~= 1
    imageOut = imresize(imageOut,1/parameters.upsample);
end

