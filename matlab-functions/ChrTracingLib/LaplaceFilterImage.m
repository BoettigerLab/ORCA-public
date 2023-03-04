function imOut = LaplaceFilterImage(imIn,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'alpha','fraction',.2};
defaults(end+1,:) = {'hsize','integer',5};
defaults(end+1,:) = {'sigma','nonnegative',0};
defaults(end+1,:) = {'invert','boolean',true};
pars = ParseVariableArguments(varargin,defaults,mfilename);

[hs,ws,zs] = size(imIn);

imOut = imIn;

for z=1:zs
    temp = imIn(:,:,z);
    if pars.sigma == 0
    temp = imfilter(double(temp),fspecial('laplacian',pars.alpha),'symmetric');  
    else
        temp = imfilter(double(temp),fspecial('log',pars.hsize,pars.sigma),'symmetric');  
    end
    temp(temp>0) = 0; 
    temp = (max(temp(:)) - temp); % figure(12); clf; imagesc(temp); colorbar;
    temp = uint16(temp);
    imOut(:,:,z) = temp;
end