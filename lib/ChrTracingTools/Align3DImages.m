function [xshift,yshift,zshift] = Align3DImages(Im1,Im2,varargin)
% [xshift,yshift,zshift] = AlignImages3D(Im1,Im2,varargin)
% defaults(end+1,:) = {'xRange', 'array', [-2,2]};
% defaults(end+1,:) = {'yRange', 'array', [-2,2]};
% defaults(end+1,:) = {'zRange', 'array', [-2,2];};
% defaults(end+1,:) = {'upsample','positive',1};
% defaults(end+1,:) = {'maxScan','integer',1E4};
% defaults(end+1,:) = {'verbose','boolean',true};

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'xRange', 'array', [-2,2]};
defaults(end+1,:) = {'yRange', 'array', [-2,2]};
defaults(end+1,:) = {'zRange', 'array', [-2,2]};
defaults(end+1,:) = {'upsample','positive',1};
defaults(end+1,:) = {'maxScan','integer',1E4};
defaults(end+1,:) = {'verbose','boolean',true};

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);


% -------------------------------------------------------------------------
% Main Function
% -------------------------------------------------------------------------
M1 = double(imresize(Im1,parameters.upsample));
M2 = double(imresize(Im2,parameters.upsample));

k = 0;
zs = parameters.zRange(1)*parameters.upsample:parameters.zRange(2)*parameters.upsample;
xs = parameters.xRange(1)*parameters.upsample:parameters.xRange(2)*parameters.upsample;
ys = parameters.yRange(1)*parameters.upsample:parameters.yRange(2)*parameters.upsample;
diffs = NaN(length(xs)*length(ys)*length(zs),1);

if length(diffs)>parameters.maxScan
    warning(['greater than ',num2str(parameters.maxScan),' combinations to scan, aborting.']);
    warning('Try increasing "maxScan" to continue');
    xshift = 0; 
    yshift = 0;
    zshift = 0;
    return
end
for z = zs
    if parameters.verbose
        disp(['3D image registration ',num2str(round(100*k/length(diffs))),'% complete']);
    end
    for x = xs
        for y= ys
            k = k+1;
            M2s = TranslateImage(M2,x,y,'zshift',z);
            Md = M1-M2s;
            temp = abs(Md);
            diffs(k) = nanmean(temp(:));
        end
    end
end
[~,i] = min(diffs);
[x,y,z] = ind2sub([length(xs),length(ys),length(zs)],i);
xshift = xs(x)/parameters.upsample;
yshift = ys(y)/parameters.upsample;
zshift = zs(z)/parameters.upsample;




