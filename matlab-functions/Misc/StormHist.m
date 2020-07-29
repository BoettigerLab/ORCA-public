function imPlot = StormHist(mlist,varargin)
% Plot dots but stack their intensities when they overlap
% this is a sort of hybrid between plotting points and plotting a rendered
% image where the dotsize is on scale of localization precision.  
% There is an option to have this 'intensity scatterplot' overlaid ontop of
% a gray scale conventional image. 
%
% I first explored this plotting method in making figures for the
% OligoSecondaries paper. 

%% Default Input Parameters
imaxes = [];
convImage = [];
trim = false;
gain = .45;

%--------------------------------------------------------------------------
% Parse variable input
%--------------------------------------------------------------------------
if nargin > 1
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;
    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            case 'imaxes'
                imaxes = CheckParameter(parameterValue,'struct','imaxes');
            case 'convImage'
                convImage = CheckParameter(parameterValue,'array','convImage');
            case 'trim'
                trim = CheckParameter(parameterValue,'boolean','trim');
            case 'gain'
                gain = CheckParameter(parameterValue,'nonnegative','gain');
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

%% Main Function

if ~isempty(imaxes)
    im = list2img(mlist,imaxes,'dotsize',1,'N',1);
else
    [im,imaxes] = list2img(mlist,'dotsize',1,'N',1);
end
imS = STORMcell2img(im,'colormap',hot(256));

if ~isempty(convImage)
    [h,w] = size(imS);
    if trim
    conv = imresize(convImage(imaxes.ymin+1:imaxes.ymax+1,...
                              imaxes.xmin:imaxes.xmax),[h,w],'nearest');
    else
        conv = imresize(convImage,[h,w],'nearest');
    end
    imPlot = ind2rgb(imS,hot(2^15))*2^16 + gain*double(repmat(conv,1,1,3));
else
    imPlot = imS;
end

if nargout == 0
   Ncolor(imPlot); 
end
