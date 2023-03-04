function imOut = IncreaseContrast(imIn,varargin)
% default stretch contast 0,1 
% equivelant to:  imadjust(imIn,stretchlim(imIn,[0,1]))
%
% updated to allow z-stack

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'low', 'nonnegative', 0};
defaults(end+1,:) = {'high', 'nonnegative', 1};
defaults(end+1,:) = {'zStack', 'boolean', false};
defaults(end+1,:) = {'showWarning', 'boolean', true};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'a 2D image matrix is required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
pars = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
%% Main Function
% -------------------------------------------------------------------------
[~,~,nChns] = size(imIn);
if ~(isa(imIn,'uint16') || isa(imIn,'uint8') ) 
    if pars.showWarning
        warning('image was not uint. May give unexpected results');
    end
    imIn = makeuint(imIn,16);
end

if ~pars.zStack
    if length(pars.high) < nChns
        high = repmat(pars.high,nChns,1);
    else
        high = pars.high;
    end
    if length(pars.low) < nChns
        low = repmat(pars.low,nChns,1);  
    else
        low = pars.low;
    end
    imOut = imIn;
    if low < high
        for i=1:nChns
            if high(i) ~= 0 
                imOut(:,:,i) = imadjust(imIn(:,:,i),stretchlim(imIn(:,:,i),[low(i),high(i)])); %   
            else
                data = imIn(:,:,i);
                minB = quantile(data(:),low(i));
                data( data<minB) = 0;
                imOut = data;
            end
        end
    else
        warning('min vlaue not less than max value');
    end

else
     imOut = double(imIn);
     bkd = quantile(imOut(:),pars.low);
     imOut = imOut - bkd;     
     satValue = quantile(imOut(:),pars.high);
     if isa(imIn,'uint16')
        imOut = uint16(2^16*imOut ./ satValue);
     elseif isa(imIn,'uint8')
         imOut = uint8(2^16*imOut ./ satValue);
     else
         imOut = imOut ./ satValue;
     end   
end

