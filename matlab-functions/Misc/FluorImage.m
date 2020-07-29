function FluorImage(varargin)
% Pretty Plotting function for fluorescent images.  Enforces that axis
% scaling is equal, sets background and colordef to black, stretches
% plotting axis of (or subplot axes) to fill the image).
% Adds a thin white box around the border (can be turned off).
%
% 
% % default parameters
% h = []; % figure handle
% trim = true; 
% box = true;
% 
% Alistair Boettiger 
% boettiger.alistair@gmail.com

% default parameters
h = []; % figure handle
trim = true; 
box = true;

if nargin > 0
    if (mod(length(varargin), 2) ~= 0 ),
        error(['Extra Parameters passed to the function ''' mfilename ''' must be passed in pairs.']);
    end
    parameterCount = length(varargin)/2;
    for parameterIndex = 1:parameterCount,
        parameterName = varargin{parameterIndex*2 - 1};
        parameterValue = varargin{parameterIndex*2};
        switch parameterName
            case 'figHandle'
                h = CheckParameter(parameterValue,'handle','figHandle');
            case 'trim'
                trim = CheckParameter(parameterValue,'boolean','trim');
            case 'box'
                box =  CheckParameter(parameterValue,'boolean','box');
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

if isempty(h);
    h = gcf;
end

colordef black;
set(gcf,'color','k');
subaxes = findobj(gcf,'type','axes');
for i=1:length(subaxes);
    try
        if box
         set(subaxes(i),'XTickLabel','','YTickLabel','','XTick',[],'YTick',[]);
        else
          axis(subaxes(i),'off'); 
        end
    catch er
        disp(er.message);
    end
    axis(subaxes(i),'image'); 
end
try
    if trim
        spaceplots(h,[.0 .0 .0 .0], [.0 .0]); 
    end
catch er
    disp(er.message);
end
    
