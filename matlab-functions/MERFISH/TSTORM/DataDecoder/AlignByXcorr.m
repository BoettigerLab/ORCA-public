function [xshift,yshift] =  AlignByXcorr(Imbead,Im647,varargin)
% [xshift,yshift] =  AlignByXcorr(Imbead,Im647)
% 

disp('function called');

%--------------------------------------------------------------------------
%% Parse variable input
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
            otherwise
                error(['The parameter ''' parameterName ''' is not recognized by the function ''' mfilename '''.']);
        end
    end
end

%--------------------------------------------------------------------------
%% Main Function
%-------------------------------------------------------------------------
 

plotsOn = true;

% ----------------- Align images by cross-correlation
D = length(Imbead); 
xshift = zeros(D,1); % computed from beadfield
yshift = zeros(D,1); 
Ibead1 =  imadjust(Imbead{1},stretchlim(Imbead{1},0));
for d=1:D
    IbeadTemp = imadjust(Imbead{d},stretchlim(Imbead{d},0));
%     IbeadTemp = imadjust(Imbead{d},stretchlim(Imbead{d},[0.01 0.02]));
%     IbeadTemp = Imbead{d};
%     IbeadTemp = imadjust(Imbead{d},stretchlim(2000,0));
    [xshift(d),yshift(d)] = CorrAlign(Ibead1,IbeadTemp);
end

% ----------------- show plots of alignment
if plotsOn
    [H,W] = size(Imbead{1});
    x1 = W/2;
    y1 = H/2;   
    Ib = zeros(2*H,2*W,D,'uint16'); % bead image
    Im = zeros(2*H,2*W,D,'uint16'); % bead image
    I0 = zeros(2*H,2*W,D,'uint16');
    for d=1:D
       I0(y1:y1-1+H,x1:x1-1+W,d) = imadjust(Imbead{d},stretchlim(Imbead{d},0));
       Ib(y1+yshift(d):y1-1+H+yshift(d),x1+xshift(d):x1-1+W+xshift(d),d) =...
           imadjust(Imbead{d},stretchlim(Imbead{d},0));
       Im(y1+yshift(d):y1-1+H+yshift(d),x1+xshift(d):x1-1+W+xshift(d),d) =...
           imadjust(Im647{d},stretchlim(Im647{d},0));
    end

    figure(1); clf;
    subplot(1,2,1); Ncolor(I0); 
    xlim([W/2,3/2*W]);
    ylim([H/2,3/2*H]);
    title('before alignment','FontSize',15,'color','w');
    subplot(1,2,2); Ncolor(Ib); 
    xlim([W/2,3/2*W]);
    ylim([H/2,3/2*H]);
    title('after alignment','FontSize',15,'color','w');
    set(gcf,'color','k');

    figure(2); clf; Ncolor(Im);
    xlim([W/2,3/2*W]);
    ylim([H/2,3/2*H]);
    colormap(hsv(D));
    chandle = colorbar;
    set(get(chandle,'ylabel'),'String','Stain Round','FontSize',15);
    set(gcf,'color','k');
end