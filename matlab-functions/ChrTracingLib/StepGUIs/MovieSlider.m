function [fig,ax] = MovieSlider(movieIn,varargin)

    defaults = cell(0,3);
    defaults(end+1,:) = {'firstFrame','integer',1}; %  
    defaults(end+1,:) = {'options','struct',struct()}; %  
    pars = ParseVariableArguments(varargin,defaults,mfilename);

    t = pars.firstFrame;
    options = pars.options;


%  [fig,ax] = MovieSlider(movieIn,t,options)
%% Inputs
% movieIn - a 4D object H x W x chn x slider-value (time)
%  t - first frame to show
%  options 
    
%% build figure and slider
    fig = uifigure;
    g = uigridlayout(fig);
    g.RowHeight = {'1x','fit'};
    g.ColumnWidth = {'1x'};
    ax = uiaxes(g);
    [h,w,~,T] = size(movieIn); 
    ax.XLim = [1,w];
    ax.YLim = [1,h];

    sld = uislider(g, ...
    "Limits",[1 T], ...
    "Value",t);
    sld.ValueChangingFcn = @(src,event) updateRange(src,event,movieIn,ax,options);

    % "MinorTicks",[],...
    % "MajorTicksMode",'manual',...
    % "MinorTicksMode",'manual',...

    % show movie
    MovieUpdate(movieIn,t,ax,options);
end


function updateRange(src,event,movieIn,ax,options) % ,imTile1,imTile1,T,ax
    t = round(event.Value);
    MovieUpdate(movieIn,t,ax,options);
end

function MovieUpdate(movieIn,t,ax,options)
    % display options
    if  isfield(options,'contrast')
        movieOut = IncreaseContrast(movieIn(:,:,:,t),'high',options.contrast(1),'low',options.contrast(2));
        frameOut = Ncolor(movieOut);
        imagesc(ax,frameOut);
    else
        frameOut = Ncolor(movieIn(:,:,:,t));
        imagesc(ax,frameOut);
    end
        % add post-plotting options 
    if  isfield(options,'labels')
        labels = options.labels;
        text(ax,labels.x,labels.y,labels.text,'color',labels.color);
    end
end