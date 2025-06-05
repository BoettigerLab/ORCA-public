function [fig,ax] = MovieSlider(imIn,spotTable,varargin)
% Overlay an input image and a generous threshold of find spots (e.g. from
% DaoFit) and allow user to adjust threshold and interactively see which
% spots are which. 

% global sliderValues;


    defaults = cell(0,3);
    defaults(end+1,:) = {'firstFrame','integer',1}; %  
    defaults(end+1,:) = {'color','colormap',[1,0,0]}; %  
    pars = ParseVariableArguments(varargin,defaults,mfilename);

    
%% build figure and slider
    fig = uifigure;
    g = uigridlayout(fig);
    g.RowHeight = {'1x','fit','fit','fit'};
    g.ColumnWidth = {'1x'};
    ax = uiaxes(g);
    [h,w,] = size(imIn); 
    ax.XLim = [1,w];
    ax.YLim = [1,h];

   

    s0 = min(spotTable.significance);
    s1 = max(spotTable.significance);
    sliderValues{1} = [s0,s1];

    h0 = min(spotTable.height);
    h1 = max(spotTable.height);
    sliderValues{2} = [h0,h1]; % current min / max

    w0 = min(spotTable.xsigma);
    w1 = max(spotTable.xsigma);
    sliderValues{3} = [w0,w1]; % current min / max

    fig.UserData.sliderValues = sliderValues;

    sld_sig = uislider(g, ...
    "range",...
    "Limits",[s0 s1], ...
    "Value",sliderValues{1});
    sld_sig.ValueChangingFcn = @(src,event) updateRange1(src,event,imIn,spotTable,sliderValues,ax,pars,fig);

    sld_height = uislider(g, ...
    "range",...
    "Limits",[h0,h1], ...
    "Value",sliderValues{2});
    sld_height.ValueChangingFcn = @(src,event) updateRange2(src,event,imIn,spotTable,sliderValues,ax,pars,fig);

    sld_spotWidth = uislider(g, ...
    "range",...
    "Limits", [w0,w1], ...
    "Value",sliderValues{3});
    sld_spotWidth.ValueChangingFcn = @(src,event) updateRange3(src,event,imIn,spotTable,sliderValues,ax,pars,fig);

    % show movie
    MovieUpdate(imIn,spotTable,sliderValues,ax,pars,fig);
end


function updateRange1(src,event,imIn,spotTable,sliderValues,ax,pars,fig) % ,imTile1,imTile1,T,ax
    sliderValues{1} = round(event.Value);
    fig.UserData.sliderValues{1}= round(event.Value);
    MovieUpdate(imIn,spotTable,sliderValues,ax,pars,fig);
end

function updateRange2(src,event,imIn,spotTable,sliderValues,ax,pars,fig) % ,imTile1,imTile1,T,ax
    sliderValues{2} = round(event.Value);
    fig.UserData.sliderValues{2}= round(event.Value);
    MovieUpdate(imIn,spotTable,sliderValues,ax,pars,fig);
end

function updateRange3(src,event,imIn,spotTable,sliderValues,ax,pars,fig) % ,imTile1,imTile1,T,ax
    sliderValues{3} = round(event.Value);
    fig.UserData.sliderValues{3}= round(event.Value);
    MovieUpdate(imIn,spotTable,sliderValues,ax,pars,fig);
end

function MovieUpdate(imIn,spotTable,sliderValues,ax,pars,fig)
    % display options 
    sliderValues = fig.UserData.sliderValues;
    keep = spotTable.significance > sliderValues{1} & spotTable.significance < sliderValues{1}(2);
    keep = keep & spotTable.height > sliderValues{2}(1)  & spotTable.height < sliderValues{2}(2);
    keep = keep & spotTable.xsigma > sliderValues{3}(1)  & spotTable.xsigma < sliderValues{3}(2);
    x = spotTable.x(keep);
    y = spotTable.y(keep);
    cla; imagesc(ax,imIn); hold(ax,'on');
    plot(ax,x,y,'o','color',pars.color);
end