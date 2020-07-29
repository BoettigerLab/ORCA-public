function ImageWithContrastSlider(im,varargin)

    defaults = cell(0,3);
    defaults(end+1,:) = {'figHandle','freeType',[]};
    defaults(end+1,:) = {'sliderLabel','string','contrast'};
    defaults(end+1,:) = {'embedSlider','boolean',true};
    defaults(end+1,:) = {'showSlider','boolean',true};
    pars = ParseVariableArguments(varargin,defaults,mfilename);

    if isempty(pars.figHandle)
        fh = figure();
    else
        fh = figure(pars.figHandle);
    end

    clf; imagesc(im); colorbar;
    if pars.showSlider
    % cannot embed a slider in a normal figure (only "uifigure")
    %    uifigures lack most of the functionality I like from real figures...
    fs = uifigure('Position',[500,500,400,50]); 
    pnl = uipanel(fs,'Position',[10,10,400,40]);
    sld = uislider(pnl,'Position',[60 30 300 3]);
    txt = uilabel(pnl,'Position',[10,10,30,30]);
    txt.Text = string(pars.sliderLabel); 
    cMax = max(im(:));
    cMin = min(im(:));
    sld.Limits = [cMin,cMax];
    sld.Value = cMax;
    sld.UserData.sH = [];
    sld.ValueChangedFcn = @(sld,event) UpdateSlider(sld,im);
    % sld.DeleteFcn = @(sld,event) DeleteSlider(sld,fs);  % disp('close when done'); 
    fh.DeleteFcn = @(sld,event) DeleteSlider(sld); % doesn't seem to work
    end

    function UpdateSlider(sld,im)
        clf; imagesc(im); colorbar;
        caxis([sld.Limits(1),sld.Value]);
    end

    function DeleteSlider(sld)
        delete(sld);
    end
end
