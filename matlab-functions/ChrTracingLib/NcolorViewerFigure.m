function figData = NcolorViewerFigure(varargin)
% mlapp does not allow calls to figure functions, imagesc, colormap 
% this is a work-around
    defaults = cell(0,3);
    defaults(end+1,:) = {'figHandle','freeType',100};
    defaults(end+1,:) = {'image','array',zeros(4)};
    defaults(end+1,:) = {'display',{'imagesc','colorbar','none'},'imagesc'};
    defaults(end+1,:) = {'names','cell',{}};
    defaults(end+1,:) = {'data','struct',struct()}; % function specific data can go here
    defaults(end+1,:) = {'clf','boolean',false};
    defaults(end+1,:) = {'callFig','boolean',true};
    pars = ParseVariableArguments(varargin,defaults,mfilename);
    
    % call figure
    if pars.callFig
        figData.handle = figure(pars.figHandle); 
    end
    % clear figure if requested
    if pars.clf
        clf;
        disp('figure cleared');
    end
    % execute requested type of display
    switch pars.display
        case 'imagesc'
            imagesc(pars.image); 
     
            nChns = pars.data.nChns;
            chnNames = pars.data.chnNames;
            colormap(pars.data.colormap);
            h = colorbar;
            set(h,'YTick',(1/nChns:1/(nChns):1)-.5/nChns,'YTickLabel', chnNames);            
        case 'none'
            % do nothing
        otherwise
            % do nothing
    end


    