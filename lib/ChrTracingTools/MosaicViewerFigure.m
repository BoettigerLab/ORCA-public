function figData = MosaicViewerFigure(varargin)
% mlapp does not allow calls to figure functions, imagesc, colormap 
% this is a work-around
    defaults = cell(0,3);
    defaults(end+1,:) = {'figHandle','freeType',100};
    defaults(end+1,:) = {'image','array',zeros(4)};
    defaults(end+1,:) = {'display',{'imagesc','Ncolor','NcolorGUI','none','drawBoxes','clearBoxes'},'imagesc'};
    defaults(end+1,:) = {'names','cell',{}};
    defaults(end+1,:) = {'data','struct',struct()}; % function specific data can go here
    defaults(end+1,:) = {'clf','boolean',false};
    defaults(end+1,:) = {'callFig','boolean',true};
    pars = ParseVariableArguments(varargin,defaults,mfilename);
    
    % call figure
    if pars.callFig
        try
            figData.handle = figure(pars.figHandle); 
        catch
           figure(100); 
        end
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
            colormap(gray);        
        case 'Ncolor'
            appHandle = NcolorApp(pars.image,'names',pars.names,...
            'figHandle',pars.figHandle,...
            'appHandle',pars.data.ncaHandle);            
            figData.appHandle = appHandle; % pass app handle back
        case 'NcolorGUI' % archived, to delete
            NcolorGUI(pars.image,'names',pars.names,...
            'figHandle',pars.figHandle,'instanceID',1,'restart',false);
        case 'drawBoxes'
            % unpack some additional plotting objects
            xy = pars.data.xy;
            f = pars.data.f;
            h = pars.data.imSize(1); 
            w = pars.data.imSize(2);
            hold on; 
            text(xy(f,1),xy(f,2),num2str(f),'color','c');
            r = rectangle('position',[xy(f,1),xy(f,2),w,h]);
            r.EdgeColor = [1 0 1];
        case 'clearBoxes'
            % handled in MosaicViewerApp
        case 'none'
            % do nothing
        otherwise
            % do nothing
    end


    