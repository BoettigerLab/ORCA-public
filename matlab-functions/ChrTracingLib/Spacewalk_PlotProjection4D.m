function figData = Spacewalk_PlotProjection4D(image4D,varargin)
% popout resizable, full functioned figure
    defaults = cell(0,3);
    defaults(end+1,:) = {'figHandle','freeType',100};
% PlotProjection4D defaults
        % defaults for projection 
        % defaults(end+1,:) = {'projection',{'xy','xz'},'xy'};
        defaults(end+1,:) = {'fixBlack','boolean',true};
        % defaults(end+1,:) = {'mode',{'subplots','single'},'single'};  % render as subplots (slow and flexible) or as a single image  
        % defaults for image tile (shared)
        defaults(end+1,:) = {'autoContrast', 'boolean', true}; 
        defaults(end+1,:) = {'tileLabels', 'cell', {}}; 
        defaults(end+1,:) = {'numRows', 'positive', 4}; 
        defaults(end+1,:) = {'colormap','colormap','hsv'};
        % Defaults for TileImageStack ('mode','subplots');
        defaults(end+1,:) = {'gap', 'nonnegative', 1}; 
        defaults(end+1,:) = {'colorTiles', 'boolean', true}; 
        defaults(end+1,:) = {'showPanelNumber', 'boolean', true}; 
        defaults(end+1,:) = {'fontSize', 'positive', 10}; 
        defaults(end+1,:) = {'imWhite','freeType',[]};
        % Defaults for TileImage  ('mode' single)
        defaults(end+1,:) = {'multicolor', 'boolean', true}; 
        defaults(end+1,:) = {'showImage', 'boolean', false}; 
        defaults(end+1,:) = {'showLabels', 'boolean', true}; 
        defaults(end+1,:) = {'zoomBox', 'array', []}; 
        defaults(end+1,:) = {'labelLocs', 'positive', []}; 
        defaults(end+1,:) = {'labelNames', 'cell', {}}; 
        defaults(end+1,:) = {'gain', 'positive', 1}; 
        defaults(end+1,:) = {'imSepColor', 'nonnegative', 2^16}; 
        % defaults for plotting fits
        defaults(end+1,:) = {'fits','freeType',[]};
        defaults(end+1,:) = {'xyzNames','cell',{}}; % {'xc','yc','zc'}
        defaults(end+1,:) = {'nmXYpix','positive',154};
        defaults(end+1,:) = {'nmZpix','positive',100};
        % parse defaults
        pars = ParseVariableArguments(varargin,defaults,mfilename);

   % call figure
   figData.handle = figure(pars.figHandle); clf;
   subplot(1,2,1);
   PlotProjection4D(image4D,'parameters',pars,'projection','xy');
   subplot(1,2,2);
   PlotProjection4D(image4D,'parameters',pars,'projection','xz');

    