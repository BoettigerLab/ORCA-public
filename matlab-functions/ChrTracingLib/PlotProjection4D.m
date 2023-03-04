function labelOffsets = PlotProjection4D(image4D,varargin)
% Plots a multicolor tile projection of a 4D image. 
% Default projection is xy
% If 'fits' (a table object containing fits.x, fits.y

defaults = cell(0,3); 
% defaults for projection 
defaults(end+1,:) = {'projection',{'xy','xz'},'xy'};
defaults(end+1,:) = {'fixBlack','boolean',true};
defaults(end+1,:) = {'mode',{'subplots','single'},'single'};  % render as subplots (slow and flexible) or as a single image  
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
defaults(end+1,:) = {'MarkerSize','positive',20};
% parse defaults
pars = ParseVariableArguments(varargin,defaults,mfilename);

% Main function
if strcmp(pars.projection,'xy')
    stackDat = squeeze(max(image4D,[],3)); 
else
    stackDat = squeeze(max(permute(image4D,[3,2,1,4]),[],3)); 
end

if pars.fixBlack
    stackDat(stackDat==0) = mode(nonzeros(stackDat(:)));
end





% figure(4); clf; 
fits = pars.fits; % a little unpacking;
if ~isempty(fits)
    if ~ismember('hybe',fits.Properties.VariableNames)
        fits.hybe = (1:height(fits))';
    end
    if ~ismember('idx',fits.Properties.VariableNames)
        fits.idx = ones(height(fits),1);
    end
    if ~isempty(pars.xyzNames)
        fits.x = fits.(pars.xyzNames{1});
        fits.y = fits.(pars.xyzNames{2});
        fits.z = fits.(pars.xyzNames{3});
    end
end

if strcmp(pars.mode,'subplots')
    TileImageStack(stackDat,'parameters',pars);
    if ~isempty(fits)
        tileDatxy = gcf;
        panelsDatxy = flipud(tileDatxy.Children);
        for h=1:length(panelsDatxy)
            if any(strcmp(fits.Properties.VariableNames,'panel'))
                isH = fits.panel == h;   % 
            else
                isH = fits.hybe == h;   
            end
           is1 = fits.idx == 1;
           is2 = fits.idx == 2;
           isN = fits.idx > 2;
           axes(panelsDatxy(h)); hold on;  
           if strcmp(pars.projection,'xy')
            plot(fits.x(isH & is1)/pars.nmXYpix,fits.y(isH & is1)/pars.nmXYpix,'yo','MarkerSize',pars.MarkerSize);
            plot(fits.x(isH & is2)/pars.nmXYpix,fits.y(isH & is2)/pars.nmXYpix,'ys','MarkerSize',pars.MarkerSize);
            plot(fits.x(isH & isN)/pars.nmXYpix,fits.y(isH & isN)/pars.nmXYpix,'y>','MarkerSize',pars.MarkerSize);
           else
            plot(fits.x(isH & is1)/pars.nmXYpix,fits.z(isH & is1)/pars.nmZpix,'yo','MarkerSize',pars.MarkerSize);
            plot(fits.x(isH & is2)/pars.nmXYpix,fits.z(isH & is2)/pars.nmZpix,'ys','MarkerSize',pars.MarkerSize);
            plot(fits.x(isH & isN)/pars.nmXYpix,fits.z(isH & isN)/pars.nmZpix,'y>','MarkerSize',pars.MarkerSize);
           end
        end
    end
else
    if pars.autoContrast
       stackDat = IncreaseContrast(stackDat);
    end
    [~, labelOffsets] = TileImage(stackDat,'parameters',pars,'showImage',true);
    if ~isempty(fits)
        hold on;
        for h=1:size(stackDat,3)
           if any(strcmp(fits.Properties.VariableNames,'panel'))
               isH = fits.panel == h;   % 
           else
               isH = fits.hybe == h;   
           end
           is1 = fits.idx == 1;
           is2 = fits.idx == 2;
           isN = fits.idx > 2;
           if any(strcmp(fits.Properties.VariableNames,'reject'))
               isRj = fits.reject;
           else
               isRj = true(size(isH));
           end
           if strcmp(pars.projection,'xy')
            plot(labelOffsets(h,1) + fits.x(isH & is1)/pars.nmXYpix,labelOffsets(h,2) + fits.y(isH & is1)/pars.nmXYpix,'yo','MarkerSize',pars.MarkerSize);
            plot(labelOffsets(h,1) + fits.x(isH & is2)/pars.nmXYpix,labelOffsets(h,2) + fits.y(isH & is2)/pars.nmXYpix,'ys','MarkerSize',pars.MarkerSize);
            plot(labelOffsets(h,1) + fits.x(isH & isN)/pars.nmXYpix,labelOffsets(h,2) + fits.y(isH & isN)/pars.nmXYpix,'y>','MarkerSize',pars.MarkerSize);
            plot(labelOffsets(h,1) + fits.x(isH & is1 & isRj)/pars.nmXYpix,labelOffsets(h,2) + fits.y(isH & is1 & isRj)/pars.nmXYpix,'ro','MarkerSize',pars.MarkerSize);
            plot(labelOffsets(h,1) + fits.x(isH & is2 & isRj)/pars.nmXYpix,labelOffsets(h,2) + fits.y(isH & is2 & isRj)/pars.nmXYpix,'rs','MarkerSize',pars.MarkerSize);
           else
            plot(labelOffsets(h,1) + fits.x(isH & is1)/pars.nmXYpix,labelOffsets(h,2) + fits.z(isH & is1)/pars.nmZpix,'yo','MarkerSize',pars.MarkerSize);
            plot(labelOffsets(h,1) + fits.x(isH & is2)/pars.nmXYpix,labelOffsets(h,2) + fits.z(isH & is2)/pars.nmZpix,'ys','MarkerSize',pars.MarkerSize);
            plot(labelOffsets(h,1) + fits.x(isH & isN)/pars.nmXYpix,labelOffsets(h,2) + fits.z(isH & isN)/pars.nmZpix,'y>','MarkerSize',pars.MarkerSize);
            plot(labelOffsets(h,1) + fits.x(isH & is1 & isRj)/pars.nmXYpix,labelOffsets(h,2) + fits.z(isH & is1 & isRj)/pars.nmZpix,'ro','MarkerSize',pars.MarkerSize);
            plot(labelOffsets(h,1) + fits.x(isH & is2 & isRj)/pars.nmXYpix,labelOffsets(h,2) + fits.z(isH & is2 & isRj)/pars.nmZpix,'rs','MarkerSize',pars.MarkerSize);
           end
        end
   end 
end