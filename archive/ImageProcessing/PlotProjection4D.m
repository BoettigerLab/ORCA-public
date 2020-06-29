function labelOffsets = PlotProjection4D(image4D,varargin)
% Plots a multicolor tile projection of a 4D image. 
% Default projection is xy
% If 'fits' (a table object containing fits.x, fits.y

defaults = cell(0,3); 
% defaults for projection 
defaults(end+1,:) = {'projection',{'xy','xz'},'xy'};
% defaults for image tile
defaults(end+1,:) = {'mode',{'subplots','single'},'single'};  % render as subplots (slow and flexible) or as a single image  
defaults(end+1,:) = {'numRows', 'positive', 4}; 
defaults(end+1,:) = {'gap', 'nonnegative', 1}; 
defaults(end+1,:) = {'colorTiles', 'boolean', true}; 
defaults(end+1,:) = {'colormap', 'colormap', 'hsv'}; 
defaults(end+1,:) = {'autoContrast', 'boolean', true}; 
defaults(end+1,:) = {'showPanelNumber', 'boolean', true}; 
defaults(end+1,:) = {'fontSize', 'positive', 10}; 
defaults(end+1,:) = {'tileLabels','array',{}};
% defaults for plotting fits
defaults(end+1,:) = {'fits','freeType',[]};
defaults(end+1,:) = {'nmXYpix','positive',154};
defaults(end+1,:) = {'nmZpix','positive',100};
% parse defaults
pars = ParseVariableArguments(varargin,defaults,mfilename);

% Main function
if strcmp(pars.projection,'xy')
    stackDat = squeeze(max(image4D,[],3)); 
else
    stackDat = squeeze(max(permute(image4D,[3,2,1,4]),[],3)); 
end

% figure(4); clf; 
fits = pars.fits; % a little unpacking;
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
           axes(panelsDatxy(h)); hold on;  
           if strcmp(pars.projection,'xy')
            plot(fits.x(isH & is1)/pars.nmXYpix,fits.y(isH & is1)/pars.nmXYpix,'yo','MarkerSize',15);
            plot(fits.x(isH & is2)/pars.nmXYpix,fits.y(isH & is2)/pars.nmXYpix,'ys','MarkerSize',15);
           else
            plot(fits.x(isH & is1)/pars.nmXYpix,fits.z(isH & is1)/pars.nmZpix,'yo','MarkerSize',15);
            plot(fits.x(isH & is2)/pars.nmXYpix,fits.z(isH & is2)/pars.nmZpix,'ys','MarkerSize',15);
           end
        end
    end
else
    
    [~, labelOffsets] = TileImage(IncreaseContrast(stackDat),'parameters',pars,'showImage',true);
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
           if strcmp(pars.projection,'xy')
            plot(labelOffsets(h,1) + fits.x(isH & is1)/pars.nmXYpix,labelOffsets(h,2) + fits.y(isH & is1)/pars.nmXYpix,'yo','MarkerSize',15);
            plot(labelOffsets(h,1) + fits.x(isH & is2)/pars.nmXYpix,labelOffsets(h,2) + fits.y(isH & is2)/pars.nmXYpix,'ys','MarkerSize',15);
           else
            plot(labelOffsets(h,1) + fits.x(isH & is1)/pars.nmXYpix,labelOffsets(h,2) + fits.z(isH & is1)/pars.nmZpix,'yo','MarkerSize',15);
            plot(labelOffsets(h,1) + fits.x(isH & is2)/pars.nmXYpix,labelOffsets(h,2) + fits.z(isH & is2)/pars.nmZpix,'ys','MarkerSize',15);
           end
        end
   end 
end