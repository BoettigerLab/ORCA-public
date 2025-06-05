function cellBorders = FindNucBorders(im2D,cntrs,radii,varargin)
% cellBorders = FindNucBorders(im2D,cntrs,radii,$NAME$,$VALUE$)
% 
% this spins the image around the cell centers, to refine the borders
% identified with a FindCircles approach based on the max derivative. 

defaults = cell(0,3);
defaults(end+1,:) = {'nAngles','integer',21};
defaults(end+1,:) = {'segLength','integer',20};
defaults(end+1,:) = {'smoothFindEdge','integer',3};
defaults(end+1,:) = {'smoothFinalBorder','integer',3};
defaults(end+1,:) = {'imageFilter','array',fspecial('gaussian',10,20)};
defaults(end+1,:) = {'showOverlayFig','integer',3};
defaults(end+1,:) = {'showExtraFig','integer',0};
defaults(end+1,:) = {'dilate','integer',1};
defaults(end+1,:) = {'xyNucWindow','positive',2.25}; % search window is the nuc radius x this muliplier
defaults(end+1,:) = {'convex','boolean',false};
defaults(end+1,:) = {'showNum','boolean',true};
pars = ParseVariableArguments(varargin,defaults,mfilename);

nCells = length(radii);
fs = pars.imageFilter;
cmap = hsv(round(1.1*nCells));
cellBorders = cell(nCells,1);
if pars.showOverlayFig > 0
    figure(pars.showOverlayFig); clf; 
    imagesc(IncreaseContrast(im2D,'high',.999)); colormap(gray);
    hold on;
end

buf = round(max(pars.xyNucWindow*round(radii)));
imD_pad = padarray(im2D,[buf,buf],0,'both');

for c= 1:nCells 
    xc1 = buf+ floor( cntrs(c,1)- pars.xyNucWindow*radii(c))+1;
    xc2 = buf+ floor( cntrs(c,1)+ pars.xyNucWindow*radii(c))  ;
    yc1 = buf+ floor( cntrs(c,2)- pars.xyNucWindow*radii(c))+1;
    yc2 = buf+ floor( cntrs(c,2)+ pars.xyNucWindow*radii(c))  ;
    nucCrop = imD_pad(yc1:yc2,xc1:xc2);
    im2 = imfilter(nucCrop,fs);

    % start spinning!
    [ny,nx] = size(im2); 
    mx = round(nx/2);
    my = round(ny/2);
    thetas = linspace(0,360,pars.nAngles);
    cellPoly = zeros(pars.nAngles+1,2);
    for n=1:pars.nAngles % n=4
        imRot = imrotate(im2,thetas(n),'crop'); % the lazy way to do this is to rotate the image. 
        x1 = round(mx+radii(c)-pars.segLength/2);
        x2 = min(round(mx+radii(c)+pars.segLength/2),nx);
        segPix = imRot(my,x1:x2);   
        segPixS = smooth(double(segPix),pars.smoothFindEdge);
        [~,idx] = min(diff(segPixS));
        % rotate back the xy data
        xr = round(radii(c)-pars.segLength/2+idx+1+pars.dilate);
        cellPoly(n,1) = mx + xr*cosd(thetas(n));
        cellPoly(n,2) = my + xr*sind(thetas(n));
        if pars.showExtraFig > 0
            figure(pars.showExtraFig); clf; 
            subplot(1,2,1); imagesc(imRot);
            hold on; plot([x1,x2],[my,my],'y-');  
            subplot(1,2,2); plot(segPix); 
            hold on; plot(segPixS);
            plot(1+idx,segPixS(1+idx),'r.');
        end
    end
    cellPoly(end,:) = cellPoly(1,:); % close the loop
    
    cellPoly(:,1) = cellPoly(:,1) + xc1 - 1 - buf;
    cellPoly(:,2) = cellPoly(:,2) + yc1 - 1 - buf;
    
    % show raw data
    if pars.showOverlayFig > 0
        figure(pars.showOverlayFig); hold on;
        plot(cellPoly(:,1),cellPoly(:,2),'.','color',cmap(c,:)); % only plot the pure z-projection outline on the z-projection data.
    end
    
    % === optional adustments to outlines
    % - smoothed, average immideate neighbors
    if pars.smoothFinalBorder > 1
        cellPoly(:,1) = smooth(cellPoly(:,1),pars.smoothFinalBorder);
        cellPoly(:,2) = smooth(cellPoly(:,2),pars.smoothFinalBorder);
    end
    % - convex hull   
    if pars.convex
        conv = convhull(cellPoly);
        cellPoly = cellPoly(conv,:);
    end
    
    cellBorders{c} = cellPoly;
    
    % -- plot final results
    if pars.showOverlayFig > 0
        figure(pars.showOverlayFig);  
        plot(cellPoly(:,1),cellPoly(:,2),'-','color',cmap(c,:)); % only plot the pure z-projection outline on the z-projection data.
        if pars.showNum
            text(cellPoly(1,1)-pars.segLength,cellPoly(1,2),num2str(c),'color',cmap(c,:))
        end
    end   
    
    if pars.showExtraFig > 0
        figure(pars.showExtraFig); clf; imagesc(xy); hold on;
        viscircles(cntrs(c,:)-[xc1,yc1] +[buf,buf],radii(c),'color','r');
        plot(cellPoly(:,1),cellPoly(:,2),'r.');
    end
end