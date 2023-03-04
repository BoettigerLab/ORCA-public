function [cellBorders,cellMesh] = FindNucBorders3D(im3D,cntrs,radii,varargin)


defaults = cell(0,3);
defaults(end+1,:) = {'nAngles','integer',21};
defaults(end+1,:) = {'segLength','integer',20};
defaults(end+1,:) = {'smoothFindEdge','integer',5};
defaults(end+1,:) = {'smoothFinalBorder','integer',3};
defaults(end+1,:) = {'imageFilter','array',fspecial('gaussian',8,16)};
defaults(end+1,:) = {'showRotFig','integer',0};
defaults(end+1,:) = {'showProjFig','integer',0};
defaults(end+1,:) = {'showOverlayFig','integer',3};
defaults(end+1,:) = {'show3DFig','integer',4};
defaults(end+1,:) = {'showTilesFig','integer',10};
defaults(end+1,:) = {'xyNucWindow','positive',2.25}; % search window is the nuc radius x this muliplier
defaults(end+1,:) = {'convex','boolean',false};
defaults(end+1,:) = {'yTiles','integer',4};
defaults(end+1,:) = {'nmXY','integer',154};
defaults(end+1,:) = {'nmZ','integer',100};
defaults(end+1,:) = {'scope',{'scope1','scope2','scope3','other','auto'},'scope1'};
pars = ParseVariableArguments(varargin,defaults,mfilename);


fs =  pars.imageFilter; % fspecial('gaussian',10,20);
segLen = round(pars.segLength*pars.nmXY/pars.nmZ);
nCells = length(radii);
cmap = hsv(round(1.1*nCells));
cellBorders = cell(nCells,1);
yTiles = pars.yTiles;
xTiles = ceil(nCells/yTiles);

if pars.showOverlayFig>0 % show plot
    figure(pars.showOverlayFig); clf; 
    im_max = max(im3D,[],3);
    imagesc(IncreaseContrast(im_max,'high',.999)); colormap(gray);
end
if pars.showTilesFig > 0
   figure(pars.showTilesFig); clf;
   figure(pars.showTilesFig+1); clf;
end
% pad image to avoid out of range errors
buf = round(max(pars.xyNucWindow*radii));
im3D_pad = padarray(im3D,[buf,buf],0,'both');

% - loop over all cells and find borders
for c= 1:nCells 
    % crop nucleus
    xc1 = buf+ floor( cntrs(c,1)- pars.xyNucWindow*radii(c))+1;
    xc2 = buf+ floor( cntrs(c,1)+ pars.xyNucWindow*radii(c))  ;
    yc1 = buf+ floor( cntrs(c,2)- pars.xyNucWindow*radii(c))+1;
    yc2 = buf+ floor( cntrs(c,2)+ pars.xyNucWindow*radii(c))  ;
    nucCrop = im3D_pad(yc1:yc2,xc1:xc2,:);

    % isotropically sample (scope1, 154 nm-xy, 100 nm)
    crop3d = imresize(nucCrop,pars.nmXY/pars.nmZ); % don't change z, downsample xy
    thetas = linspace(0,360,pars.nAngles);
    phis = linspace(0,180,pars.nAngles);   
    [nx,ny,nz] = size(crop3d); 
    imx = round(nx/2);
    imy = round(ny/2);
    imz = round(nz/2);
    cellPoly = zeros(pars.nAngles,pars.nAngles+1,3);
    for m=1:pars.nAngles % m=8
        % rotate image
        nuc3D = RotateMatrix3D(crop3d,phis(m),0); 
        xy = max(nuc3D,[],3);   
        im2 = imfilter(xy,fs);  
        [ny,nx] = size(xy);
        mx = round(nx/2);
        my = round(ny/2);
        tempPoly = zeros(pars.nAngles,2);
        for n=1:pars.nAngles % n=1
            imRot = imrotate(im2,thetas(n),'crop'); % the lazy way to do this is to rotate the image. 
            x1 = round(mx+radii(c)-segLen/2);
            x2 = min(round(mx+radii(c)+segLen/2),nx);
            segPix = imRot(my,x1:x2);   
            segPixS = smooth(double(segPix),pars.smoothFindEdge);
            [~,idx] = min(diff(segPixS));
            % rotate back the xy data
            xr = round(radii(c)-segLen/2+idx+1);
            % circle in current coordinates
            xc = mx + xr*cosd(thetas(n));
            yc = my + xr*sind(thetas(n));
            tempPoly(n,:) = [xc,yc];
            % circle rotated back into original coordinates
            cellPoly(m,n,1) = imx + xr*cosd(thetas(n));
            cellPoly(m,n,2) = imy + xr*sind(thetas(n))*cosd(phis(m)); %  yc*cosd(phis(m));
            cellPoly(m,n,3) = imz + (yc-my)*sind(phis(m));
            % cellPoly(m,n,1) = imx + xr*cosd(thetas(n));
            % cellPoly(m,n,2) = imy + xr*sind(thetas(n))*cosd(phis(m));
            % cellPoly(m,n,3) = imz + xr*sind(thetas(n))*sind(phis(m));
            if pars.showRotFig > 0
                figure(pars.showRotFig); clf; 
                subplot(1,2,1); imagesc(imRot);
                hold on; plot([x1,x2],[my,my],'y-');  
                subplot(1,2,2); plot(segPix); 
                hold on; plot(segPixS);
                plot(1+idx,segPixS(1+idx),'r.');
            end
        end
        cellPoly(m,n+1,:) = cellPoly(m,n,:); % close loop 
        if pars.showProjFig > 0
            nuc3D = RotateMatrix3D(crop3d,phis(m),0); 
            rot2D = max(nuc3D,[],3);
            [xy,xz] = ProjectIm3D(crop3d); % 
            figure(pars.showProjFig); clf; 
            subplot(1,3,1); imagesc(rot2D); % in current coords
            hold on; plot(tempPoly(:,1),tempPoly(:,2),'y-');  
            title('rotated coordinates');
            subplot(1,3,2); imagesc(xy); 
            hold on; plot(cellPoly(m,:,1),cellPoly(m,:,2),'y-');
            title('orig coordinates, xy');
            subplot(1,3,3); imagesc(xz); 
            hold on; plot(cellPoly(m,:,1),cellPoly(m,:,3),'y-');
            title('orig coordinates, xz');
            colormap(gray);
        end
        
        % === optional adustments to outlines
        % - smoothed, average immideate neighbors
        if pars.smoothFinalBorder > 1
            cellPoly(m,:,1) = smooth(cellPoly(m,:,1),pars.smoothFinalBorder);
            cellPoly(m,:,2) = smooth(cellPoly(m,:,2),pars.smoothFinalBorder);
        end
        % - convex hull   
        if pars.convex
            conv = convhull(squeeze(cellPoly(m,:,:)));
            cellPoly(m,:,:) = cellPoly(m,conv,:);
        end
    end
    
    
    % --- plot outline on cells
    %  note, these plots are overlayed on cubic voxels
    if pars.showTilesFig > 0
        %-- xy-proj
        [xy,xz] = ProjectIm3D(crop3d);
        figure(pars.showTilesFig); 
        subplot(yTiles,xTiles,c); imagesc(xy); hold on;
        plot(squeeze(cellPoly(1,:,1)),squeeze(cellPoly(1,:,2)),'.-','color',cmap(c,:));
        title(c); 
        colormap(gray);
        % -- near xz proj
        a = floor((pars.nAngles+1)/2);    
        figure(pars.showTilesFig+1);
        subplot(yTiles,xTiles,c); imagesc(xz); hold on;
        plot(squeeze(cellPoly(a,:,1)),squeeze(cellPoly(a,:,3)),'.-','color',cmap(c,:));
        title(c);
        colormap(gray);
    end
    
    % convert back into FOV coordinates and pixels
    cellPoly(:,:,1) = pars.nmZ/pars.nmXY*cellPoly(:,:,1) + xc1 - 1 -buf;
    cellPoly(:,:,2) = pars.nmZ/pars.nmXY*cellPoly(:,:,2) + yc1 - 1 -buf;
    
    cellBorders{c} = cellPoly;
    
    % plot outlines on image (original, non-cubic voxels)
    if pars.showOverlayFig > 0
        figure(pars.showOverlayFig); hold on;
        plot(cellPoly(1,:,1),cellPoly(1,:,2),'.-','color',cmap(c,:)); % only plot the pure z-projection outline on the z-projection data.
        hold on; text(cellPoly(1,1,1)-20,cellPoly(1,1,2),num2str(c),'color','c');
    end
    
    if pars.show3DFig > 0 || nargout > 1
        x = cellPoly(:,:,1); x=x(:);
        y = cellPoly(:,:,2); y=y(:);
        z = cellPoly(:,:,3); z=z(:);
        dt = delaunay(x,y,z);
        % t = tsearchn(cellPoly,dt,spots) % convex
        if pars.show3DFig > 0
            figure(pars.show3DFig); clf;
            trisurf(dt,x,y,z,'FaceAlpha',.3);
            trisurf(dt,x,y,z,'EdgeColor','none','FaceColor',[.1,.3,1],'FaceAlpha',.3);
            material dull;
            camlight left;
        end
    end
    pause(.1);
end