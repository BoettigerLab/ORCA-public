function tubeFig = PlotPolymerTube(xyz,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'fillmissing','boolean',true}; % 
defaults(end+1,:) = {'shortestPath','boolean',false}; 
defaults(end+1,:) = {'showInterp','boolean',false}; 
defaults(end+1,:) = {'tubeRadius','positive',20};
defaults(end+1,:) = {'sphereRadius','positive',30};
defaults(end+1,:) = {'showSpheres','boolean',true};
defaults(end+1,:) = {'showTube','boolean',true};
defaults(end+1,:) = {'w','positive',1200}; % box size
defaults(end+1,:) = {'numColors','positive',[]}; % box size
defaults(end+1,:) = {'cent','array',[]}; % box size
defaults(end+1,:) = {'colormap','colormap','hsvCut'}; 
defaults(end+1,:) = {'showLabels','boolean',false}; 
defaults(end+1,:) = {'labels','cell',{}}; 
defaults(end+1,:) = {'alpha','fraction',.5}; 
defaults(end+1,:) = {'lightOn','boolean',true}; 
defaults(end+1,:) = {'fit','boolean',true};
defaults(end+1,:) = {'fontSize','positive',20};
defaults(end+1,:) = {'autoAxisLabels','boolean',false};
defaults(end+1,:) = {'method',{'pchip','spline','skip'},'skip'};
defaults(end+1,:) = {'center','boolean',true}; 
defaults(end+1,:) = {'interpPts','positive',10}; 
defaults(end+1,:) = {'applyColorMap','boolean',true};  
defaults(end+1,:) = {'number','boolean',false}; % obsolete 
defaults(end+1,:) = {'view','freeType',[]};
% parameters for removeJumps
defaults(end+1,:) = {'maxJump','positive',inf}; % box size
defaults(end+1,:) = {'localRegion','integer',4}; % number of points in front and behind to use for estimating expected position 
defaults(end+1,:) = {'maxDiff','positive',2}; % fold change greater than the median step size of other points from their local area which is acceptable
defaults(end+1,:) = {'maxAbsStep','positive',inf}; % maximum absolute step size
defaults(end+1,:) = {'removeLoners','boolean',false}; % if a point is the only one within its local region, drop it. 

pars = ParseVariableArguments(varargin,defaults,mfilename);

w= pars.w; 
tubeFig = gcf;

if pars.center
    xyz = CenterPolymer(xyz);
end

if isempty(xyz)
    warning('data is empty'); 
    return
end

if sum(isnan(xyz(:,1)),1) == size(xyz,1)
    warning('data is all nan');
    return
else

    % interpolate spline data
    % xyz_raw has nans
    % xyz     is interpolated, no NaNs
    % polyData 
    xyz_raw = xyz;
    polyData = xyz;
    % truncate out NaNs (not sure this is desired)
    skip = isnan(xyz(:,1));
    polyData(skip,:) = []; 

    if ~isinf(pars.maxJump)    
        xyz = RemoveJumps(xyz,'maxAbsStep',pars.maxJump,'localRegion',pars.localRegion,'maxDiff',pars.maxDiff);
        xyz_raw = xyz;
        % jumps = nanmean(squareform(pdist(xyz))) > pars.maxJump;% [false; sqrt(sum(diff(xyz,1,1).^2,2)) > pars.maxJump]; % not ideal, both the step out and the step back are too big 
        % xyz(jumps,:) = nan;
        % xyz = fillmissing(xyz,'linear');
        % xyz_raw(jumps,:) = nan; %
    end

    % fill missing (only done for the tube, xyz, not xyz_raw)
    if pars.fillmissing
        xyz(1,:) = polyData(1,:)-.1; 
        xyz(end,:) = polyData(end,:)+.1;
        xyz = fillmissing(xyz,'linear');
    end


    % setup colormap
    if isempty(pars.numColors)
        numHybes = size(xyz_raw,1);
    else
        numHybes = pars.numColors;
    end
    cmap = GetColorMap(pars.colormap,numHybes); 
    cmapOut = cmap;
    % cmapOut(skip,:) = [];
    
   % plot    
    if pars.showTube && ~pars.shortestPath
        PlotTube(xyz,'r',pars.tubeRadius,'interpPts',pars.interpPts,'method',pars.method,...
                 'colormap',cmapOut,'lightingOn',false,'verbose',false,'applyColorMap',pars.applyColorMap,'showInterp',pars.showInterp); hold on;   
        set(gca,'color','w'); 
        shading flat;
        hold on;  
        alpha(pars.alpha);  
        if pars.showSpheres
            if pars.showTube
                freezeColors; 
            end
        end
    end
    
    if pars.showSpheres
        PlotSpheres(xyz_raw,'r',pars.sphereRadius,'color',cmap,'alpha',pars.alpha,'lightingOn',pars.lightOn);
    end
    
 
    
    
    % add labels on spheres (optional); 
    if pars.showLabels
        if isempty(pars.labels)
            labels = cellstr(num2str( (1:size(xyz,1))' ))';
        else
            labels = pars.labels;
        end
        text(xyz_raw(:,1),xyz_raw(:,2),xyz_raw(:,3),labels,'color','w');
    end
    if isempty(pars.cent)
        cent = nanmean(xyz,1);
    else
        cent = pars.cent;
    end
    

    if ~isempty(pars.view)
        view(pars.view);
    end
    
    if pars.lightOn
        material dull;
        camlight left;
        lighting gouraud;
    end
    
    tubeFig.Color = 'w';
    
    if pars.shortestPath && pars.showTube
        hold on;
        x = polyData(:,1); y=polyData(:,2); z=polyData(:,3); 
        col=permute(cat(3,cmapOut,cmapOut),[3,1,2]);
        surface([x,x]',[y,y]',[z,z]',col,...
            'facecol','no',...
            'edgecol', 'flat',... %
            'linew',pars.tubeRadius); 
        set(gca,'color','w');
        if pars.applyColorMap
            colormap(cmapOut);
        end
    end
    

    
end

if pars.fit
    axis tight;
else
     xlim([cent(1)-w,cent(1)+w]);  
     ylim([cent(2)-w,cent(2)+w]);  
     zlim([cent(3)-w*.8,cent(3)+w*.8]); 
    view(-23,8);   
end
if pars.autoAxisLabels
    xl = 100*round(get(gca,'Xlim')/100)+[50,-50];
    yl = 100*round(get(gca,'Ylim')/100)+[50,-50];
    zl = 100*round(get(gca,'Zlim')/100)+[50,-50];
    xl2 = xl - xl(1) ; 
    yl2 = yl - yl(1) ;
    zl2 = zl - zl(1) ;
    set(gca,'FontSize',pars.fontSize,'Xtick',xl,'XtickLabel',xl2,...
        'Ytick',yl,'YtickLabel',yl2,'Ztick',zl,'ZtickLabel',zl2);
end
