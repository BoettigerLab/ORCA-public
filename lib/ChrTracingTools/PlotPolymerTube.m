function tubeFig = PlotPolymerTube(xyz,varargin)

defaults = cell(0,3);
defaults(end+1,:) = {'fillmissing','boolean',true};
defaults(end+1,:) = {'tubeRadius','positive',20};
defaults(end+1,:) = {'sphereRadius','positive',30};
defaults(end+1,:) = {'showSpheres','boolean',true};
defaults(end+1,:) = {'showTube','boolean',true};
defaults(end+1,:) = {'w','positive',1200}; % box size
defaults(end+1,:) = {'maxJump','positive',inf}; % box size
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
defaults(end+1,:) = {'method',{'spline','skip'},'skip'};
defaults(end+1,:) = {'center','boolean',true};
defaults(end+1,:) = {'number','boolean',false}; % obsolete
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
    xyz_raw = xyz;
    polyData = xyz(~isnan(xyz(:,1)),:); 
    xyz(1,:) = polyData(1,:)-.1; 
    xyz(end,:) = polyData(end,:)+.1;
    xyz = fillmissing(xyz,'linear');
    jumps = [false; sqrt(sum(diff(xyz,1,1).^2,2)) > pars.maxJump];
    xyz(jumps,:) = nan;
    xyz = fillmissing(xyz,'linear');
    xyz_raw(jumps,:) = nan;
    
    % setup colormap
    if isempty(pars.numColors)
        numHybes = size(xyz_raw,1);
    else
        numHybes = pars.numColors;
    end
    cmap = GetColorMap(pars.colormap,numHybes); 
    cmapOut = cmap;
    cmapOut(isnan(xyz(:,1)),:) = [];
    
    % plot 
    if pars.showTube
        PlotTube(xyz,'r',pars.tubeRadius,'interpPts',10,'method',pars.method,...
                 'colormap',cmapOut,'lightingOn',false,'verbose',false); hold on;   
        set(gca,'color','k'); 
        shading flat;
        hold on;  
        alpha(pars.alpha); 
        
    end
    if pars.showSpheres
        if pars.showTube
            freezeColors; 
        end
        PlotSpheres(xyz_raw,'r',pars.sphereRadius,'color',cmap,'alpha',pars.alpha);
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
    
    material dull;
    if pars.lightOn
        camlight left;
    end
    lighting gouraud;
    tubeFig.Color = 'w';
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
