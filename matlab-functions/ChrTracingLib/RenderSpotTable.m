function [imOut,ax,ay] = RenderSpotTable(datSpots,varargin)
% [imOut,ax,ay] = RenderSpotTable(datSpots)
% 
%% Inputs
% datSpots 
% a table containing variable names: x, y, z, and hybe, which are the
% position coordinates (x,y,z) for each spot and the hybe in which it was
% detected.  Ideally also contains upper and lower bounds for x,y,z
% positions, labelled 'xL', 'xU', 'yL', 'yU', 'zL' and 'zU'.  
%
%% Outputs
% 

%%
defaults = cell(0,3);
defaults(end+1,:) = {'gain','positive',1};
defaults(end+1,:) = {'boxWidth','integer',1000};
defaults(end+1,:) = {'showPlots','boolean',true}; % supress display
defaults(end+1,:) = {'vecLength','positive',100}; % length of x,y,z coordinate axis 
defaults(end+1,:) = {'vecOffset','positive',100}; % offset from UL of the x,y,z coordinate axis
defaults(end+1,:) = {'scale','positive',.5}; % 
defaults(end+1,:) = {'border','integer',100};
defaults(end+1,:) = {'spotWidth','integer',100};
defaults(end+1,:) = {'numColors','integer',[]};
defaults(end+1,:) = {'recenter','array',[500;500;0]};
defaults(end+1,:) = {'thetaX','float',15}; % rotation around x in degrees. 
defaults(end+1,:) = {'thetaY','float',-10}; % rotation around y in degrees. 
defaults(end+1,:) = {'thetaZ','float',0}; % rotation around z in degrees. 
defaults(end+1,:) = {'err','nonnegative',0}; % additional positional uncertainties to propagate
defaults(end+1,:) = {'showHybeNum','boolean',true};
defaults(end+1,:) = {'verbose','boolean',true}; 
defaults(end+1,:) = {'veryverbose','boolean',false}; 
parameters = ParseVariableArguments(varargin,defaults,mfilename);
% parameters = ParseVariableArguments([],defaults,mfilename);
%

% datSpots = spotTable(~logical(spotTable.isfid),:);

% short hand parsing of optional parameters

recenter =parameters.recenter;
boxWidth = parameters.boxWidth;
sc = parameters.scale;
border = parameters.border;
tx = parameters.thetaX; % rotate degrees around x
ty = parameters.thetaY; % rotate degrees around y
tz = parameters.thetaZ; % rotate degrees around z 
err = parameters.err; 

%% compute some more parameters and initialize variables
numSpots = height(datSpots);
numHybes = max(datSpots.hybe);
cmap = jet(numHybes); 
imSize = sc*boxWidth+border;
imOut = zeros(imSize,imSize,numHybes,'uint16');
imH = zeros(numSpots,1,'uint16'); 
xyz_cols = StringFind(datSpots.Properties.VariableNames,{'x','y','z'},'exactly',true);
spotCenter = mean(datSpots{:,xyz_cols},1)';
xys = NaN(numSpots,2);
warn = false; 
warn2 = false;
for i=1:numSpots
    xyz = datSpots{i,xyz_cols}' - spotCenter;
    h = datSpots.hybe(i);
    [xyz,axisVecs] = RotateVector3D(xyz,tx,ty,tz,[0;0;0]); 
    xyz = xyz + recenter;
    xys(i,:) = sc*[xyz(1),xyz(2)]+ border/2;
    % determine spot width based on uncertainty;
    % a few metrics are allowed, which are tested in order of preference
    try
        w = sc*mean([(datSpots.xU(i) - datSpots.xL(i)),(datSpots.yU(i) - datSpots.yL(i))]);
    catch
        try
            w= sc/2*mean([datSpots.wx(i),datSpots.wy(i)]);
        catch
            w = 2;
            if parameters.verbose
               cprintf([1 .5 0],'warning: did not find confidence limits or spot-width data');  
            end
        end
    end
        
    w = sqrt(w^2 + (sc*err)^2);
    w = min(imSize/10,w); 
    spotWidth = floor(6*w);
    im = uint16(10000*sqrt(datSpots.h(i))*fspecial('gaussian',spotWidth,w));
    x1 = floor(sc*xyz(1) - spotWidth/2) + border/2;
    x2 = floor(sc*xyz(1) + spotWidth/2) + border/2;
    y1 = floor(sc*xyz(2) - spotWidth/2) + border/2;
    y2 = floor(sc*xyz(2) + spotWidth/2) + border/2;
    % if image is off the edge, just trim it. 
    [rows,cols] = size(im);
    imx1 = 1;
    imx2 = cols;
    imy1 = 1;
    imy2 = rows;
    if x1 < 0
        imx1 = -x1+1;
        x1 = 0;
        warn = true;
    end
    if y1 < 0
        imy1 = -y1+1;
        y1 = 0;
        warn = true;
    end
    if x2 > imSize
        imx2 = cols - (x2-imSize);
        x2 = imSize;
        warn = true;
    end
    if y2 > imSize
        imy2 = rows - (y2-imSize);
        y2 = imSize;
        warn = true;
    end
    if imx1 > cols || imx2 < 0 || imy1 > rows || imy2 < 0
       warn2 = true; 
    end
    imOut(y1+1:y2,x1+1:x2,h) =  imOut(y1+1:y2,x1+1:x2,h) + im(imy1:imy2,imx1:imx2);
    imH(i) = max(im(:));
end

if warn && parameters.veryverbose
   cprintf([1 .5 0], 'warning: at least 1 spot is partially off screen'); 
   cprintf([1 .5 0], 'Try increasing "boxWidth" or decreasing "scale".');
end
if warn2 && parameters.verbose
    cprintf([1 .5 0], 'warning: at least 1 spot is completely off screen'); 
    cprintf([1 .5 0], 'Try increasing "boxWidth" or decreasing "scale".');
end
    

c = 2^16/max(imH);
imOut = imOut*c;
ax = cat(1,zeros(1,3),axisVecs(1,:))*parameters.vecLength*sc+parameters.vecOffset*sc;
ay = cat(1,zeros(1,3),axisVecs(2,:))*parameters.vecLength*sc+parameters.vecOffset*sc;
if parameters.showPlots
    Ncolor(parameters.gain*imOut,cmap); hold on;
    colormap(cmap); colorbar;   
    plot(ax,ay,'w');
    text(ax(2,:)+2,ay(2,:)-2,{'x','y','z'},'color','w');
    if parameters.showHybeNum
        hNums = cellstr(num2str( (1:numSpots)' ));
       text(xys(:,1),xys(:,2),hNums,'color','w');
    end
end
