function [X,Y,Z,V,surfHandle] = PlotTube(d,varargin)
% 
% PlotTube(d,'interpPts',4,'subDivisions',20,'r',5);
% PlotTube(d,'lightingOn',true,'showInterp','false');
% 
% can have continuous color change
% can have 3D lighting effects
% 
% Requires / wraps matlab function tubeplot.
% allows an additional interpolation of the data. 

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'interpPts', 'positive', 1};
defaults(end+1,:) = {'subDivisions', 'positive', 30};
defaults(end+1,:) = {'r', 'positive', 5};
defaults(end+1,:) = {'lightingOn', 'boolean', true};
defaults(end+1,:) = {'showInterp', 'boolean', false};
defaults(end+1,:) = {'method',{'pchip','spline','linear','nearest','next','previous','cubic','skip'},'pchip'};
defaults(end+1,:) = {'colormap','colormap',[]};
defaults(end+1,:) = {'verbose','boolean',true};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 1
    error('matlabSTORM:invalidArguments', 'a nx3 data array is required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% parameters = ParseVariableArguments([], defaults, mfilename);

%% Main Function

n = length(d);
CS = cat(1,0,cumsum(sqrt(sum(diff(d,[],1).^2,2))));  % From matlab central help. Thanks Sven! http://www.mathworks.com/matlabcentral/answers/21849-how-to-generate-a-3d-spline-curve-can-interp3-be-of-help
% dd = interp1(CS, d, unique([CS(:)' linspace(0,CS(end),parameters.interpPts*n)]),'pchip'); % Interpolate at 100 equally spaced locations from the start to end of your curve. Use 'spline' interpolation. Also, throw in the original points as part of your output since this was requested.  
try
dd = interp1(CS, d, unique([linspace(0,CS(end),parameters.interpPts*n)]),parameters.method); % Interpolate at 100 equally spaced locations from the start to end of your curve. Use 'spline' interpolation. Also, throw in the original points as part of your output since this was requested.


catch er
    if parameters.verbose
        disp(er.message);
        disp('skipping interpolation'); 
    end
    dd = d;
end
if size(dd,1) > 1
    [X,Y,Z,V,surfHandle] = tubeplot(dd(:,1),dd(:,2),dd(:,3),parameters.r,(1:length(dd))',parameters.subDivisions);
else 
    warning('data contains only 1 point');
    return
end

% Some default display options
shading flat; 
if isempty(parameters.colormap)
    colormap(jet(length(dd)));
else
    cmap = GetColorMap(parameters.colormap); 
    if ~strcmp(parameters.method,'skip')
        try
            cmap1 = interp1(CS, cmap(:,1), unique([linspace(0,CS(end),parameters.interpPts*n)]),parameters.method); 
            cmap2 = interp1(CS, cmap(:,2), unique([linspace(0,CS(end),parameters.interpPts*n)]),parameters.method); 
            cmap3 = interp1(CS, cmap(:,3), unique([linspace(0,CS(end),parameters.interpPts*n)]),parameters.method); 
            cmap = [cmap1',cmap2',cmap3'];
            cmap(cmap<0) = 0;
            cmap(cmap>1) = 1;
        catch
        end
    end
    colormap(cmap);
end
if parameters.lightingOn
    material dull;
    camlight left;
    lighting gouraud;
end
   
if parameters.showInterp
    figure; clf; hold on
    plot3(dd(:,1),dd(:,2),dd(:,3),'.r-')
    plot3(d(:,1),d(:,2),d(:,3),'ob-')
    axis image, view(3), legend({'Original','Interp. Spline'})
end
    
    
  