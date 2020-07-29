function PlotSpheres(d,varargin)
% plot 3 data as 3D spheres. 

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'subDivisions', 'positive', 30};
defaults(end+1,:) = {'r', 'positive', 5};
defaults(end+1,:) = {'lightingOn', 'boolean', true};
defaults(end+1,:) = {'colormap', 'colormap', 'jet'};

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

%% Main Function

% parameters.lightingOn = true;
% parameter.subDivisions = 30; 
% parameters.r = 10; 
% parameters.colormap = 'jet';

r = parameters.r;
[sx,sy,sz] = sphere(parameters.subDivisions);

try
    clr = eval([parameters.colormap,'(size(d,1))']);
    for k =1:length(d); 
        surf(r*sx+d(k,1),r*sy+d(k,2),r*sz+d(k,3),...
            'EdgeColor','none','FaceColor',clr(k,:)); hold on;
    end
catch
    clr = parameters.colormap; 
    for k =1:length(d); 
        surf(r*sx+d(k,1),r*sy+d(k,2),r*sz+d(k,3),...
            'EdgeColor','none','FaceColor',clr); hold on;
    end
end
    
if parameters.lightingOn
    material dull;
    camlight left;
    lighting gouraud;
end
