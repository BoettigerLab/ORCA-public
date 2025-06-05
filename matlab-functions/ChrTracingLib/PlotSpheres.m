function PlotSpheres(d,varargin)
% plot 3 data as 3D spheres. 
% defaults(end+1,:) = {'subDivisions', 'positive', 30};
% defaults(end+1,:) = {'r', 'positive', 5};
% defaults(end+1,:) = {'lightingOn', 'boolean', true};
% defaults(end+1,:) = {'color', 'colormap', [.3 .3 .3]};
% defaults(end+1,:) = {'alpha', 'nonnegative', 1};


% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'subDivisions', 'positive', 30};
defaults(end+1,:) = {'r', 'positive', 5};
defaults(end+1,:) = {'sphereRadius', 'positive', []};
defaults(end+1,:) = {'lightingOn', 'boolean', true};
defaults(end+1,:) = {'color', 'colormap', [.3 .3 .3]};
defaults(end+1,:) = {'alpha', 'nonnegative', 1};
defaults(end+1,:) = {'view', 'freeType', []};

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

if ~isempty(parameters.sphereRadius) % allow "SphereRadius" input to override earlier behavior 
    parameters.r = parameters.sphereRadius;
end

r = parameters.r;
if numel(r) == 1
    r = repmat(r,size(d,1),1);
end

clr = parameters.color;
if size(clr,1) == 1
    clr = repmat(clr,size(d,1),1);
end

alphas = parameters.alpha;
if size(alphas,1) == 1
    alphas = repmat(alphas,size(d,1),1);
end
    
[sx,sy,sz] = sphere(parameters.subDivisions);

for k =1:size(d,1)
    if ~isnan(d(k,1))
        surf(r(k)*sx+d(k,1),r(k)*sy+d(k,2),r(k)*sz+d(k,3),...
            'EdgeColor','none','FaceColor',clr(k,:),'FaceAlpha',alphas(k)); hold on;
    end
end
   
if ~isempty(parameters.view)
    view(parameters.view);
end

if parameters.lightingOn
    material dull;
    camlight right;
    lighting gouraud;
end

% For fast plotting, take surf out of a loop
% would be best to initialize X,Y,Z as size(repmat(sx,size(d,1))).

% X = [];
% Y = [];
% Z = [];
% for k =1:size(d,1)
%     X = [X; r(k)*sx+d(k,1)];
%     Y = [Y; r(k)*sy+d(k,2)];
%     Z = [Z; r(k)*sz+d(k,3)];
% end
% figure(100); hold on; surf(X,Y,Z,'FaceColor','r','FaceAlpha',.3,'EdgeColor','none');

