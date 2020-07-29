function PlotPolymer(B,varargin)
% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'nodeColor', 'colormap', [0 0 0]};
defaults(end+1,:) = {'linkColor', 'colormap', [0 0 1]};
defaults(end+1,:) = {'bondColor', 'colormap', [1 0 0]};
defaults(end+1,:) = {'nodeSize', 'nonnegative', 10};
defaults(end+1,:) = {'linkWidth', 'nonnegative', 2};
defaults(end+1,:) = {'bondWidth', 'nonnegative', 1};
defaults(end+1,:) = {'labels', 'array', []};
defaults(end+1,:) = {'showNodes', 'boolean', true};
defaults(end+1,:) = {'showLinks', 'boolean', true};
defaults(end+1,:) = {'showBonds', 'boolean', true};
defaults(end+1,:) = {'Attached', 'array', []};
defaults(end+1,:) = {'setLims','boolean',false};

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


N = size(B,1);
if parameters.showNodes
    plot3(B(:,1),B(:,2),B(:,3),'.','color',parameters.nodeColor,'MarkerSize',parameters.nodeSize);
    hold on;
end

if ~isempty(parameters.labels)
    text(B(:,1),B(:,2),B(:,3),parameters.labels,'color',parameters.nodeColor);
end

if parameters.showLinks && size(parameters.linkColor,1) < 2
    plot3(B(:,1),B(:,2),B(:,3),'color',parameters.linkColor,'LineWidth',parameters.linkWidth); hold on;
else
    hasData = find(~isnan(B(:,1)));
    for n=1:length(hasData)-1      
        i1 = hasData(n);
        i2 = hasData(n+1);
        c = floor(mean([i1,i2]));
        plot3(B([i1,i2],1),B([i1,i2],2),B([i1,i2],3),'color',parameters.linkColor(c,:),'LineWidth',parameters.linkWidth);
        hold on;       
    end    
end
  
if parameters.setLims
    xlim([min(B(:,1)),max(B(:,1))]);
    ylim([min(B(:,2)),max(B(:,2))]);
    zlim([min(B(:,3)),max(B(:,3))]);
end
if parameters.showBonds && ~isempty(parameters.Attached)
    for i=1:N
        plot3([B(i,1); B(parameters.Attached(i,:),1)],...
              [B(i,2); B(parameters.Attached(i,:),2)],...
              [B(i,3); B(parameters.Attached(i,:),3)],...
              'color',parameters.bondColor,'LineWidth',parameters.bondWidth);
    end
end