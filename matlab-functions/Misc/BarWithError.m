function handles = BarWithError(barvalues, errors, varargin)
%
%
% BarWithError is the m-by-n matrix of barvalues to be plotted.
% BarWithError calls the MATLAB bar function and plots m groups of n bars using the width and bw_colormap parameters.
% If you want all the bars to be the same color, then set bw_colormap equal to the RBG matrix value ie. (bw_colormap = [1 0 0] for all red bars)
% BarWithError then calls the MATLAB errorbar function to draw barvalues with error bars of length error.
% groupnames is an m-length cellstr vector of groupnames (i.e. groupnames = {'group 1'; 'group 2'}).  For no groupnames, enter [] or {}
% The errors matrix is of the same form of the barvalues matrix, namely m group of n errors.
% Gridstatus is either 'x','xy', 'y', or 'none' for no grid.
% No legend will be shown if the legend paramter is not provided
% 'error_sides = 2' plots +/- std while 'error_sides = 1' plots just + std
% legend_type = 'axis' produces the legend along the x-axis while legend_type = 'plot' produces the standard legend.  See figure for more details
%
%--------------------------------------------------------------------------
% Variable Inputs (Flag/ data type /(default)):
%--------------------------------------------------------------------------
% The following default values are used if parameters are left out.
% width = 1 (0 < width < 1; widths greater than 1 will produce overlapping bars)
% groupnames = '1', '2', ... number_of_groups
% bw_title, bw_xlabel, bw_ylabel = []
% bw_color_map = jet
% gridstatus = 'none'
% bw_legend = []
% error_sides = 2;
% legend_type = 'plot';
%
% A list of handles are returned so that the user can change the properties of the plot
% handles.ax: handle to current axis
% handles.bars: handle to bar plot
% handles.errors: a vector of handles to the error plots, with each handle corresponding to a column in the error matrix
% handles.legend: handle to legend
%
%
% See the MATLAB functions bar and errorbar for more information
%
% Alistair Boettiger (adapted from code by Bolu Ajiboye)
% boettiger.alistair@gmail.com
% July 01 2014
% 
% Original BarWithError author Author: Bolu Ajiboye
% Created: October 18, 2005 (ver 1.0)
% Updated: Dec 07, 2006 (ver 2.1)
% Updated: July 21, 2008 (ver 2.3)
%
%--------------------------------------------------------------------------
% Creative Commons License CC BY
%--------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------



defaults = cell(0,3);
defaults(end+1,:) = {'width', 'positive', 1};
defaults(end+1,:) = {'groupnames', 'array', 1:length(barvalues)};
defaults(end+1,:) = {'bw_title', 'string', ''};
defaults(end+1,:) = {'bw_xlabel', 'string', ''};
defaults(end+1,:) = {'bw_ylabel', 'string', ''};
defaults(end+1,:) = {'bw_legend', 'string', ''};
defaults(end+1,:) = {'bw_colormap', 'colormap', jet};
defaults(end+1,:) = {'gridstatus', 'string', 'none'};
defaults(end+1,:) = {'error_sides', 'nonnegative', 2};
defaults(end+1,:) = {'legend_type', 'string', 'plot'};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2
    error('matlabSTORM:invalidArguments', 'Must have at least the first two arguments:  BarWithError(barvalues, errors)');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

width = parameters.width;
groupnames = parameters.groupnames;
bw_title = parameters.bw_title;
bw_xlabel = parameters.bw_xlabel;
bw_ylabel = parameters.bw_ylabel;
bw_colormap = parameters.bw_colormap;
gridstatus = parameters.gridstatus;
bw_legend = parameters.bw_legend;
error_sides = parameters.error_sides;
legend_type = parameters.legend_type;


change_axis = 0;
ymax = 0;
if size(barvalues,1) ~= size(errors,1) || size(barvalues,2) ~= size(errors,2)
	error('barvalues and errors matrix must be of same dimension');
else
	if size(barvalues,2) == 1
		barvalues = barvalues';
		errors = errors';
	end
	if size(barvalues,1) == 1
		barvalues = [barvalues; zeros(1,length(barvalues))];
		errors = [errors; zeros(1,size(barvalues,2))];
		change_axis = 1;
	end
	numgroups = size(barvalues, 1); % number of groups
	numbars = size(barvalues, 2); % number of bars in a group
	if isempty(width)
		width = 1;
	end
	
	% Plot bars
	handles.bars = bar(barvalues, width,'edgecolor','k', 'linewidth', 2);
	hold on
	if ~isempty(bw_colormap)
		colormap(bw_colormap);
	else
		colormap(jet);
	end
	if ~isempty(bw_legend) && ~strcmp(legend_type, 'axis')
		handles.legend = legend(bw_legend, 'location', 'best', 'fontsize',12);
		legend boxoff;
	else
		handles.legend = [];
	end
	
	% Plot erros
	for i = 1:numbars
		x =get(get(handles.bars(i),'children'), 'xdata');
		x = mean(x([1 3],:));
		handles.errors(i) = errorbar(x, barvalues(:,i), errors(:,i), 'k', 'linestyle', 'none', 'linewidth', 2);
		ymax = max([ymax; barvalues(:,i)+errors(:,i)]);
	end
	
	if error_sides == 1
		set(gca,'children', flipud(get(gca,'children')));
	end
	
	ylim([0 ymax*1.1]);
	xlim([0.5 numgroups-change_axis+0.5]);
	
	if strcmp(legend_type, 'axis')
		for i = 1:numbars
			xdata = get(handles.errors(i),'xdata');
			for j = 1:length(xdata)
				text(xdata(j),  -0.03*ymax*1.1, bw_legend(i), 'Rotation', 60, 'fontsize', 12, 'HorizontalAlignment', 'right');
			end
		end
		set(gca,'xaxislocation','top');
	end
	
	if ~isempty(bw_title)
		title(bw_title, 'fontsize',14);
	end
	if ~isempty(bw_xlabel)
		xlabel(bw_xlabel, 'fontsize',14);
	end
	if ~isempty(bw_ylabel)
		ylabel(bw_ylabel, 'fontsize',14);
	end
	
	set(gca, 'xticklabel', groupnames, 'box', 'off', 'ticklength', [0 0], 'fontsize', 12, 'xtick',1:numgroups, 'linewidth', 2,'xgrid','off','ygrid','off');
	if ~isempty(gridstatus) && any(gridstatus == 'x')
		set(gca,'xgrid','on');
	end
	if ~isempty(gridstatus) && any(gridstatus ==  'y')
		set(gca,'ygrid','on');
	end
	
	handles.ax = gca;
	
	hold off
end