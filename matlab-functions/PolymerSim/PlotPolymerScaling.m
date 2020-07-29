function dfit = PlotPolymerScaling(locuslengths,r_g,varargin)
% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'colormap', 'freeType', ''};
defaults(end+1,:) = {'lineWidth', 'positive', 3};
defaults(end+1,:) = {'markerSize', 'positive', 30};
defaults(end+1,:) = {'confInt','array',[]};
defaults(end+1,:) = {'showFit','boolean',true};

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



%%
if ~iscell(locuslengths)
    locuslengths = {locuslengths};
end
if ~iscell(r_g)
    r_g = {r_g};
end

numDatasets = length(r_g);
if isempty(parameters.colormap);
    parameters.colormap = hsv(numDatasets);
end
if ischar(parameters.colormap)
    parameters.colormap = eval([ parameters.colormap, '(',num2str(numDatasets),')'] );
end


%% Original Version
% Plot points
for i=1:numDatasets
    plot(log10(locuslengths{i}),log10(r_g{i}),'.','color',parameters.colormap(i,:),'MarkerSize',parameters.markerSize); hold on;
end

% Add trend lines
if parameters.showFit
    ftype = fittype('c*x+b','coeff',{'c','b'},'ind','x');
    dfit = cell(numDatasets,1); 
    legendEq = cell(numDatasets,1); 
    for i=1:numDatasets
        dfit{i} = fit(log10(locuslengths{i}),log10(r_g{i}),ftype,'StartPoint',[1,0]);
        h(i) = plot(dfit{i}); hold on;
        legendEq{i} = ['y=',num2str(dfit{i}.c,2),'*L+',num2str(dfit{i}.b,2)];
    end
    
     % Recolor trend lines
    for i=1:numDatasets
        set(h(i),'color',parameters.colormap(i,:),'LineWidth',parameters.lineWidth); hold on;
    end
end


% Add legend, axis labels, 
legend(legendEq,'Location','Best');
title('median locus radius of gyration vs length','FontSize',15);
set(gcf,'color','w'); set(gca,'FontSize',13);
xlabel('log_1_0(locus length (kb))','FontSize',15);
ylabel('log_1_0(locus radius (nm))','FontSize',15);
xspan = log10(  cat(1,locuslengths{:}) );
x=linspace(min(xspan),max(xspan),10);
xlim([min(x),max(x)]);