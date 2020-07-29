function [spots,subplotNum,axesHandles] = GetPixelInfo(varargin)

% get coordinates and subplot number from pixel position in select figure

defaults = cell(0,3);
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'markSpots','boolean',true};
defaults(end+1,:) = {'cleanup','boolean',true};
defaults(end+1,:) = {'markerColor','colormap','w'}; 
pars = ParseVariableArguments(varargin,defaults,mfilename);

% prompt user to select locus of interest
if pars.verbose
   disp('Select figure to edit, then press any key'); 
end
pause;
f = gcf;
if pars.verbose
    disp(['Selected figure: ',f.Name]);
end

dcmObj = datacursormode(f);
figure(f); 
spots = zeros(0,2);
axesHandles = cell(0,1);
subplotNum = zeros(0,1);
marks = cell(1,0);

set(dcmObj,'DisplayStyle','window','SnapToDataVertex','off','Enable','on')
if pars.verbose
    disp('Select spots by left clicking with the cursor. Right click when done.')
end
% pause;  % Wait while the user does this by hand
keepCollecting = true; 
while keepCollecting
    [x,y,button] = myginput(1,'arrow');
    if button == 1
        keepCollecting = true;
    else 
        keepCollecting = false; %#ok<NASGU>
        break;
    end
        spots(end+1,:) = [x,y]; %#ok<*AGROW>
        axesHandles{end+1} = gca;
    try
        aNum = find(f.Children == axesHandles{end});
        subplotNum(end+1,1) = length(f.Children)+1 -aNum ; % subplots are number backwards
    catch
        if pars.verbose
            disp('figure has no subplots.');
        end
        subplotNum(end+1,1) = 0;
        axesHandles{end+1} = gca;
    end
    if pars.verbose
        disp(['Position Selected. subplot: ', num2str(subplotNum(end)), ' pixel: ',num2str(spots(end,1)),', ',num2str(spots(end,2))]);
    end
    if pars.markSpots
       subplot(axesHandles{end}); hold on;
       marks{1,end+1} = plot(spots(end,1),spots(end,2),'+','color',pars.markerColor);
    end
end
% clean up the added marks; 
if pars.cleanup
    for i=1:length(marks)
        delete(marks{i});
    end
end