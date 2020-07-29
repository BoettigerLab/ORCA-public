function addedSpots = AddSpotsToFig(figSpots,varargin)
% Use mouse to select pixels to add as data points
% 
% defaults(end+1,:) = {'markSpots','boolean',true};
% defaults(end+1,:) = {'cleanup','boolean',true};
% defaults(end+1,:) = {'markerColor','colormap','r'}; 
% 
% see also SelectSpotsInFig
defaults = cell(0,3);
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'markSpots','boolean',true};
defaults(end+1,:) = {'cleanup','boolean',true};
defaults(end+1,:) = {'markerColor','colormap','r'}; 
pars = ParseVariableArguments(varargin,defaults,mfilename);


% directions
if pars.verbose
    disp('Select spots by left clicking with the cursor. Right click when done.')
end

% set up data cursor mode
clear dcm_obj;
dcm_obj = datacursormode(figSpots);
set(dcm_obj,'DisplayStyle','window','SnapToDataVertex','off','Enable','on');

% record data
addedSpots = zeros(0,2);
marks = cell(0,1);
keepCollecting = true; 
while keepCollecting
    [x,y,button] = myginput(1,'arrow');
    if button == 1
        keepCollecting = true;
        addedSpots(end+1,:) = [x,y]; %#ok<AGROW>
        if pars.markSpots
            figure(figSpots); hold on; 
            marks{end+1} = plot(x,y,'o','color',pars.markerColor); %#ok<AGROW>
        end
    else 
        keepCollecting = false; %#ok<NASGU>
        break;
    end
    
end
set(dcm_obj,'Enable','off');


% clean up the added marks; 
if pars.cleanup
    for i=1:length(marks)
        delete(marks{i});
    end
end
