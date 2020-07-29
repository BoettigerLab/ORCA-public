function [selectedIdx,targetHandle] = SelectSpotsInFig(figSpots,varargin)
% Use mouse to select datapoints in figure
% 
% defaults(end+1,:) = {'markSpots','boolean',true};
% defaults(end+1,:) = {'cleanup','boolean',true};
% defaults(end+1,:) = {'markerColor','colormap','r'}; 
% 
% see also AddSpotsToFig
global figData;

defaults = cell(0,3);
defaults(end+1,:) = {'verbose','boolean',true};
defaults(end+1,:) = {'markSpots','boolean',true};
defaults(end+1,:) = {'cleanup','boolean',true};
defaults(end+1,:) = {'markerColor','colormap','r'}; 
pars = ParseVariableArguments(varargin,defaults,mfilename);


% if spots are passed instead of a figure handle
if ~ishandle(figSpots) && length(figSpots) > 1 
    plot(figSpots(:,1),figSpots(:,2),'k.');
end

if pars.verbose
    disp('Select spots by left clicking with the cursor. Right click when done.')
end

% set up data cursor mode
clear dcm_obj;
dcm_obj = datacursormode(figSpots);
set(dcm_obj,'DisplayStyle','window','SnapToDataVertex','on','Enable','on','UpdateFcn',@RecordData);

% record data
selectedIdx = zeros(0,1);
marks = cell(0,1);
keepCollecting = true; 
while keepCollecting
    [~,~,button] = myginput(1,'arrow');
    if button == 1
        keepCollecting = true;
        if pars.markSpots
            figure(figSpots); hold on; 
            marks{end+1} = plot(figData.pos(1),figData.pos(2),'o','color',pars.markerColor); %#ok<AGROW>
        end
        selectedIdx(end+1,:) = figData.idx; %#ok<AGROW>
    else 
        keepCollecting = false; %#ok<NASGU>
        break;
    end
    
end
set(dcm_obj,'Enable','off');
targetHandle = figData.target;


% clean up the added marks; 
if pars.cleanup
    for i=1:length(marks)
        delete(marks{i});
    end
end

function output_txt = RecordData(~,event_obj)
global figData
% event_obj
figData.pos = get(event_obj,'Position');
figData.idx = get(event_obj,'DataIndex');
figData.target = get(event_obj,'Target');
output_txt = 'f';