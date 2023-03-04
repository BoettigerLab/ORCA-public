 function [spotTable,tforms,fixDist,fsrUsed] = ChromaticAlignFromTable(spotTable,varargin)

% all distance/position values in the input table are now assumed to be 
% in the same units (nanometers). No mixing pixels and nanometers.
% 
% optionally, allows filtering of spots on brightness (ideally this was
% done well the first time). This option also allows explority fits.
% Just specify 'minHeightRef' or 'minHeightOt'.
% 
% ------------------ Notes -------------------------------------
% FOR IMPROVEMENT: we could use some outlier detection -- there is a
% strong prior of how the chromatic shift should look.
% FOR IMPROVEMENT: is polynomial fitting really best? Is there a lower
% dimensional / more accurate neural net representation of chromatic
% distortion?  
% 
% ----------------- Updates ----------------------------------------
% 2020-12-09: revised to use automatic detection of replicate hybs based on
% the labels for readout and channel.  Any readout which is detected in
% both directions will be used. 
% -----------------------------------------------------------------
% Alistair Boettiger
% CC BY Jan 2019

defaults = cell(0,3);
% local parameters
defaults(end+1,:) = {'refChn','integer',647};
defaults(end+1,:) = {'minHeightRef','nonnegative',0};
defaults(end+1,:) = {'minHeightOt','nonnegative',0};
% parameters for Polymap3D
defaults(end+1,:) = {'chnColors','colormap',hsv(3)};
defaults(end+1,:) = {'polyOrder','integer',2};
defaults(end+1,:) = {'showPlots','boolean',true};
defaults(end+1,:) = {'max2D','nonnegative',500}; % in nm 
defaults(end+1,:) = {'max3D','nonnegative',0}; % in nm
pars = ParseVariableArguments(varargin,defaults,mfilename);


% enforce numerical array for channel
chns = spotTable.chn;
if isstring(chns) || ischar(chns)
    spotTable.chn = cellstr(chns);
    chns = spotTable.chn;
end
if iscell(chns)
    spotTable.chn = cellfun(@str2double,chns);
end
datChns = unique(spotTable.chn);
if ~any(datChns == pars.refChn)
    disp('data channels: '); disp(datChns); 
    error(['refChn ', num2str(pars.refChn), ' not found in channel list']);
end
otherChns = datChns(datChns ~= pars.refChn);
nChns = length(otherChns);

% this is where the chromatic shifts will be saved
spotTable.xcShift = zeros(height(spotTable),1);
spotTable.ycShift = zeros(height(spotTable),1);
spotTable.zcShift = zeros(height(spotTable),1);

tforms = cell(nChns,1);
for c=1:nChns % loop over all non-reference colors
    % Match images of the same spot taken in different colors
    otChn = otherChns(c);
    tempTable = spotTable ; % ( comboTable.fov==f,:);
    nReads = max(tempTable.readout);
    % scan for color swap readouts
    xyzRefs = cell(nReads,1);
    xyzAlts = cell(nReads,1);
    for r=1:nReads  % r =4; %  r=6
        r1c1 = tempTable.readout== r & tempTable.chn == pars.refChn;
        r1c2 = tempTable.readout== r & tempTable.chn == otChn;
        if sum(r1c1) > 0 && sum(r1c2) > 0 % we have a color correct chn
           c1Table = tempTable(r1c1,:);
           c2Table = tempTable(r1c2,:);
           [fsrUsed,ia,ir] = intersect(c2Table.fs,c1Table.fs); % match spots by FOV and spot number
           xyzRefs{r} = c1Table{ir,1:3} + [c1Table.locusX(ir),c1Table.locusY(ir),0*ir];
           xyzAlts{r} = c2Table{ia,1:3} + [c2Table.locusX(ia),c2Table.locusY(ia),0*ia];     
           % [tform3D,~,fixDist] = Polymap3D(xyzRefs{r},xyzAlts{r},'figDist',r,'figMap',r+1);
           % pause();
        end
    end
    
    a = cat(1,xyzRefs{:});
    b = cat(1,xyzAlts{:});
    if ~isempty(a) && ~isempty(b)
        [tform3D,~,fixDist] = Polymap3D(a,b,'parameters',pars);
        spotTable = ApplyShiftToTable(spotTable,otChn,tform3D);
        tforms{c} = tform3D; 
    else
        error('no matches found to align with');
    end

    % this is now a function
%     %% ======= Save chromatic shifts to data table =======% 
%     % Add chromatic shifts as xyz shifts for the speicfic points into the table 
%     % (probably best for a new loop over otChns) here, more modular
%     isOt = spotTable.chn == otChn;
%     stageXYZ = [spotTable.x(isOt) + spotTable.locusX(isOt),...
%                 spotTable.y(isOt) + spotTable.locusY(isOt),...
%                 spotTable.z(isOt)];
%     newStageXYZ = tforminv(tform3D,stageXYZ(:,1),stageXYZ(:,2),stageXYZ(:,3));
%     spotTable.xcShift(isOt) = newStageXYZ(:,1) - stageXYZ(:,1);
%     spotTable.ycShift(isOt) = newStageXYZ(:,2) - stageXYZ(:,2);
%     spotTable.zcShift(isOt) = newStageXYZ(:,3) - stageXYZ(:,3);
end 

   